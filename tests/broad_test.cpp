#ifdef WIN32
	#include <intrin.h>
#endif
#include "test2D.hpp"
#include "bc_roundrobin.hpp"
#include <unordered_set>
#include "../ntree/ntree.hpp"
#include "../ntree/dim2.hpp"
#include "../ntree/arrayentry.hpp"

namespace boom {
	namespace test {
		using geo2d::PointM;
		using geo2d::Circle;
		using geo2d::CircleM;
		using geo2d::AABB;
		using geo2d::AABBM;
		using geo2d::ModelMgr;
		template <class T>
		class BroadC_Dim2 : public test2d::Narrow {
			public:
				using base_t = test2d::Narrow;
				using ColMgr = typename std::tuple_element<0, T>::type;
				constexpr static int N_AObj = std::tuple_element<1, T>::type::value,
									N_BObj = std::tuple_element<2, T>::type::value;
				constexpr static bool B_CheckCollision = std::tuple_element<3, T>::type::value;
				ModelMgr	_mmgr;
				ColMgr		_cmgr;
			public:
				using HCol = typename ColMgr::SHdl;
				using HLCol = typename ColMgr::LHdl;
				using HLCol_V = std::vector<HLCol>;
				using MMgr = typename ColMgr::MMgr;
				using HLMdl = typename ColMgr::HLMdl;
				BroadC_Dim2(): _cmgr(1.f) {}
				auto& getColMgr() {
					return _cmgr;
				}
		};
		using ColMgr_NTree = ::boom::ColMgr<::boom::ntree::NTree<
												::boom::ntree::CTDim_2D,
												::boom::ntree::CTEnt_Array,
											1>
											, geo2d::Types, uint32_t>;
		using ColMgr_RR = ::boom::ColMgr<BroadC_RoundRobin<Circle>,
											geo2d::Types, uint32_t>;
		template <int N>
		using IConst = std::integral_constant<int, N>;
		using BC_RR_t = std::tuple<ColMgr_RR, IConst<32>, IConst<32>, std::true_type>;
		using BC_NT_t = std::tuple<ColMgr_NTree, IConst<32>, IConst<0>, std::true_type>;
		using BroadCTypeList2D = ::testing::Types<BC_RR_t,
												BC_NT_t>;
		TYPED_TEST_CASE(BroadC_Dim2, BroadCTypeList2D);

		namespace {
			//! 形状をランダムに動かす
			/*! 動かすかどうかもランダム */
			template <class RD>
			void MoveShape(RD& rd, geo2d::TfLeafBase* p) {
				// 50%の確率で動かす = 新たに形状を作成
				if(rd.template getUniform<int>({0,1}) == 0) {
					void* ms = p->getModelSource()->getCore();
					switch(p->getCID()) {
						case CircleM::GetCID(): {
							auto* pt = reinterpret_cast<Circle*>(ms);
							*pt = test2d::GenRCircle(rd);
							break; }
						case AABBM::GetCID(): {
							auto* pt = reinterpret_cast<AABB*>(ms);
							*pt = test2d::GenRAABB(rd);
							break; }
						default:
							Assert(Trap, "unknown shape Id")
					}
					// 動かしたらsetAsChangedを呼ぶ
					p->setAsChanged();
				}
			}
			//! ランダムな階層構造の形状定義
			template <class CT, class RD, class CM>
			auto AddRandomTree(RD& rd, CM& cm, int n, CMask mask, typename CM::user_t ud) {
				std::vector<typename CM::HLCol> v(n);
				for(int i=0 ; i<n ; i++) {
					auto sp = test2d::MakeRandomTree<CT,CT>(rd, 4, 1);
					v[i] = cm.addCol(mask, mgr_model2d.emplace(sp), ud++);
				}
				return std::move(v);
			}
		}
		TYPED_TEST(BroadC_Dim2, CheckCollision) {
			auto rd = this->getRand();
			auto& cm = this->getColMgr();
			auto fnN = [&](const spn::RangeI& r){ return rd.template getUniform<int>(r); };
			using CT = spn::CType<geo2d::Circle, geo2d::AABB>;
			using this_t = BroadC_Dim2<TypeParam>;
			// (NTreeのルーチンがまだ用意出来てない為)
			if(!this_t::B_CheckCollision)
				return;

			// TypeAを適当に追加
			auto vA = AddRandomTree<CT>(rd, cm, fnN({0,this_t::N_AObj}), 0x00000001, 0);
			// TypeBも適当に追加
			auto vB = AddRandomTree<CT>(rd, cm, fnN({0,this_t::N_BObj}), 0x80000001, 1000);

			using HCol = typename this_t::HCol;
			using Narrow_t = typename this_t::Narrow_t;
			using ColMgr = typename this_t::ColMgr;
			using MMgr = typename ColMgr::MMgr;
			using HLMdl = typename ColMgr::HLMdl;
			HLMdl hlMdl = MMgr::_ref().emplace(test2d::MakeRandomTree<CT,CT>(rd, 4, 1));
			// RoundRobinクラスによる判定
			// -> TypeA and TypeBと判定
			using HCSet = std::unordered_set<HCol>;
			HCSet result_rb[2];
			cm.checkCollision(0x00000001, hlMdl, [&](HCol hc){
				result_rb[0].insert(hc);
			});
			// -> TypeAと判定
			cm.checkCollision(0x80000001, hlMdl, [&](HCol hc){
				result_rb[1].insert(hc);
			});

			// 自前で判定
			HCSet result_diy[2];
			for(auto& a : vA) {
				if(Narrow_t::Hit(hlMdl->get(), a->getModel(), 0)) {
					result_diy[0].insert(a);
					result_diy[1].insert(a);
				}
			}
			for(auto& b : vB) {
				if(Narrow_t::Hit(hlMdl->get(), b->getModel(), 0)) {
					result_diy[0].insert(b);
				}
			}

			// 2つの結果を比較
			for(int i=0 ; i<2 ; i++)
				ASSERT_EQ(result_rb[i], result_diy[i]);
		}
		TYPED_TEST(BroadC_Dim2, BroadCollision) {
			auto rd = this->getRand();
			auto& cm = this->getColMgr();
			using this_t = BroadC_Dim2<TypeParam>;
			using HCol = typename this_t::HCol;
			using Narrow_t = typename this_t::Narrow_t;

			auto fnN = [&](const spn::RangeI& r){ return rd.template getUniform<int>(r); };
			using CT = spn::CType<geo2d::Circle, geo2d::AABB>;
			// MSBが0ならTypeA
			// 目印としてUserData=0x00
			auto v0 = AddRandomTree<CT>(rd, cm, fnN({0,this_t::N_AObj}), 0x00000001, 0x0000);
			// MSBが1ならTypeB
			// 目印としてUserData=0x01
			auto v1 = AddRandomTree<CT>(rd, cm, fnN({0,this_t::N_BObj}), 0x80000001, 0x1000);

			auto fnCollect = [](auto& vleaf, auto* mdl) {
				auto v = dynamic_cast<geo2d::TfBase*>(mdl)->shared_from_this();
				auto res = test2d::CollectLeaf(v, (geo2d::TfLeafBase*)nullptr);
				for(auto* r : res)
					vleaf.push_back(r);
			};
			test2d::TfLeaf2DPtr_V	leaf0,
									leaf1;
			for(auto& h : v0)
				fnCollect(leaf0, h->getModel());
			for(auto& h : v1)
				fnCollect(leaf1, h->getModel());

			// RoundRobinを使わずに判定した結果の格納
			auto fnCheck = [](auto& fcmCur, auto& fcmPrev, auto& cm, int id0, int id1, HCol hc0, HCol hc1, auto t){
				auto idc0 = (hc0.getIndex() << 16) | hc1.getIndex();
				auto idc1 = (hc1.getIndex() << 16) | hc0.getIndex();
				auto *p0 = hc0->getModel(),
					 *p1 = hc1->getModel();
				auto idpair = std::make_tuple(id0,id1);
				auto itr = fcmPrev.find(idpair);
				ASSERT_EQ(fcmCur.count(idpair), 0);
				if(Narrow_t::Hit(p0, p1, t)) {
					int count = 0;
					if(itr != fcmPrev.end()) {
						// 衝突カウンタ値の更新
						count = itr->second + 1;
					}
					// エントリ作成
					fcmCur.emplace(idpair, count);
					cm[idc0] = count;
					cm[idc1] = count;
				} else {
					if(itr != fcmPrev.end()) {
						// EndCollision
						cm[idc0] = -1;
						cm[idc1] = -1;
					}
				}
			};
			// チェック用アルゴリズムによる衝突フレーム数のカウント
			using FCMap = std::unordered_map<std::tuple<int,int>, int, spn::TupleHash>;
			FCMap fcmap[2];		// 衝突履歴用に2つ使用
			int fcsw = 0;

			// [MyIndex : OtherIndex]
			using CHMap = std::unordered_map<uint32_t, int>;
			CHMap chmap[2];		// [0]=ColMgr用, [1]=チェック用アルゴリズム
			auto printch = [](CHMap& c) {
				std::cout << "<" << c.size() << " items>" << std::endl;
				for(auto& c2 : c)
					std::cout << std::hex << c2.first << ": " << c2.second << std::endl;
			};
			auto fnCheckR = [](auto& vl, auto& c){
				for(auto& h : vl) {
					uint32_t id = h.get().getIndex() << 16;
					auto hs = h.get();
					h->getCollision([&c, id](auto& hist){
						auto id2 = id | hist.hCol.getIndex();
						ASSERT_EQ(c.count(id2), 0);
						c[id2] = hist.nFrame;
					});
					h->getEndCollision([&c, id](auto& hist){
						auto id2 = id | hist.hCol.getIndex();
						ASSERT_EQ(c.count(id2), 0);
						c[id2] = -1;
					});
				}
			};
			// 形状の変数値をシャッフルしながら何回か比較
			int nShuffle = fnN({5,50});
			int nA = v0.size(),
				nB = v1.size();
			while(nShuffle-- > 0) {
				// RoundRobinによる判定
				cm.update();
				chmap[0].clear();
				fnCheckR(v0, chmap[0]);
				fnCheckR(v1, chmap[0]);

				// 自前判定
				chmap[1].clear();
				auto& fcmPrev = fcmap[fcsw];
				fcsw ^= 1;
				auto& fcmCur = fcmap[fcsw];
				fcmCur.clear();
				auto tm = cm.getAccum();
				for(int i=0 ; i<nA ; i++) {
					// A -> A
					for(int j=i+1 ; j<nA ; j++) {
						fnCheck(fcmCur, fcmPrev, chmap[1], i,j, v0[i], v0[j], tm);
					}
					// A -> B
					for(int j=0 ; j<nB ; j++) {
						fnCheck(fcmCur, fcmPrev, chmap[1], i, j+0x10000, v0[i], v1[j], tm);
					}
				}
				// ---- コンソール出力による確認 ----
				// LogOutput("----- ChMap0 -----");
				// printch(chmap[0]);
				// LogOutput("----- ChMap1 -----");
				// printch(chmap[1]);
				ASSERT_EQ(chmap[0].size(), chmap[1].size());
				for(auto& ch0 : chmap[0]) {
					auto itr = chmap[1].find(ch0.first);
					ASSERT_NE(itr, chmap[1].end());
					ASSERT_EQ(ch0.second, itr->second);
				}

				for(auto* p : leaf0)
					MoveShape(rd, p);
				for(auto* p : leaf1)
					MoveShape(rd, p);
			}
		}
	}
}
