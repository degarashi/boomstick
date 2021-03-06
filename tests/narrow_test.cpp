#ifdef WIN32
	#include <intrin.h>
#endif
#include "test2D.hpp"

namespace boom {
	namespace test2d {
		void Narrow::SetUp() {
			base_t::SetUp();
			Narrow_t::Initialize();
		}
		TfSP GenRNode(int id) {
			#define MAKENODE(typ)	case typ::GetCID(): return MakeAsNode<typ>();
			switch(id) {
				MAKENODE(geo2d::AABB)
				MAKENODE(geo2d::Circle)
				default:
					Assert(Trap, false, "unknown collision-id")
			}
			#undef MAKENODE
			return nullptr;
		}
	}
	namespace test {
		using geo2d::PointM;
		using Narrow = test2d::Narrow;
		// 階層構造を含めたNarrowPheseテスト (2D)
		TEST_F(Narrow, RandomTree2D) {
			auto rd = getRand();
			auto rdf = rd.template getUniformF<float>();

			// ランダムに当たり判定階層構造を作る
			using CT = spn::CType<geo2d::Circle, geo2d::AABB>;
			auto c0 = test2d::MakeRandomTree<CT,CT>(rdf, 64, 8);
			auto c1 = test2d::MakeRandomTree<CT,CT>(rdf, 64, 8);
			auto v0 = test2d::CollectLeaf(c0),
				 v1 = test2d::CollectLeaf(c1);
			// 個別に判定した結果
			bool b0 = false;
			for(auto* p0 : v0) {
				for(auto* p1 : v1) {
					if((b0 |= Narrow_t::Hit(p0, p1, 0)))
						break;
				}
			}
			// ツリー構造で判定した結果
			bool b1 = Narrow_t::Hit(c0.get(), c1.get(), 0);
			// 両者は一致する筈
			ASSERT_EQ(b0, b1);
		}

		template <class T>
		class TfNode : public Narrow {
			public:
				using Narrow::Narrow;
		};
		using TfNodeTypeList = ::testing::Types<geo2d::Circle,
												geo2d::Segment,
												geo2d::Poly,
												geo2d::Convex>;
		TYPED_TEST_CASE(TfNode, TfNodeTypeList);

		using namespace spn::test;

		template <class T>
		using BVolume2D = TfNode<T>;
		using BVolume2DTL = ::testing::Types<geo2d::Point,
											geo2d::Segment,
											geo2d::AABB,
											geo2d::Poly,
											geo2d::Convex>;
		TYPED_TEST_CASE(BVolume2D, BVolume2DTL);
		// BVolume, BBoxのテスト
		TYPED_TEST(BVolume2D, Test) {
			auto rd = this->getRand();
			auto rdf = rd.template getUniformF<float>();
			using Shape = TypeParam;
			using ShapeM = geo2d::Model<Shape>;

			constexpr spn::RangeF c_range{-1e2f, 1e2f};
			ShapeM s;
			test2d::GenRShape(s, rdf, c_range);
			geo2d::CircleM bv;
			s.im_getBVolume(bv);
			geo2d::AABBM ab;
			s.im_getBVolume(ab);

			// 元形状とBVolume, BBox形状の衝突判定比較
			using this_t = TfNode<TypeParam>;
			constexpr int nCheck = 100;
			for(int i=0 ; i<nCheck ; i++) {
				// 比較対象はPolygon
				auto rp = test2d::GenRPoly(rdf, c_range);

				bool b_base = this_t::Narrow_t::Hit(&s, &rp, 0),
					 b_bv = this_t::Narrow_t::Hit(&bv, &rp, 0),
					 b_ab = this_t::Narrow_t::Hit(&ab, &rp, 0);
				// 元形状が衝突しているならBVolume, BBoxも衝突している
				if(b_base) {
					ASSERT_TRUE(b_bv);
					ASSERT_TRUE(b_ab);
				}
				// BVolumeが衝突していないなら元形状もなし
				if(!b_bv) {
					ASSERT_FALSE(b_base);
				}
				// BBoxが衝突していないなら元形状もなし
				if(!b_ab) {
					ASSERT_FALSE(b_base);
				}
			}
		}
		// 単一ノードによる姿勢変換テスト
		TYPED_TEST(TfNode, TFNode2D) {
			auto rd = this->getRand();
			auto rdf = rd.template getUniformF<float>();
			using Shape = TypeParam;
			using ShapeM = geo2d::Model<Shape>;
			// 基本の形状に姿勢変換を掛けた物(=A)と
			// 変換後の座標で直接生成した物(=B)の2種類を用意
			ShapeM s;
			test2d::GenRShape(s, rdf, spn::RangeF{-1e2f, 1e2f});
			spn::Pose2D ps(Vec2::Random(rdf, {-1e2f, 1e2f}),
							spn::DegF(rdf({-180.f, 180.f})),
							Vec2(rdf({1e-1f, 1e1f})));
			auto spA = std::make_shared<geo2d::TfLeaf<>>(std::make_shared<ShapeM>(s * ps.getToWorld())),
				 spB = std::make_shared<geo2d::TfLeaf<>>(std::make_shared<ShapeM>(s));
			spB->setPose(ps);

			// サポート写像が一致するか確認
			constexpr int nCheck = 100;
			for(int i=0 ; i<nCheck ; i++) {
				auto dir = spn::Vec2::RandomDir(rdf);
				auto p0 = spA->im_support(dir);
				auto p1 = spB->im_support(dir);
				ASSERT_LE(p0.distance(p1), 5e-3f);
			}
			for(int i=0 ; i<nCheck ; i++) {
				// 衝突試験用の形状を1つ用意
				ShapeM check;
				test2d::GenRShape(check, rdf);
				// 結果は同じになる筈
				using this_t = TfNode<TypeParam>;
				bool bA = this_t::Narrow_t::Hit(spA.get(), &check, 0),
					 bB = this_t::Narrow_t::Hit(spB.get(), &check, 0);
				ASSERT_EQ(bA, bB);
			}
		}
		// ツリーノードによる姿勢変換テスト
		TEST_F(Narrow, TfTree2D) {
			auto rd = this->getRand();
			auto rdf = rd.template getUniformF<float>();

			using CTLeaf = spn::CType<geo2d::Circle, geo2d::AABB>;
			using CTNode = spn::CType<geo2d::Circle>;
			auto spA = test2d::MakeRandomTree<CTNode, CTNode>(rdf, 64, 8),
				spB = spA->cloneTree();
			auto vA = test2d::CollectLeaf(spA),
				 vB = test2d::CollectLeaf(spB);

			// 基本の形状に姿勢変換を掛けた物(=A)と
			// 変換後の座標で直接生成した物(=B)の2種類を用意
			spn::Pose2D ps(Vec2::Random(rdf, {-1e3f, 1e3f}),
							spn::DegF(rdf({-180.f, 180.f})),
							Vec2(rdf({1e-2f, 1e2f})));
			for(auto& vs : vA) {
				auto* vsp = static_cast<geo2d::TfLeaf<>*>(vs);
				vsp->setPose(ps);
			}
			auto& toWorld = ps.getToWorld();
			for(auto& vs : vB) {
				auto* vsp = static_cast<geo2d::TfLeaf<>*>(vs);
				auto& mdl = *static_cast<geo2d::Circle*>(vsp->getModelSource()->getCore());
				mdl = mdl * toWorld;
			}

			constexpr int nCheck = 5;
			for(int i=0 ; i<nCheck ; i++) {
				auto t0 = test2d::MakeRandomTree<CTNode, CTLeaf>(rdf, 64, 8);
				bool bA = Narrow_t::Hit(spA.get(), t0.get(), 0),
					bB = Narrow_t::Hit(spB.get(), t0.get(), 0);
				ASSERT_EQ(bA, bB);
			}
		}
	}
}
