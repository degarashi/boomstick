#ifdef WIN32
	#include <intrin.h>
#endif
#include "test2D.hpp"

namespace boom {
	namespace test {
		using geo2d::PointM;
		class Narrow : public spn::test::RandomTestInitializer {
			protected:
				using base_t = spn::test::RandomTestInitializer;
				using Types = ::boom::geo2d::Types;
				using Narrow_t = Types::Narrow;
				void SetUp() override {
					base_t::SetUp();
					Narrow_t::Initialize();
				}
		};
		namespace {
			using SP = std::shared_ptr<geo2d::TfBase>;
			template <class T, class A>
			SP MakeAsLeaf(const A& s) {
				return std::make_shared<geo2d::TfLeaf<T>>(s);
			}
			template <class T>
			SP MakeAsNode() {
				return std::make_shared<geo2d::TfNode_Static<T>>();
			}
			template <class RD>
			SP MakeRandomTree(RD& rd, int nIteration, int maxDepth) {
				using geo2d::Circle;
				using geo2d::AABB;
				using test2d::GenRCircle;
				using test2d::GenRAABB;
				enum Manipulation {
					MNP_Add,			//!< 現在の階層にリーフノードを加える
					MNP_Up,				//!< 階層を上がる
					MNP_MakeChild,		//!< 子ノードを作ってそこにカーソルを移動
					N_Manipulation
				};
				auto fnI = [&rd](const spn::RangeI& r){ return rd.template getUniform<int>(r); };
				SP spRoot = (fnI({0,1}) == 0) ?
								MakeAsNode<Circle>() : MakeAsNode<AABB>();
				auto spCursor = spRoot;
				int cursorDepth = 0;
				int nIter = fnI({0, nIteration});
				for(int i=0 ; i<nIter ; i++) {
					int m = fnI({0, N_Manipulation-1});
					switch(m) {
						case MNP_Add: {
							spCursor->addChild((fnI({0,1})==0) ?
									MakeAsLeaf<Circle>(GenRCircle(rd)) : MakeAsLeaf<AABB>(GenRAABB(rd)));
							break; }
						case MNP_Up:
							// 深度が0の時は何もしない
							if(cursorDepth > 0) {
								spCursor = spCursor->getParent();
								--cursorDepth;
							}
							break;
						case MNP_MakeChild:
							// 最大深度を超えている時は何もしない
							if(cursorDepth < maxDepth) {
								auto c = (fnI({0,1}) == 0) ?
											MakeAsNode<Circle>() : MakeAsNode<AABB>();
								spCursor->addChild(c);
								spCursor = c;
								++cursorDepth;
							}
							break;
					}
					Assert(Trap, cursorDepth <= maxDepth)
				}
				return spRoot;
			}
		}
		// 階層構造を含めたNarrowPheseテスト (2D)
		TEST_F(Narrow, Dim2) {
			auto rd = getRand();

			// ランダムに当たり判定階層構造を作る
			auto c0 = MakeRandomTree(rd, 64, 8);
			auto c1 = MakeRandomTree(rd, 64, 8);
			std::vector<const geo2d::TfBase*> v0,v1;
			c0->template iterateDepthFirst<false>([&v0](auto& node, int depth){
				if(node.isLeaf())
					v0.push_back(&node);
				return geo2d::TfBase::Iterate::StepIn;
			});
			c1->template iterateDepthFirst<false>([&v1](auto& node, int depth){
				if(node.isLeaf())
					v1.push_back(&node);
				return geo2d::TfBase::Iterate::StepIn;
			});
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
	}
}
