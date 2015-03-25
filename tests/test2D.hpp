#pragma once
#include "test.hpp"
#include "geom2D.hpp"

namespace boom {
	namespace test2d {
		template <class RD>
		geo2d::PointM GenRPoint(RD& rd, const spn::RangeF& r={-1e4f, 1e4f}) {
			return geo2d::PointM(spn::test::GenR2Vec(rd, r));
		}
		template <class RD>
		geo2d::CircleM GenRCircle(RD& rd, const spn::RangeF& rC={-1e4f, 1e4f},
											const spn::RangeF& rR={0, 1e3f})
		{
			return geo2d::CircleM(spn::test::GenR2Vec(rd, rC),
							rd.template getUniform<float>(rR));
		}
		template <class RD>
		geo2d::SegmentM GenRSegment(RD& rd, const spn::RangeF& rV={-1e3f, 1e3f}) {
			auto rv = [&](){ return spn::test::GenR2Vec(rd, rV); };
			return geo2d::SegmentM(rv(), rv());
		}
		template <class RD>
		geo2d::PolyM GenRPoly(RD& rd, const spn::RangeF& rV={-1e3f, 1e3f}) {
			auto rv = [&](){ return spn::test::GenR2Vec(rd, rV); };
			return geo2d::PolyM(rv(), rv(), rv());
		}
		template <class RD>
		geo2d::ConvexM GenRConvex(RD& rd, int n, const spn::RangeF& rV={-1e3f, 1e3f}) {
			return geo2d::ConvexM(
					geo2d::Convex::FromConcave(
						test::GenRVectors<2>(rd, n, rV, NEAR_THRESHOLD_SQ)
					)
			);
		}
		template <class RD>
		geo2d::AABBM GenRAABB(RD& rd, const spn::RangeF& rV={-1e3f, 1e3f}) {
			auto rv = [&](){ return spn::test::GenR2Vec(rd, rV); };
			Vec2 v0 = rv(),
				 v1 = rv(),
				 tmp = v0;
			v0.selectMin(v1);
			v1.selectMax(tmp);
			return geo2d::AABBM(v0, v1);
		}
		class Narrow : public spn::test::RandomTestInitializer {
			protected:
				using base_t = spn::test::RandomTestInitializer;
				using Types = ::boom::geo2d::Types;
				using Narrow_t = Types::Narrow;
				void SetUp() override;
		};

		using TfSP = std::shared_ptr<geo2d::TfBase>;
		using TfSP_V = std::vector<TfSP>;
		using TfBase2DPtr_V = std::vector<const geo2d::TfBase*>;
		TfBase2DPtr_V CollectLeaf(const TfSP& spRoot);
		template <class T, class A>
		TfSP MakeAsLeaf(const A& s) {
			return std::make_shared<geo2d::TfLeaf<T>>(s);
		}
		template <class T>
		TfSP MakeAsNode() {
			return std::make_shared<geo2d::TfNode_Static<T>>();
		}
		template <class RD>
		TfSP MakeRandomTree(RD& rd, int nIteration, int maxDepth) {
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
			TfSP spRoot = (fnI({0,1}) == 0) ?
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
}
