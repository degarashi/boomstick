#pragma once
#include "test.hpp"
#include "geom2D.hpp"
#include "tfleaf2D.hpp"

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
		geo2d::RayM GenRRay(RD& rd, const spn::RangeF& rV={-1e3f, 1e3f}) {
			return geo2d::RayM(spn::test::GenR2Vec(rd, rV),
								spn::test::GenR2Dir(rd));
		}
		template <class RD>
		geo2d::LineM GenRLine(RD& rd, const spn::RangeF& rV={-1e3f, 1e3f}) {
			return GenRRay(rd, rV).asLine();
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
		geo2d::ConvexM GenRConvex(RD& rd, int n=-1, const spn::RangeF& rV={-1e3f, 1e3f}) {
			if(n < 0)
				n = rd.template getUniform<int>({3, 32});

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
		template <class RD>
		void GenRShape(geo2d::PointM& p, RD& rd) { p = GenRPoint(rd); }
		template <class RD>
		void GenRShape(geo2d::CircleM& p, RD& rd) { p = GenRCircle(rd); }
		template <class RD>
		void GenRShape(geo2d::RayM& p, RD& rd) { p = GenRRay(rd); }
		template <class RD>
		void GenRShape(geo2d::LineM& p, RD& rd) { p = GenRLine(rd); }
		template <class RD>
		void GenRShape(geo2d::SegmentM& p, RD& rd) { p = GenRSegment(rd); }
		template <class RD>
		void GenRShape(geo2d::PolyM& p, RD& rd) { p = GenRPoly(rd); }
		template <class RD>
		void GenRShape(geo2d::AABBM& p, RD& rd) { p = GenRAABB(rd); }
		template <class RD>
		void GenRShape(geo2d::ConvexM& p, RD& rd) { p = GenRConvex(rd); }

		class Narrow : public spn::test::RandomTestInitializer {
			protected:
				using base_t = spn::test::RandomTestInitializer;
				using Types = ::boom::geo2d::Types;
				using Narrow_t = Types::Narrow;
				void SetUp() override;
		};

		using TfSP = std::shared_ptr<geo2d::TfBase>;
		using TfSP_V = std::vector<TfSP>;
		using TfBase2DPtr_V = std::vector<geo2d::TfBase*>;
		using TfBase2DPtrC_V = std::vector<const geo2d::TfBase*>;
		template <class T>
		auto CollectLeaf(T& spRoot) {
			std::vector<decltype(spRoot.get())> v;
			spRoot->template iterateDepthFirst<false>([&v](auto& node, int depth){
				if(node.isLeaf())
					v.push_back(&node);
				return geo2d::TfBase::Iterate::StepIn;
			});
			return std::move(v);
		}
		template <class T, class A>
		TfSP MakeAsLeaf(const A& s) {
			return std::make_shared<geo2d::TfLeaf<T>>(s);
		}
		template <class T>
		TfSP MakeAsNode() {
			return std::make_shared<geo2d::TfNode_Static<T>>();
		}
		TfSP GenRNode(int id);
		template <class RD>
		TfSP GenRLeaf(RD& rd, int id) {
			#define MAKELEAF(typ)	case geo2d::typ::GetCID(): return MakeAsLeaf<geo2d::typ>(GenR##typ(rd));
			switch(id) {
				MAKELEAF(Point)
				MAKELEAF(Line)
				MAKELEAF(Ray)
				MAKELEAF(Segment)
				MAKELEAF(AABB)
				MAKELEAF(Poly)
				MAKELEAF(Circle)
				MAKELEAF(Convex)
				default:
					Assert(Trap, false, "unknown collision-id")
			}
			#undef MAKELEAF
			return nullptr;
		}
		template <template <class...> class CT>
		void SetCid(int* dst, CT<>*) {}
		template <template <class...> class CT,
					class T0, class... Ts>
		void SetCid(int* dst, CT<T0, Ts...>*) {
			*dst = T0::GetCID();
			SetCid(dst+1, (CT<Ts...>*)nullptr);
		}
		template <class CTNode, class CTLeaf, class RD>
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
			constexpr int NLeaf = CTLeaf::size,
							NNode = CTNode::size;
			int CidId_Leaf[NLeaf],
				CidId_Node[NNode];
			SetCid(CidId_Leaf, (CTLeaf*)nullptr);
			SetCid(CidId_Node, (CTNode*)nullptr);

			TfSP spRoot = GenRNode(CidId_Node[fnI({0,NNode-1})]);
			auto spCursor = spRoot;
			int cursorDepth = 0;
			int nIter = fnI({0, nIteration});
			for(int i=0 ; i<nIter ; i++) {
				int m = fnI({0, N_Manipulation-1});
				switch(m) {
					case MNP_Add: {
						spCursor->addChild(GenRLeaf(rd, CidId_Leaf[fnI({0,NLeaf-1})]));
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
							auto c = GenRNode(CidId_Node[fnI({0,NNode-1})]);
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
