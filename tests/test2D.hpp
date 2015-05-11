#pragma once
#include "test.hpp"
#include "geom2D.hpp"
#include "tfleaf2D.hpp"

namespace boom {
	namespace test2d {
		using RangeF_OP = spn::Optional<spn::RangeF>;
		using Int_OP = spn::Optional<int>;
		namespace defval {
			const spn::RangeF point_pos{-1e4f, 1e4f},
								circle_center(point_pos),
								circle_radius{0, 1e3f},
								ray_pos{-1e3f, 1e3f},
								line_pos(ray_pos),
								segment_pos(ray_pos),
								capsule_pos(ray_pos),
								capsule_radius(circle_radius),
								poly_pos(ray_pos),
								convex_pos(ray_pos),
								aabb_pos(ray_pos);
		}
		template <class RD>
		geo2d::PointM GenRPoint(RD& rd, const spn::RangeF& r=defval::point_pos) {
			return geo2d::PointM(spn::test::GenR2Vec(rd, r));
		}
		template <class RD>
		geo2d::CircleM GenRCircle(RD& rd, const spn::RangeF& rC=defval::circle_center,
										const spn::RangeF& rR=defval::circle_radius)
		{
			return geo2d::CircleM(spn::test::GenR2Vec(rd, rC),
							rd.template getUniform<float>(rR));
		}
		template <class RD>
		geo2d::RayM GenRRay(RD& rd, const spn::RangeF& rV=defval::ray_pos) {
			return geo2d::RayM(spn::test::GenR2Vec(rd, rV),
								spn::test::GenR2Dir(rd));
		}
		template <class RD>
		geo2d::LineM GenRLine(RD& rd, const spn::RangeF& rV=defval::line_pos) {
			return GenRRay(rd, rV).asLine();
		}
		template <class RD>
		geo2d::SegmentM GenRSegment(RD& rd, const spn::RangeF& rV=defval::segment_pos) {
			auto rv = [&](){ return spn::test::GenR2Vec(rd, rV); };
			return geo2d::SegmentM(rv(), rv());
		}
		template <class RD>
		geo2d::CapsuleM GenRCapsule(RD& rd, const spn::RangeF& rV=defval::capsule_pos,
											const spn::RangeF& rR=defval::capsule_radius) {
			auto rv = [&](){ return spn::test::GenR2Vec(rd, rV); };
			return geo2d::CapsuleM(rv(), rv(),
								rd.template getUniform<float>(rR));
		}
		template <class RD>
		geo2d::PolyM GenRPoly(RD& rd, const spn::RangeF& rV=defval::poly_pos) {
			auto rv = [&](){ return spn::test::GenR2Vec(rd, rV); };
			auto p = geo2d::PolyM(rv(), rv(), rv());
			// 頂点が時計回りになっているかチェック
			if(!p.isCW())
				p.invert();
			return p;
		}
		template <class RD>
		geo2d::ConvexM GenRConvex(RD& rd, const spn::RangeF& rV=defval::convex_pos, int n=-1) {
			if(n < 0)
				n = rd.template getUniform<int>({3, 32});

			return geo2d::ConvexM(
					geo2d::Convex::FromConcave(
						test::GenRVectors<2>(rd, n, rV, NEAR_THRESHOLD_SQ)
					)
			);
		}
		template <class RD>
		geo2d::AABBM GenRAABB(RD& rd, const spn::RangeF& rV=defval::aabb_pos) {
			auto rv = [&](){ return spn::test::GenR2Vec(rd, rV); };
			Vec2 v0 = rv(),
				 v1 = rv(),
				 tmp = v0;
			v0.selectMin(v1);
			v1.selectMax(tmp);
			return geo2d::AABBM(v0, v1);
		}
		template <class... Args>
		void GenRShape(geo2d::PointM& p, Args&&... args) { p = GenRPoint(std::forward<Args>(args)...); }
		template <class... Args>
		void GenRShape(geo2d::CircleM& p, Args&&... args) { p = GenRCircle(std::forward<Args>(args)...); }
		template <class... Args>
		void GenRShape(geo2d::RayM& p, Args&&... args) { p = GenRRay(std::forward<Args>(args)...); }
		template <class... Args>
		void GenRShape(geo2d::LineM& p, Args&&... args) { p = GenRLine(std::forward<Args>(args)...); }
		template <class... Args>
		void GenRShape(geo2d::CapsuleM& p, Args&&... args) { p = GenRCapsule(std::forward<Args>(args)...); }
		template <class... Args>
		void GenRShape(geo2d::SegmentM& p, Args&&... args) { p = GenRSegment(std::forward<Args>(args)...); }
		template <class... Args>
		void GenRShape(geo2d::PolyM& p, Args&&... args) { p = GenRPoly(std::forward<Args>(args)...); }
		template <class... Args>
		void GenRShape(geo2d::AABBM& p, Args&&... args) { p = GenRAABB(std::forward<Args>(args)...); }
		template <class... Args>
		void GenRShape(geo2d::ConvexM& p, Args&&... args) { p = GenRConvex(std::forward<Args>(args)...); }

		class Narrow : public spn::test::RandomTestInitializer {
			protected:
				using base_t = spn::test::RandomTestInitializer;
				using Types = ::boom::geo2d::Types;
				using Narrow_t = Types::Narrow;
				void SetUp() override;
		};

		using TfSP = std::shared_ptr<geo2d::TfBase>;
		using TfSP_V = std::vector<TfSP>;
		using TfLeaf2DPtr_V = std::vector<geo2d::TfLeafBase*>;
		using TfLeaf2DPtrC_V = std::vector<const geo2d::TfLeafBase*>;
		template <class T, class D=decltype(std::declval<T>().get())>
		auto CollectLeaf(T& spRoot, D=nullptr) {
			std::vector<D> v;
			spRoot->template iterateDepthFirst<false>([&v](auto& node, int depth){
				if(node.isLeaf())
					v.push_back(static_cast<D>(&node));
				return geo2d::TfBase::Iterate::StepIn;
			});
			return std::move(v);
		}
		template <class T, class A>
		TfSP MakeAsLeaf(const A& s) {
			return std::make_shared<geo2d::TfLeaf<>>(std::make_shared<A>(s));
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
			Assert(Trap, nIteration>0)
			TfSP spRoot;
			do {
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

				spRoot = GenRNode(CidId_Node[fnI({0,NNode-1})]);
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
			} while(!spRoot->imn_refresh(0));
			return std::move(spRoot);
		}
	}
}
