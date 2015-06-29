#pragma once
#include "test.hpp"
#include "geom2D.hpp"
#include "tfleaf2D.hpp"

namespace boom {
	namespace test2d {
		template <class T>
		using FRand = std::function<T (const spn::Range<T>&)>;
		using FRandF = FRand<float>;
		using FRandI = FRand<int>;
		// Monotoneの内部距離をとっとけば点との判定可能
		class MonoPolygon {
			using DistV = std::vector<Vec2>;
			private:
				spn::Vec2		_vOrigin,
								_vDir;		//!< 軸方向
				// 主軸で左右に分けた場合の頂点数と原点からの距離
				DistV			_distL,
								_distR;
				float			_widthOffset;

				MonoPolygon(DistV&& dL, DistV&& dR, const Vec2& ori, const Vec2& dir, float wofs);
			public:
				geo2d::PointL getPoints() const;
				const spn::Vec2& getDir() const;
				bool hit(const spn::Vec2& p, float threshold=NEAR_THRESHOLD) const;
				static MonoPolygon Random(const FRandF& rff, const FRandI& rfi, const spn::RangeF& rV, const spn::RangeF& rLen, int nMaxV);
		};
		/*! 3角ポリゴンを1つずつ加えていって最後に外周をとればConcaveになる
			3角ポリゴンの集合体なので
			当たり判定のチェックは容易に出来る */
		class ConcavePolygon {
			template <class T>
			struct Neighbor {
				T*	pNeighbor;
				int	edgeId;		//!< pNeighbor側のエッジ
			};
			template <class T>
			using PtrT = std::array<Neighbor<T>, 3>;
			using IdxT = IdxTriangleDataR<PtrT>;
			using PolyPV = std::vector<IdxT*>;
			private:
				geo2d::PointL	_vtx;
				PolyPV			_poly;

				ConcavePolygon() = default;
			public:
				geo2d::PointL getPoints() const;
				bool hit(const spn::Vec2& p, float threshold=NEAR_THRESHOLD) const;
				static ConcavePolygon Random(const FRandF& rff, const FRandI& rfi, const spn::RangeF& rV, int nPoly);
		};

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
		template <class RDF>
		geo2d::PointM GenRPoint(const RDF& rdf, const spn::RangeF& r=defval::point_pos) {
			return geo2d::PointM(spn::Vec2::Random(rdf, r));
		}
		template <class RDF>
		geo2d::CircleM GenRCircle(const RDF& rdf, const spn::RangeF& rC=defval::circle_center,
										const spn::RangeF& rR=defval::circle_radius)
		{
			return geo2d::CircleM(spn::Vec2::Random(rdf, rC),
									rdf(rR));
		}
		template <class RDF>
		geo2d::RayM GenRRay(const RDF& rdf, const spn::RangeF& rV=defval::ray_pos) {
			return geo2d::RayM(spn::Vec2::Random(rdf, rV),
								spn::Vec2::RandomDir(rdf));
		}
		template <class RDF>
		geo2d::LineM GenRLine(const RDF& rdf, const spn::RangeF& rV=defval::line_pos) {
			return GenRRay(rdf, rV).asLine();
		}
		template <class RDF>
		geo2d::SegmentM GenRSegment(const RDF& rdf, const spn::RangeF& rV=defval::segment_pos) {
			auto rv = [&](){ return spn::Vec2::Random(rdf, rV); };
			return geo2d::SegmentM(rv(), rv());
		}
		template <class RDF>
		geo2d::CapsuleM GenRCapsule(const RDF& rdf, const spn::RangeF& rV=defval::capsule_pos,
											const spn::RangeF& rR=defval::capsule_radius) {
			auto rv = [&](){ return spn::Vec2::Random(rdf, rV); };
			return geo2d::CapsuleM(rv(), rv(), rdf(rR));
		}
		template <class RDF>
		geo2d::PolyM GenRPoly(const RDF& rdf, const spn::RangeF& rV=defval::poly_pos) {
			auto rv = [&](){ return spn::Vec2::Random(rdf, rV); };
			auto p = geo2d::PolyM(rv(), rv(), rv());
			// 頂点が時計回りになっているかチェック
			if(!p.isCW())
				p.invert();
			return p;
		}
		template <class RDF>
		geo2d::ConvexM GenRConvex(const RDF& rdf, const spn::RangeF& rV=defval::convex_pos, int n=-1) {
			if(n < 0)
				n = static_cast<int>(rdf({3, 32+1}));

			return geo2d::ConvexM(
					geo2d::Convex::FromConcave(
						test::GenRVectors<2>(rdf, n, rV, NEAR_THRESHOLD_SQ)
					)
			);
		}
		template <class RDF>
		geo2d::AABBM GenRAABB(const RDF& rdf, const spn::RangeF& rV=defval::aabb_pos) {
			auto rv = [&](){ return spn::Vec2::Random(rdf, rV); };
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
				return spn::Iterate::StepIn;
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
		template <class RDF>
		TfSP GenRLeaf(const RDF& rdf, int id) {
			#define MAKELEAF(typ)	case geo2d::typ::GetCID(): return MakeAsLeaf<geo2d::typ>(GenR##typ(rdf));
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
		template <class CTNode, class CTLeaf, class RDF>
		TfSP MakeRandomTree(const RDF& rdf, int nIteration, int maxDepth) {
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
				auto fnI = [&rdf](const spn::RangeI& r){ return static_cast<int>(rdf({r.from, r.to+1})); };
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
							spCursor->addChild(GenRLeaf(rdf, CidId_Leaf[fnI({0,NLeaf-1})]));
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
