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
	}
}
