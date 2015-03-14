#pragma once
#include "spinner/tests/test.hpp"
#include "geom2D.hpp"

namespace boom {
	namespace test {
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
			auto rv = [&](){ return spn::test::GenR2Vec(rd, rV); };
			geo2d::PointL pl(n);
			for(int i=0 ; i<n ; i++) {
				Vec2 p;
				bool bLoop;
				do {
					bLoop = false;
					p = rv();
					for(int j=0 ; j<i ; j++) {
						if(pl[j].dist_sq(p) <= geo2d::NEAR_THRESHOLD_SQ) {
							bLoop = true;
							break;
						}
					}
				} while(bLoop);
				pl[i] = p;
			}
			return geo2d::ConvexM(geo2d::Convex::FromConcave(std::move(pl)));
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
