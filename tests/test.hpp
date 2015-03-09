#pragma once
#include "spinner/tests/test.hpp"
#include "geom2D.hpp"

namespace boom {
	namespace test {
		template <class RD>
		geo2d::PointM GenRPoint(RD& rd, const spn::RangeF& r={-1e4f, 1e4f}) {
			auto rc = [&](){ return rd.template getUniform<float>(r); };
			return geo2d::PointM(spn::test::GenRVec<2,false>(rc));
		}
		template <class RD>
		geo2d::CircleM GenRCircle(RD& rd, const spn::RangeF& rC={-1e4f, 1e4f},
											const spn::RangeF& rR={0, 1e3f})
		{
			auto rc = [&](){ return rd.template getUniform<float>(rC); };
			return geo2d::CircleM(spn::test::GenRVec<2,false>(rc),
							rd.template getUniform<float>(rR));
		}
		template <class RD>
		geo2d::SegmentM GenRSegment(RD& rd, const spn::RangeF& rV={-1e3f, 1e3f}) {
			auto rc = [&](){ return rd.template getUniform<float>(rV); };
			auto rv = [&](){ return spn::test::GenRVec<2,false>(rc); };
			return geo2d::SegmentM(rv(), rv());
		}
		template <class RD>
		geo2d::PolyM GenRPoly(RD& rd, const spn::RangeF& rV={-1e3f, 1e3f}) {
			auto rc = [&](){ return rd.template getUniform<float>(rV); };
			auto rv = [&](){ return spn::test::GenRVec<2,false>(rc); };
			return geo2d::PolyM(rv(), rv(), rv());
		}
		template <class RD>
		geo2d::ConvexM GenRConvex(RD& rd, int n, const spn::RangeF& rV={-1e3f, 1e3f}) {
			auto rc = [&](){ return rd.template getUniform<float>(rV); };
			auto rv = [&](){ return spn::test::GenRVec<2,false>(rc); };
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
			auto rc = [&](){ return rd.template getUniform<float>(rV); };
			auto rv = [&](){ return spn::test::GenRVec<2,false>(rc); };
			Vec2 v0 = rv(),
				 v1 = rv(),
				 tmp = v0;
			v0.selectMin(v1);
			v1.selectMax(tmp);
			return geo2d::AABBM(v0, v1);
		}
		template <class RD>
		spn::Vec2 GenR2Vec(RD& rd, const spn::RangeF& rV={-1e4f, 1e4f}) {
			auto rc = [&](){ return rd.template getUniform<float>(rV); };
			return spn::test::GenRVec<2,false>(rc);
		}
		template <class RD>
		spn::Vec2 GenR2Dir(RD& rd) {
			auto rc = [&](){ return rd.template getUniform<float>({-1.f, 1.f}); };
			return spn::test::GenRDir<2,false>(rc);
		}
	}
}
