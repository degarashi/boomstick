#pragma once
#include "test.hpp"
#include "geom3D.hpp"

namespace boom {
	namespace test3d {
		template <class RD>
		geo3d::PointM GenRPoint(RD& rd, const spn::RangeF& r={-1e4f, 1e4f}) {
			return geo3d::PointM(spn::test::GenR3Vec(rd, r));
		}
		template <class RD>
		geo3d::SphereM GenRSphere(RD& rd, const spn::RangeF& rC={-1e4f, 1e4f},
											const spn::RangeF& rR={0, 1e3f})
		{
			return geo3d::SphereM(spn::test::GenR3Vec(rd, rC),
								rd.template getUniform<float>(rR));
		}
		template <class RD>
		geo3d::LineM GenRLine(RD& rd, const spn::RangeF& r={-1e4f, 1e4f}) {
			return geo3d::LineM(spn::test::GenR3Vec(rd, r),
								spn::test::GenR3Dir(rd));
		}
		template <class RD>
		geo3d::RayM GenRRay(RD& rd, const spn::RangeF& r={-1e4f, 1e4f}) {
			return geo3d::RayM(spn::test::GenR3Vec(rd, r),
								spn::test::GenR3Dir(rd));
		}
		template <class RD>
		geo3d::SegmentM GenRSegment(RD& rd, const spn::RangeF& r={-1e4f, 1e4f}) {
			using spn::test::GenR3Vec;
			return geo3d::SegmentM(GenR3Vec(rd, r),
									GenR3Vec(rd, r));
		}
		template <class RD>
		geo3d::CapsuleM GenRCapsule(RD& rd, const spn::RangeF& rV={-1e4f, 1e4f},
										const spn::RangeF& rR={0, 1e3f})
		{
			return geo3d::CapsuleM(GenRSegment(rd, rV),
									rd.template getUniform<float>(rR));
		}
		template <class RD>
		geo3d::ConeM GenRCone(RD& rd, const spn::RangeF& rV={-1e4f, 1e4f},
									const spn::RangeF& rR={0, 1e3f},
									const spn::RangeF& rL={0, 1e3f})
		{
			return geo3d::ConeM(spn::test::GenR3Vec(rd, rV),
								spn::test::GenR3Dir(rd),
								rd.template getUniform<float>(rR),
								rd.template getUniform<float>(rL));
		}
		template <class RD>
		geo3d::AABBM GenRAABB(RD& rd, const spn::RangeF& r={-1e4f, 1e4f}) {
			using spn::test::GenR3Vec;
			Vec3 v0 = GenR3Vec(rd, r),
				 v1 = GenR3Vec(rd, r),
				 tmp = v0;
			v0.selectMin(v1);
			v1.selectMax(tmp);
			return geo3d::AABBM(v0, v1);
		}
		template <class RD>
		geo3d::FrustumM GenRFrustum(RD& rd, const spn::RangeF& rV={-1e4f, 1e4f},
											const spn::RangeF& rFov={0, spn::DEGtoRAD(130.f)},
											const spn::RangeF& rDist={0, 1e3f},
											const spn::RangeF& rAspect={1e-2f, 1e1f})
		{
			auto rf = [&](auto& r){ return rd.template getUniform<float>(r); };
			using spn::test::GenR3Vec;
			using spn::test::GenR3Dir;
			return geo3d::FrustumM(
						GenR3Vec(rd, rV),
						GenR3Dir(rd),
						GenR3Dir(rd),
						spn::RadF(rf(rFov)),
						rf(rDist),
						rf(rAspect)
					);
		}
		template <class RD>
		geo3d::ConvexPM GenRConvex(RD& rd, const int n, const spn::RangeF& r={-1e4f, 1e4f}) {
			geo3d::ConvexP cnv(test::GenRVectors<3>(rd, n, r, NEAR_THRESHOLD_SQ));
			cnv.quickHull();
			return geo3d::ConvexPM(std::move(cnv));
		}
		template <class RD>
		geo3d::Cylinder GenRCylinder(RD& rd, const spn::RangeF& r={-1e4f, 1e4f}) {
			return geo3d::Cylinder(GenRCapsule(rd, r));
		}
	}
}
