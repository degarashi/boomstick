#pragma once
#include "test.hpp"
#include "geom3D.hpp"

namespace boom {
	namespace test3d {
		template <class RDF>
		geo3d::PointM GenRPoint(const RDF& rdf, const spn::RangeF& r={-1e4f, 1e4f}) {
			return geo3d::PointM(spn::Vec3::Random(rdf, r));
		}
		template <class RDF>
		geo3d::SphereM GenRSphere(const RDF& rdf, const spn::RangeF& rC={-1e4f, 1e4f},
											const spn::RangeF& rR={0, 1e3f})
		{
			return geo3d::SphereM(spn::Vec3::Random(rdf, rC),
								rdf(rR));
		}
		template <class RDF>
		geo3d::LineM GenRLine(const RDF& rdf, const spn::RangeF& r={-1e4f, 1e4f}) {
			return geo3d::LineM(spn::Vec3::Random(rdf, r),
								spn::Vec3::RandomDir(rdf));
		}
		template <class RDF>
		geo3d::RayM GenRRay(const RDF& rdf, const spn::RangeF& r={-1e4f, 1e4f}) {
			return geo3d::RayM(spn::Vec3::Random(rdf, r),
								spn::Vec3::RandomDir(rdf));
		}
		template <class RDF>
		geo3d::SegmentM GenRSegment(const RDF& rdf, const spn::RangeF& r={-1e4f, 1e4f}) {
			using spn::Vec3;
			return geo3d::SegmentM(Vec3::Random(rdf, r),
									Vec3::Random(rdf, r));
		}
		template <class RDF>
		geo3d::CapsuleM GenRCapsule(const RDF& rdf, const spn::RangeF& rV={-1e4f, 1e4f},
										const spn::RangeF& rR={0, 1e3f})
		{
			return geo3d::CapsuleM(GenRSegment(rdf, rV),
									rdf(rR));
		}
		template <class RDF>
		geo3d::ConeM GenRCone(const RDF& rdf, const spn::RangeF& rV={-1e4f, 1e4f},
									const spn::RangeF& rR={0, 1e3f},
									const spn::RangeF& rL={0, 1e3f})
		{
			return geo3d::ConeM(spn::Vec3::Random(rdf, rV),
								spn::Vec3::RandomDir(rdf),
								rdf(rR),
								rdf(rL));
		}
		template <class RDF>
		geo3d::AABBM GenRAABB(const RDF& rdf, const spn::RangeF& r={-1e4f, 1e4f}) {
			using spn::Vec3;
			Vec3 v0 = Vec3::Random(rdf, r),
				 v1 = Vec3::Random(rdf, r),
				 tmp = v0;
			v0.selectMin(v1);
			v1.selectMax(tmp);
			return geo3d::AABBM(v0, v1);
		}
		template <class RDF>
		geo3d::FrustumM GenRFrustum(const RDF& rdf, const spn::RangeF& rV={-1e4f, 1e4f},
											const spn::RangeF& rFov={0, spn::DEGtoRAD(130.f)},
											const spn::RangeF& rDist={0, 1e3f},
											const spn::RangeF& rAspect={1e-2f, 1e1f})
		{
			using spn::Vec3;
			return geo3d::FrustumM(
						Vec3::Random(rdf, rV),
						Vec3::RandomDir(rdf),
						Vec3::RandomDir(rdf),
						spn::RadF(rdf(rFov)),
						rdf(rDist),
						rdf(rAspect)
					);
		}
		template <class RDF>
		geo3d::ConvexPM GenRConvex(const RDF& rdf, const int n, const spn::RangeF& r={-1e4f, 1e4f}) {
			geo3d::ConvexP cnv(test::GenRVectors<3>(rdf, n, r, NEAR_THRESHOLD_SQ));
			cnv.quickHull();
			return geo3d::ConvexPM(std::move(cnv));
		}
		template <class RDF>
		geo3d::Cylinder GenRCylinder(const RDF& rdf, const spn::RangeF& r={-1e4f, 1e4f}) {
			return geo3d::Cylinder(GenRCapsule(rdf, r));
		}
	}
}
