#pragma once
#include "test.hpp"
#include "geom3D.hpp"

namespace boom {
	namespace test3d {
		using spn::RangeF;
		using spn::DegF;
		using spn::RadF;
		using spn::Vec3;
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
		geo3d::FrustumM GenRFrustum(const RDF& rdf,
									const RangeF& r={-1e4f, 1e4f},
									const RangeF& rDist={1e-2f, 1e2f},
									const RangeF& rAngDeg={5.f, 170.f},
									const RangeF& rAsp={1e-1f, 1e1f})
		{
			const auto pos = Vec3::Random(rdf, r);
			auto q = Quat::Random(rdf);
			const auto ang = rdf(rAngDeg);
			const auto dist = rdf(rDist);
			const auto asp = rdf(rAsp);
			return geo3d::FrustumM(pos, Vec3{0,0,1}*q, Vec3{0,1,0}*q, DegF(ang), dist, asp);
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
