#include "geom3D.hpp"

namespace boom {
	namespace geo3d {
		// --------------------- Capsule ---------------------
		Capsule::Capsule(const Vec3& p0, const Vec3& p1, float r):
			Segment(p0, p1), radius(r) {}
		Capsule::Capsule(const Segment& s, float r):
			Segment(s), radius(r) {}
		Vec3 Capsule::bs_getCenter() const {
			return bs_getCenter();
		}
		Vec3 Capsule::bs_getGCenter() const {
			return bs_getGCenter();
		}
		Sphere Capsule::bs_getBVolume() const {
			auto s = bs_getBVolume();
			s.radius += radius;
			return s;
		}
		float Capsule::bs_getArea() const {
			INVOKE_ERROR
		}
		Mat33 Capsule::bs_getInertia() const {
			INVOKE_ERROR
		}
		Vec3 Capsule::support(const Vec3& dir) const {
			if((to - from).dot(dir) > 0)
				return Sphere(to, radius).support(dir);
			return Sphere(from, radius).support(dir);
		}
		Capsule Capsule::operator * (const AMat43& m) const {
			Capsule c;
			static_cast<Segment&>(c) = static_cast<Segment&>(c) * m;
			c.radius = radius;
			return c;
		}

		bool Capsule::hit(const Vec3& p) const {
			Vec3 cp = NearestPoint(asLine(), p, [](float f){return spn::Saturate(f,0.f,1.f);});
			return cp.dist_sq(p) <= spn::Square(radius);
		}
		namespace {
			const auto fnSat = [](float f){ return spn::Saturate(f, 0.f,1.f); };
		}
		bool Capsule::hit(const Segment& s) const {
			Vec3x2 res = NearestPoint(asLine(), s.asLine(), fnSat, fnSat);
			return res.first.dist_sq(res.second) <= spn::Square(radius);
		}
		bool Capsule::hit(const Capsule& c) const {
			Vec3x2 res = NearestPoint(asLine(), c.asLine(), fnSat, fnSat);
			return res.first.dist_sq(res.second) <= spn::Square(radius + c.radius);
		}
		std::ostream& operator << (std::ostream& os, const Capsule& c) {
			return os << "Capsule(3d) [ from: " << c.from << std::endl
						<< "to: " << c.to << std::endl
						<< "radius: " << c.radius << ']';
		}
	}
}
