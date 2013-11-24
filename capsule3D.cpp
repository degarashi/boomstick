#include "geom3D.hpp"

namespace boom {
	namespace geo3d {
		// --------------------- Capsule ---------------------
		Vec3 Capsule::bs_getCenter() const {
			return seg.bs_getCenter();
		}
		Vec3 Capsule::bs_getGCenter() const {
			return seg.bs_getGCenter();
		}
		Sphere Capsule::bs_getBSphere() const {
			auto s = seg.bs_getBSphere();
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
			if((seg.to - seg.from).dot(dir) > 0)
				return Sphere(seg.to, radius).support(dir);
			return Sphere(seg.from, radius).support(dir);
		}
		Capsule Capsule::operator * (const AMat43& m) const {
			Capsule c;
			c.seg = seg * m;
			c.radius = radius;
			return c;
		}
	}
}
