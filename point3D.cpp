#include "geom3D.hpp"

namespace boom {
	namespace geo3d {
		const Vec3& Point::bs_getGCenter() const {
			return *this;
		}
		const Vec3& Point::bs_getCenter() const {
			return *this;
		}
		Sphere Point::bs_getBSphere() const {
			return Sphere(*this, NEAR_THRESHOLD);
		}
		const float& Point::bs_getArea() const {
			INVOKE_ERROR }
		const Mat33& Point::bs_getInertia() const {
			INVOKE_ERROR }
		const Vec3& Point::support(const Vec3& dir) const {
			return *this;
		}

		bool Point::hit(const Point& p) const {
			return distance(p) <= NEAR_THRESHOLD;
		}
		bool Point::hit(const Line& ls) const {
			Vec3 cp = NearestPoint(ls, static_cast<const Vec3&>(*this), [](float f){return f;});
			return distance(cp) <= NEAR_THRESHOLD;
		}
	}
}
