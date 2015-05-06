#include "geom2D.hpp"

namespace boom {
	namespace geo2d {
		float Point::distance(const Segment& s) const {
			return std::fabs(((*this) - s.from).ccw(s.getDir()));
		}
		LNear Point::nearest(const Segment& s) const {
			return s.nearest(*this);
		}
		bool Point::hit(const Vec2& p, float t) const {
			return distance(p) <= t;
		}
		const Vec2& Point::support(const Vec2& dir) const { return *this; }
		const float& Point::bs_getArea() const { INVOKE_ERROR }
		const Vec2& Point::bs_getCenter() const { return *this; }
		Circle Point::bs_getBVolume() const {
			// 円の半径が0だと点同士の時にヒットしないので微量含める
			return Circle(*this, NEAR_THRESHOLD);
		}
		AABB Point::bs_getBBox() const {
			const Vec2 tmp(NEAR_THRESHOLD);
			return AABB(*this - tmp,
						*this + tmp);
		}
		const float& Point::bs_getInertia() const { INVOKE_ERROR }
		Point Point::operator * (const AMat32& m) const {
			return Point(asVec3(1) * m);
		}
		std::ostream& operator << (std::ostream& os, const Point& c) {
			return os << "Point(2d) [ pos: " << static_cast<const Vec2&>(c) << ']';
		}
	}
}
