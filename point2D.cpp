#include "geom2D.hpp"

namespace boom {
	namespace geo2d {
		float Point::distance(const Segment& s) const {
			return std::fabs(((*this) - s.from).ccw(s.getDir()));
		}
		LNear Point::nearest(const Segment& s) const {
			return s.nearest(*this);
		}
		bool Point::hit(const Point& p) const {
			return distance(p) <= NEAR_THRESHOLD;
		}
		const Vec2& Point::support(const Vec2& dir) const { return *this; }
		const float& Point::bs_getArea() const { INVOKE_ERROR }
		const Vec2& Point::bs_getCenter() const { return *this; }
		Circle Point::bs_getBVolume() const {
			// 円の半径が0だと点同士の時にヒットしないので微量含める
			return Circle(*this, NEAR_THRESHOLD);
		}
		const float& Point::bs_getInertia() const { INVOKE_ERROR }
	}
}