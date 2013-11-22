//! 3D形状同士の当たり判定に関する実装 (hit関数)
#include "geom3D.hpp"

namespace boom {
	namespace geo3d {
		bool Sphere::hit(const Sphere& s) const {
			return center.dist_sq(s.center) <= spn::Square(radius + s.radius);
		}
		bool Sphere::hit(const Line& l) const {
			return l.dist_sq(center) <= spn::Square(radius);
		}
		bool Sphere::hit(const Capsule& c) const {
			return c.seg.dist_sq(center) <= spn::Square(radius);
		}
		bool Sphere::hit(const Ray& r) const {
			return r.nearest(center).dist_sq(center) <= spn::Square(radius);
		}

		bool Segment::hit(const Capsule& c) const {
			return false;
		}
		bool Frustum::hit(const Frustum& c) const {
			return false;
		}
	}
}
