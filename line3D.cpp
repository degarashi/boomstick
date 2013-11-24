#include "geom3D.hpp"

namespace boom {
	namespace geo3d {
		Line::Line(const Vec3& p, const Vec3& d): pos(p), dir(d) {}
		Line Line::FromPoints(const Vec3& p0, const Vec3& p1) {
			return Line(p0, (p1-p0).normalization());
		}
		Line Line::operator * (const AMat43& m) const {
			Line ls;
			ls.pos = pos.asVec4(1) * m;
			ls.dir = dir.asVec4(0) * m;
			return ls;
		}
		Vec3x2 Line::nearestPoint(const Line& s) const {
			auto fn = [](float f) { return f; };
			return NearestPoint(*this, s, fn, fn);
		}
		Vec3 Line::nearest(const Vec3& p) const {
			return NearestPoint(*this, p, [](float f){ return f; });
		}
		float Line::dist_sq(const Vec3& p) const {
			return nearest(p).dist_sq(p);
		}

		bool Line::hit(const Line& ls) const {
			auto fn = [](float f){return f;};
			Vec3x2 res = NearestPoint(*this, ls, fn, fn);
			return res.first.distance(res.second) <= Point::NEAR_THRESHOLD;
		}
	}
}
