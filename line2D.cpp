#include "geom2D.hpp"

namespace boom {
	namespace geo2d {
		Line::Line(const Vec2& p, const Vec2& d): pos(p), dir(d) {}
		Vec2x2 Line::nearest(const Line& st) const {
			auto fn = [](float f) { return f; };
			return NearestPoint(*this, st, fn, fn);
		}
		Vec2 Line::nearest(const Vec2& p) const {
			return pos + (p-pos).dot(dir);
		}
		float Line::distance(const Vec2& p) const {
			return std::fabs(dir.ccw(p - pos));
		}
		Vec2 Line::placeOnLine(const Vec2& p) const {
			return pos + dir*posDot(p);
		}
		float Line::posDot(const Vec2& p) const {
			Vec2 tv(p-pos);
			return dir.dot(tv);
		}
		Line Line::operator * (const AMat32& m) const {
			return Line{pos.asVec3(1)*m, dir.asVec3(0)*m};
		}
		std::ostream& operator << (std::ostream& os, const Line& l) {
			return os << "Line [ pos: " << l.pos << std::endl
						<< "dir: " << l.dir << ']';
		}
		bool Line::hit(const Vec2& p, float t) const {
			Vec2 v = p - pos;
			float d = std::abs(v.dot(dir));
			return (pos + dir*d).dist_sq(p) <= spn::Square(t);
		}
	}
}
