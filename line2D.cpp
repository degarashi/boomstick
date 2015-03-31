#include "geom2D.hpp"

namespace boom {
	namespace geo2d {
		Line::Line(const Vec2& p, const Vec2& d): pos(p), dir(d) {}
		Vec2x2 Line::nearest(const Line& st) const {
			auto fn = [](float f) { return f; };
			return NearestPoint(*this, st, fn, fn);
		}
		Vec2 Line::nearest(const Vec2& p) const {
			return pos + dir * (p-pos).dot(dir);
		}
		float Line::bs_getArea() const { INVOKE_ERROR }
		float Line::bs_getInertia() const { INVOKE_ERROR }
		const Vec2& Line::bs_getCenter() const {
			return pos;
		}
		const Circle& Line::bs_getBVolume() const { INVOKE_ERROR }

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
		LineDivision Line::checkSide(const Vec2& p, float t) const {
			auto p2 = p - pos;
			float d = dir.cw(p2);
			if(d < -t)
				return LineDivision::Ccw;
			if(d > t)
				return LineDivision::Cw;
			return LineDivision::OnLine;
		}
		Line Line::operator * (const AMat32& m) const {
			return Line{pos.asVec3(1)*m, dir.asVec3(0)*m};
		}
		std::ostream& operator << (std::ostream& os, const Line& l) {
			return os << "Line [ pos: " << l.pos << std::endl
						<< "dir: " << l.dir << ']';
		}
		Vec2 Line::support(const Vec2& dir) const {
			if(dir.dot(dir) > 0)
				return pos + dir*INFINITY_LENGTH;
			return pos + dir*-INFINITY_LENGTH;
		}
		bool Line::hit(const Vec2& p, float t) const {
			Vec2 v = p - pos;
			float d = std::abs(v.dot(dir));
			return (pos + dir*d).dist_sq(p) <= spn::Square(t);
		}
	}
}
