#include "geom2D.hpp"

namespace boom {
	namespace geo2d {
		Segment::Segment(const Vec2& v0, const Vec2& v1): from(v0), to(v1) {}
		float Segment::distance(const Segment& ls) const {
			if(hit(ls))
				return 0;

			auto ln = nearest(ls.from);
			float distsq = ls.from.dist_sq(ln.first);
			ln = nearest(ls.to);
			distsq = std::min(distsq, ls.to.dist_sq(ln.first));
			ln = ls.nearest(from);
			distsq = std::min(distsq, from.dist_sq(ln.first));
			ln = ls.nearest(to);
			distsq = std::min(distsq, to.dist_sq(ln.first));
			return spn::Sqrt(distsq);
		}
		float Segment::length() const {
			return from.distance(to);
		}
		float Segment::len_sq() const {
			return from.dist_sq(to);
		}
		bool Segment::hit(const Vec2& p, float t) const {
			Vec2 v = p - from,
				 dir = to-from;
			float len = dir.normalize();
			if(len < ZEROVEC_LENGTH)
				return Point(from).hit(p);

			float d = v.dot(dir);
			if(d < -DOT_THRESHOLD)
				return false;
			if(d > len+DOT_THRESHOLD)
				return false;
			v -= dir*d;
			return v.len_sq() <= t*SQUARE_RATIO;
		}
		bool Segment::hit(const Segment& l, float t) const {
			return IsCrossing(asLine(), l.asLine(), length(), l.length(), t);
		}
		LNear Segment::nearest(const Vec2& p) const {
			Vec2 toP(p-from),
				toV1(to-from);
			if(toV1.len_sq() < ZEROVEC_LENGTH_SQ)
				return LNear(from, LinePos::Begin);
			float lenV1 = toV1.normalize();
			float d = toV1.dot(toP);
			if(d <= 0)
				return LNear(from, LinePos::Begin);
			else if(d >= lenV1)
				return LNear(to, LinePos::End);
			else {
				float t = d / lenV1;
				return LNear(from*(1-t) + to*t, LinePos::OnLine);
			}
		}
		float Segment::ratio(const Vec2& p) const {
			Vec2 toP(p-from),
				toV1(to-from);
			float len = toV1.length();
			return toV1.dot(toP) * spn::Rcp22Bit(len);
		}
		Line Segment::asLine() const {
			return Line(from, (to-from).normalization());
		}
		Vec2 Segment::support(const Vec2& dir) const {
			float d[2] = {dir.dot(from), dir.dot(to)};
			if(d[0] > d[1])
				return from;
			return to;
		}
		bool Segment::online(const Vec2& p) const {
			Vec2 toV1(to-from),
				toP(p-from);
			float len = toV1.normalize();
			float d = toV1.dot(toP);
			return spn::IsInRange(d, 0.f, len+1e-5f);
		}
		Vec2_OP Segment::crossPoint(const Segment& s) const {
			if(auto cp = crossPoint(s.asLine())) {
				if(s.online(*cp))
					return cp;
			}
			if(auto cp = s.crossPoint(asLine())) {
				if(online(*cp))
					return cp;
			}
			return spn::none;
		}
		Vec2_OP Segment::crossPoint(const Line& l) const {
			float c0 = l.dir.dot(from - l.pos),
					c1 = l.dir.dot(to - l.pos);
			if(c0*c1 <= 0) {
				Vec2 diff(to-from);
				c0 = std::abs(c0);
				c1 = std::abs(c1);
				float d = c0 / (c0+c1);
				return from + diff*d;
			}
			return spn::none;
		}
		Segment Segment::operator * (const AMat32& m) const {
			return Segment{from.asVec3(1)*m, to.asVec3(1)*m};
		}
		Vec2 Segment::bs_getCenter() const {
			return (from + to) * 0.5f;
		}
		AABB Segment::bs_getBBox() const {
			AABB ab(from, from);
			ab.minV.selectMin(to);
			ab.maxV.selectMax(to);
			return ab;
		}
		Circle Segment::bs_getBVolume() const {
			return Circle((from + to) * 0.5f,
							from.distance(to));
		}
		Vec2 Segment::getDir() const {
			return (to-from).normalization();
		}
		std::ostream& operator << (std::ostream& os, const Segment& s) {
			return os << "Segment(2d) [ from: " << s.from << std::endl
						<< "to: " << s.to << ']';
		}
	}
}
