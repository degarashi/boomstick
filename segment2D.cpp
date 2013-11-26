#include "geom2D.hpp"

namespace boom {
	namespace geo2d {
		Segment::Segment(const Vec2& v0, const Vec2& v1): from(v0), to(v1) {}
		float Segment::distance(const Segment& ls) const {
			auto fn = [](float f) { return spn::Saturate(f, 0.f, 1.f); };
			Vec2x2 vp = NearestPoint(asLine(), ls.asLine(), fn,fn);
			return vp.first.distance(vp.second);
		}
		float Segment::length() const {
			return from.distance(to);
		}
		float Segment::len_sq() const {
			return from.dist_sq(to);
		}
		bool Segment::hit(const Segment& l) const {
			return distance(l) < Point::NEAR_THRESHOLD;
		}
		LNear Segment::nearest(const Vec2& p) const {
			Vec2 toP(p-from),
				toV1(to-from);
			if(toV1.len_sq() <1e-6f)
				return LNear(Vec2(), LinePos::NotHit);
			float lenV1 = toV1.length();
			toV1 *= spn::Rcp22Bit(lenV1);
			float d = toV1.dot(toP);
			if(d <= 0)
				return LNear(from, LinePos::Begin);
			else if(d >= lenV1)
				return LNear(to, LinePos::End);
			else
				return LNear(from+toV1*d, LinePos::OnLine);
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
			toV1.normalize();
			return spn::IsNear(toV1.dot(toP), toP.length(), 1e-5f);
		}
		LNear Segment::crossPoint(const Segment& ls) const {
			auto fn = [](float f){ return f; };
			Vec2x2 v2 = NearestPoint(asLine(), ls.asLine(), fn, fn);
			return LNear(v2.first, online(v2.first) ? LinePos::OnLine : LinePos::NotHit);
		}
		LNear Segment::crossPoint(const Line& l) const {
			float c0 = l.dir.ccw(from-l.pos),
				c1 = l.dir.ccw(to-l.pos);
			if(c0*c1 <= 0) {
				Vec2 diff(to-from);
				c0 = std::fabs(c0);
				float d = c0 / (c0+std::fabs(c1));
				return LNear(from + diff*d, LinePos::OnLine);
			}
			return LNear(Vec2(), LinePos::NotHit);
		}
		Segment Segment::operator * (const AMat32& m) const {
			return Segment{from.asVec3(1)*m, to.asVec3(1)*m};
		}
		Vec2 Segment::bs_getCenter() const {
			return (from + to) * 0.5f;
		}
		Circle Segment::bs_getBCircle() const {
			return Circle(from + to * 0.5f,
							from.distance(to));
		}
		Vec2 Segment::getDir() const {
			return (to-from).normalization();
		}
	}
}
