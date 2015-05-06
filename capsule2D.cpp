#include "geom2D.hpp"

namespace boom {
	namespace geo2d {
		Capsule::Capsule(const Vec2& f, const Vec2& t, float r):
			from(f), to(t), radius(r)
		{}

		float Capsule::bs_getArea() const { INVOKE_ERROR }
		float Capsule::bs_getInertia() const { INVOKE_ERROR }
		Vec2 Capsule::bs_getCenter() const {
			return (from+to) * 0.5f;
		}
		Circle Capsule::bs_getBVolume() const {
			return Circle(bs_getCenter(),
							from.distance(to) + radius*2);
		}
		AABB Capsule::bs_getBBox() const {
			auto dir = (to-from).normalization();
			return AABB(from - dir*radius,
						to + dir*radius);
		}
		Vec2 Capsule::support(const Vec2& dir) const {
			if(from.dot(dir) > to.dot(dir)) {
				return from + dir*radius;
			}
			return to + dir*radius;
		}
		Capsule Capsule::operator * (const AMat32& m) const {
			Capsule c;
			c.from = from.asVec3(1) * m;
			c.to = to.asVec3(1) * m;
			c.radius = radius * m.convertA22().calcDeterminant();
			return c;
		}
		Capsule& Capsule::operator += (const Vec2& ofs) {
			from += ofs;
			to += ofs;
			return *this;
		}
		bool Capsule::hit(const Vec2& p, float t) const {
			return Segment(from, to).nearest(p).first.distance(p) <= radius+t;
		}
		bool Capsule::hit(const Segment& s, float t) const {
			Segment s0(from, to);
			return s0.distance(s) <= radius+t;
		}
		bool Capsule::hit(const Capsule& c, float t) const {
			Segment s0(from, to),
					s1(c.from, c.to);
			return s0.distance(s1) <= radius+c.radius+t;
		}
		std::ostream& operator << (std::ostream& os, const Capsule& c) {
			return os << "Capsule(2d) [ from: " << c.from << std::endl
						<< "to: " << c.to << std::endl
						<< "radius: " << c.radius << ']';
		}
	}
}
