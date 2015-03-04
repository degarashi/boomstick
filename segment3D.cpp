#include "geom3D.hpp"

namespace boom {
	namespace geo3d {
		Segment::Segment(const Vec3& p0, const Vec3& p1): from(p0), to(p1) {}
		Vec3 Segment::bs_getGCenter() const {
			return (from + to)/2;
		}
		Vec3 Segment::bs_getCenter() const {
			return bs_getGCenter();
		}
		Sphere Segment::bs_getBVolume() const {
			auto v = bs_getCenter();
			return Sphere(v, (to-from).length()/2);
		}
		const Vec3& Segment::support(const Vec3& dir) const {
			return (to-from).dot(dir) > 0 ? to : from;
		}
		float Segment::dist_sq(const Vec3& p) const {
			return nearest(p).first.dist_sq(p);
		}
		float Segment::getRatio(const Vec3& pos) const {
			Vec3 tv = pos - from;
			float rlen = spn::Rcp22Bit(from.distance(to));
			Vec3 dir = (to-from) * rlen;
			float r = dir.dot(tv);
			return r * rlen;
		}
		LNear Segment::nearest(const Vec3& p) const {
			LinePos lpos;
			float len = from.distance(to);
			Vec3 res = NearestPoint(asLine(), p, [&lpos, len](float f) ->float{
				if(f < 0.0f) {
					lpos = LinePos::Begin;
					return 0;
				}
				if(f > len) {
					lpos = LinePos::End;
					return len;
				}
				lpos = LinePos::OnLine;
				return f;
			});
			return LNear(res, lpos);
		}
		Vec3x2 Segment::nearest(const Segment& l) const {
			auto fn = [](float f){ return spn::Saturate(f, 0.f, 1.f); };
			return NearestPoint(asLine(), l.asLine(), fn, fn);
		}
		float Segment::getLength() const {
			return from.distance(to);
		}
		Vec3 Segment::getDir() const {
			return (to - from).normalization();
		}
		Line Segment::asLine() const {
			return Line(from, getDir());
		}
		Ray Segment::asRay() const {
			return Ray(from, getDir());
		}
		Segment Segment::operator * (const AMat43& m) const {
			Segment s;
			s.from = from.asVec4(1) * m;
			s.to = to.asVec4(1) * m;
			return s;
		}
		std::tuple<bool,float,float> Segment::_crossPoint(const Plane& plane) const {
			float f0 = plane.dot(from),
					f1 = plane.dot(to);
			return std::make_tuple(f0*f1>0, std::fabs(f0), std::fabs(f1));
		}
		spn::Optional<Vec3> Segment::crossPoint(const Plane& plane) const {
			bool b;
			float f0,f1;
			std::tie(b,f0,f1) = _crossPoint(plane);
			if(b)
				return spn::none;
			return getDir() * (f0 / (f0+f1));
		}
		spn::Optional<Vec3> Segment::crossPointFit(const Plane& plane, float threshold) const {
			bool b;
			float f0,f1;
			std::tie(b,f0,f1) = _crossPoint(plane);
			if(b) {
				int flag = 0;
				// 面上に点があるか
				if(f0 <= threshold)
					flag |= 0x01;
				if(f1 <= threshold)
					flag |= 0x02;
				switch(flag) {
					case 0x00: return spn::none;
					case 0x01: return from;
					case 0x02: return to;
					case 0x03: return (from+to)/2;
					default: return spn::none;
				}
			}
			return getDir() * (f0 / (f0+f1));
		}
		bool Segment::hit(const Plane& p) const {
			return crossPointFit(p, 1e-5f);
		}
		std::ostream& operator << (std::ostream& os, const Segment& s) {
			return os << "Segment(3d) [ from: " << s.from << std::endl
						<< "to: " << s.to << ']';
		}
	}
}
