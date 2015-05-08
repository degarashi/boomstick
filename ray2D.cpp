#include "geom2D.hpp"

namespace boom {
	namespace geo2d {
		Ray::Ray(const Vec2& p, const Vec2& d): pos(p), dir(d) {}
		const Line& Ray::asLine() const { return *reinterpret_cast<const Line*>(this); }
		Vec2x2 Ray::nearest(const Ray& r) const {
			auto fn = [](float f) { return std::max(0.f, f); };
			if(auto v2 = NearestPoint(asLine(), r.asLine(), fn, fn))
				return *v2;
			auto p = nearest(r.pos);
			return {p, r.pos};
		}
		Vec2 Ray::nearest(const Vec2& p) const {
			return NearestPoint(asLine(), p, [](float f){ return std::max(0.f,f); });
		}
		Ray Ray::operator * (const AMat32& m) const {
			return Ray{pos.asVec3(1)*m, dir.asVec3(0)*m};
		}
		float Ray::bs_getArea() const { INVOKE_ERROR }
		float Ray::bs_getInertia() const { INVOKE_ERROR }
		const Vec2& Ray::bs_getCenter() const {
			return pos;
		}
		const Circle& Ray::bs_getBCircle() const { INVOKE_ERROR }
		const AABB& Ray::bs_getBBox() const { INVOKE_ERROR }
		std::ostream& operator << (std::ostream& os, const Ray& r) {
			return os << "Ray [ pos: " << r.pos << std::endl
						<< "dir: " << r.dir << ']';
		}
		Vec2 Ray::support(const Vec2& dir) const {
			if(dir.dot(dir) > 0)
				return pos + dir*INFINITY_LENGTH;
			return pos;
		}
		bool Ray::hit(const Vec2& p, float t) const {
			Vec2 v = p - pos;
			float d = v.dot(dir);
			if(d < 0)
				return false;
			return (pos + dir*d).dist_sq(p) <= t*SQUARE_RATIO;
		}
	}
}
