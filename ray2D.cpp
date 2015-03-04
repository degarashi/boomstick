#include "geom2D.hpp"

namespace boom {
	namespace geo2d {
		Ray::Ray(const Vec2& p, const Vec2& d): pos(p), dir(d) {}
		const Line& Ray::asLine() const { return *reinterpret_cast<const Line*>(this); }
		Vec2x2 Ray::nearest(const Ray& r) const {
			auto fn = [](float f) { return std::max(0.f, f); };
			return NearestPoint(asLine(), r.asLine(), fn, fn);
		}
		Vec2 Ray::nearest(const Vec2& p) const {
			return NearestPoint(asLine(), p, [](float f){ return std::max(0.f,f); });
		}
		Ray Ray::operator * (const AMat32& m) const {
			return Ray{pos.asVec3(1)*m, dir.asVec3(0)*m};
		}
		std::ostream& operator << (std::ostream& os, const Ray& r) {
			return os << "Ray [ pos: " << r.pos << std::endl
						<< "dir: " << r.dir << ']';
		}
	}
}
