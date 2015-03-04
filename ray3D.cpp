#include "geom3D.hpp"
namespace boom {
	namespace geo3d {
		// --------------------- Ray ---------------------
		Ray::Ray(const Vec3& p, const Vec3& d): pos(p), dir(d) {}
		Ray Ray::FromPoints(const Vec3& from, const Vec3& to) {
			return Ray(from, (to-from).normalization());
		}
		Ray Ray::operator * (const AMat43& m) const {
			return Ray(pos.asVec4(1) * m,
						dir.asVec4(0) * m);
		}
		Vec3 Ray::getPt(float len) const {
			return pos + dir*len;
		}
		Vec3x2 Ray::nearest(const Ray& r) const{
			auto fn = [](float f){ return std::max(0.f,f); };
			return NearestPoint(asLine(), r.asLine(), fn, fn);
		}
		Vec3 Ray::nearest(const Vec3& p) const{
			return NearestPoint(asLine(), p, [](float f){return std::max(0.f,f); });
		}
		const Line& Ray::asLine() const{
			return *reinterpret_cast<const Line*>(this);
		}
		std::ostream& operator << (std::ostream& os, const Ray& r) {
			return os << "Ray(3d) [ pos: " << r.pos << std::endl
						<< "dir: " << r.dir << ']';
		}
	}
}
