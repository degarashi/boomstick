#include "geom3D.hpp"

namespace boom {
	namespace geo3d {
		const Vec3& Sphere::bs_getCenter() const { return center; }
		const Vec3& Sphere::bs_getGCenter() const { return center; }
		const Sphere& Sphere::bs_getBVolume() const { return *this; }
		float Sphere::bs_getArea() const { INVOKE_ERROR }
		Mat33 Sphere::bs_getInertia() const { INVOKE_ERROR }
		Vec3 Sphere::support(const Vec3& dir) const {
			return center + dir*radius;
		}
		void Sphere::setBoundary(const IModel* m) {
			*this = m->im_getBVolume();
		}
		void Sphere::appendBoundary(const IModel* m) {
			Sphere s2 = m->im_getBVolume();
			Vec3 toS2 = s2.center - center;
			float lensq = toS2.len_sq(),
				  rdsq = spn::Square(radius  - s2.radius);
			// 現在の球で収まるか？
			if(lensq > rdsq) {
				toS2 *= spn::RSqrt(lensq);
				// 新たな球を算出
				Vec3 tv0 = s2.center + toS2 * s2.radius,
					 tv1 = center - toS2 * radius;
				center = (tv0 + tv1) * .5f;
				radius = tv0.distance(tv1) * .5;
			}
		}
		Sphere Sphere::operator * (const AMat43& m) const {
			return Sphere(center.asVec4(1) * m,
						radius);
		}

		bool Sphere::hit(const Vec3& p) const {
			return center.dist_sq(p) <= spn::Square(radius);
		}
		bool Sphere::hit(const Sphere& s) const {
			return center.dist_sq(s.center) <= spn::Square(radius + s.radius);
		}
		bool Sphere::hit(const Line& l) const {
			return l.dist_sq(center) <= spn::Square(radius);
		}
		bool Sphere::hit(const Capsule& c) const {
			return c.dist_sq(center) <= spn::Square(radius);
		}
		bool Sphere::hit(const Ray& r) const {
			return r.nearest(center).dist_sq(center) <= spn::Square(radius);
		}
		bool Sphere::hit(const Segment& s) const {
			Vec3 cp = NearestPoint(s.asLine(), center, [](float f){return f;});
			return cp.dist_sq(center) <= spn::Square(radius);
		}
		std::ostream& operator << (std::ostream& os, const Sphere& s) {
			return os << "Sphere(3d) [ center: " << s.center << std::endl
					<< "radius: " << s.radius << ']';
		}
	}
}
