#include "geom3D.hpp"

namespace boom {
	namespace geo3d {
		const Vec3& Sphere::bs_getCenter() const { return center; }
		const Vec3& Sphere::bs_getGCenter() const { return center; }
		const Sphere& Sphere::bs_getBSphere() const { return *this; }
		float Sphere::bs_getArea() const { INVOKE_ERROR }
		Mat33 Sphere::bs_getInertia() const { INVOKE_ERROR }
		Vec3 Sphere::support(const Vec3& dir) const {
			return center + dir*radius;
		}
		Sphere Sphere::Cover(MdlItr mI, MdlItr mE) {
			Sphere sp(Vec3(0),0);
			int n = mE - mI;
			if(n == 0)
				return sp;

			auto mI2 = mI;
			while(mI2 != mE)
				sp.center += mI2->im_getCenter();
			sp.center /= n;

			while(mI != mE) {
				Sphere s2 = mI->im_getBSphere();
				float d = sp.center.distance(mI->im_getCenter()) + s2.radius/2;
				if(sp.radius < d)
					sp.radius = d;
			}
			return sp;
		}
		Sphere Sphere::operator * (const AMat43& m) const {
			return Sphere(center.asVec4(1) * m,
						radius);
		}

		bool Sphere::hit(const Sphere& s) const {
			return center.dist_sq(s.center) <= spn::Square(radius + s.radius);
		}
		bool Sphere::hit(const Line& l) const {
			return l.dist_sq(center) <= spn::Square(radius);
		}
		bool Sphere::hit(const Capsule& c) const {
			return c.seg.dist_sq(center) <= spn::Square(radius);
		}
		bool Sphere::hit(const Ray& r) const {
			return r.nearest(center).dist_sq(center) <= spn::Square(radius);
		}
	}
}
