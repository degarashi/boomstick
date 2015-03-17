#include "geom3D.hpp"

namespace boom {
	namespace geo3d {
		// --------------------- AABB ---------------------
		AABB::AABB(const Vec3& v_min, const Vec3& v_max): vmin(v_min), vmax(v_max) {}
		AABB AABB::FromPoints(const Vec3* v, int n) {
			Assert(Trap, n>0);
			AABB ab;
			ab.vmin = ab.vmax = v[0];
			for(int i=1 ; i<n ; i++) {
				for(int j=0 ; j<3 ; j++) {
					if(ab.vmin.m[j] > v[1].m[j])
						ab.vmin.m[j] = v[1].m[j];
					else if(ab.vmax.m[j] < v[1].m[j])
						ab.vmax.m[j] = v[1].m[j];
				}
			}
			return ab;
		}

		Mat33 AABB::bs_getInertia() const {
			auto dv = vmax - vmin;
			const float r12 = 1.f/12.f,
						x2 = spn::Square(dv.x),
						y2 = spn::Square(dv.y),
						z2 = spn::Square(dv.z);
			return Mat33(r12*(y2+z2), 0, 0,
						0, r12*(z2+x2), 0,
						0, 0, r12*(x2+y2));
		}
		Vec3 AABB::bs_getGCenter() const {
			return (vmin + vmax) /2;
		}
		Vec3 AABB::bs_getCenter() const {
			return bs_getGCenter();
		}
		Sphere AABB::bs_getBVolume() const {
			auto tmp = vmax-vmin;
			int idx;
			if(tmp.x > tmp.y)
				idx = (tmp.x > tmp.z) ? 0 : 2;
			else
				idx = (tmp.y > tmp.z) ? 1 : 2;
			Sphere sp;
			sp.radius = tmp.m[idx]/2;
			sp.center = vmin + tmp*sp.radius;
			return sp;
		}
		float AABB::bs_getArea() const {
			auto tmp = vmax - vmin;
			return tmp.x * tmp.y * tmp.z;
		}
		Vec3 AABB::support(const Vec3& dir) const {
			return Vec3(dir.x > 0 ? vmax.x : vmin.x,
						dir.y > 0 ? vmax.y : vmin.y,
						dir.z > 0 ? vmax.z : vmin.z);
		}
		AABB AABB::operator * (const AMat43& m) const {
			Vec3 tmp[2] = {vmin.asVec4(1) * m,
							vmax.asVec4(1) * m};
			return FromPoints(tmp, 2);
		}

		bool AABB::hit(const AABB& ab) const {
			for(int i=0 ; i<3 ; i++) {
				if(ab.vmax.m[i] < vmin.m[i] || ab.vmin.m[i] > vmax.m[i])
					return false;
			}
			return true;
		}
		bool AABB::hit(const Plane& plane) const {
			Vec3 v = support(-plane.getNormal());
			return plane.dot(v) <= 0;
		}
		std::ostream& operator << (std::ostream& os, const AABB& a) {
			return os << "AABB(3d) [ min: " << a.vmin << std::endl
						<< "max: " << a.vmax << ']';
		}
	}
}
