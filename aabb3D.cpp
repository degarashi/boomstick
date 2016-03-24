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
		AABB AABB::expand(const float w) const {
			auto proc = [w](auto& f0, auto& f1){
				if(w < 0 && (f1-f0)<std::abs(w)) {
					f0 = f1 = (f1+f0)/2;
				} else {
					f0 -= w;
					f1 += w;
				}
			};
			AABB ret(vmin, vmax);
			proc(ret.vmin.x, ret.vmax.x);
			proc(ret.vmin.y, ret.vmax.y);
			proc(ret.vmin.z, ret.vmax.z);
			return ret;
		}
		bool AABB::hit(const Vec3& p) const {
			using spn::IsInRange;
			return IsInRange(p.x, vmin.x, vmax.x) &&
					IsInRange(p.y, vmin.y, vmax.y) &&
					IsInRange(p.z, vmin.z, vmax.z);
		}
		bool AABB::hit(const AABB& ab) const {
			for(int i=0 ; i<3 ; i++) {
				if(ab.vmax.m[i] < vmin.m[i] || ab.vmin.m[i] > vmax.m[i])
					return false;
			}
			return true;
		}
		const EdgeList AABB::cs_edge = {
			{0,1}, {0,2}, {1,3}, {2,3},
			{4,5}, {4,6}, {5,7}, {6,7},
			{0,4}, {2,6}, {1,5}, {3,7}
		};
		Vec3List AABB::getPoints() const {
			Vec3List tmp(8);
			auto* dst = tmp.data();
			iteratePoints([&dst](const Vec3& v){
				*dst++ = v;
			});
			return tmp;
		}
		bool AABB::isBackside(const Plane& p, const float t) const {
			const Vec3 pt[2] = {vmin, vmax};
			const Vec3 v(
				pt[static_cast<int>(p.a > 0.f)].x,
				pt[static_cast<int>(p.b > 0.f)].y,
				pt[static_cast<int>(p.c > 0.f)].z
			);
			return p.dot(v) < t;
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
