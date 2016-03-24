#include "geom3D.hpp"

namespace boom {
	Vec3 TriangleRatio(const Vec3& v0, const Vec3& v1, const Vec3& v2, const Vec3& vt) {
		const Vec3 	toV1(v1-v0),
					toV2(v2-v0),
					vC(toV1.cross(toV2)),
					toVT(vt-v0);
		const float det = spn::CramerDet(toV1, toV2, vC);
		return spn::CramersRule(toV1, toV2, vC, toVT, spn::Rcp22Bit(det));
	}
	namespace geo3d {
		Crossing CheckCrossing(const Vec3List& vl, const Plane& p, float t) {
			int flag = 0;
			for(auto& v : vl) {
				const auto d = p.dot(v);
				if(std::abs(d) >= t) {
					if(d > 0) flag |= 0b01;
					else flag |= 0b10;
				}
			}
			if(flag == 0)
				return Crossing::OnPlane;
			if(flag == 0b01)
				return Crossing::Front;
			if(flag == 0b10)
				return Crossing::Back;
			return Crossing::Cross;
		}
		bool HitCheck(const Vec3List& vl, const Vec3List& vel, const EdgeList& el, float t) {
			for(auto& e : el) {
				for(auto& v : vl) {
					const auto plane = Plane::FromPts(vel[e.first], vel[e.second], v);
					const auto c0 = CheckCrossing(vl, plane, t),
								c1 = CheckCrossing(vel, plane, t);
					if(c0 != Crossing::Cross &&
						c1 != Crossing::Cross)
					{
						if(c0 == Crossing::OnPlane ||
							c1 == Crossing::OnPlane ||
							c0 != c1)
							return false;
					}
				}
			}
			return true;
		}
		bool HitCheck(const Vec3List& vl0, const EdgeList& el0,
						const Vec3List& vl1, const EdgeList& el1,
						float t)
		{
			if(!HitCheck(vl0, vl1, el1, t))
				return false;
			return HitCheck(vl1, vl0, el0, t);
		}
	}
}
