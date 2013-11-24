#include "geom3D.hpp"

namespace boom {
	namespace geo3d {
		void VFrusFloat::init(float nearZ, float farZ, float aspect, float fov) {
			nearDistance = nearZ;
			farDistance = farZ;

			float a = aspect,
				e = 1.0f / tan(fov/2),
				f = sqrt(e*e + 1),
				f2 = sqrt(e*e + a*a);
			leftRightX = e / f2;
			leftRightZ = a / f2;

			topBottomY = -e / f;
			topBottomZ = 1.0f / f;
		}
		Frustum::Points::Points(const Points& pts) {
			std::memcpy(this, &pts, sizeof(*this)); }
		Frustum::Points& Frustum::Points::operator = (const Points& pts) {
			std::memcpy(this, &pts, sizeof(*this)); return *this; }
		bool Frustum::Points::hit(const Plane& p) const {
			// 視点から四隅への線分が交差しているかを調べる
			Segment seg;
			seg.from = center;
			for(int i=0 ; i<countof(pt) ; i++) {
				seg.to = pt[i];
				if(seg.hit(p))
					return true;
			}
			return false;
		}
		Frustum::Frustum(const Vec3& ori, const Vec3& dir, const Vec3& up, float fov, float dist, float aspect) {
			setOfsVec(ori);
			setRot(AQuat::FromAxis(up % dir, up, dir));
			float h = std::tan(fov/2);
			setScale(h*aspect, h, dist);
		}
		Frustum::Frustum(const Vec3& ori, const Vec3& at, const Vec3& up, float fov, float aspect):
			Frustum(ori, (at-ori).normalization(), up, fov, (at-ori).length(), aspect)
		{}
		Vec3 Frustum::bs_getGCenter() const {
			auto pts = getPoints(true);
			Vec3 c(0,0,0);
			for(int i=0 ; i<5 ; i++)
				c += pts.ar[i];
			return c / 5;
		}
		Vec3 Frustum::bs_getCenter() const {
			return bs_getGCenter();
		}
		Sphere Frustum::bs_getBSphere() const {
			auto pts = getPoints(true);
			Vec3 c = bs_getCenter();
			float len = c.dist_sq(pts.center);
			for(int i=0 ; i<4 ; i++)
				len = std::max(c.dist_sq(pts.pt[i]), len);
			return Sphere(c, spn::Sqrt(len));
		}
		float Frustum::bs_getArea() const { INVOKE_ERROR }
		Mat33 Frustum::bs_getInertia() const { INVOKE_ERROR }

		Frustum::Points Frustum::getPoints(bool bWorld) const {
			auto& sc = getScale();
			float width = sc.x,
					height = sc.y,
					length = sc.z;
			AVec3 vR = getRight() * width/2,
					vU = getUp() * height/2,
					fp = getOffset() + getDir() * length;
			Points pts;
			pts.center = getOffset();
			const float c_coeff[][2] = {{-1,1}, {1,1}, {-1,-1}, {1,-1}};
			for(int i=0 ; i<4 ; i++)
				pts.pt[i] = fp + vR*c_coeff[i][0] + vU*c_coeff[i][1];
			if(bWorld) {
				auto& m = getToWorld();
				for(auto& p : pts.ar)
					p = p.asVec4(1) * m;
			}
			return pts;
		}
		//TODO: 遅い実装なので後で何とかする
		Vec3 Frustum::support(const Vec3& dir) const {
			// dirをローカルに変換
			Vec3 ldir = dir.asVec4(0) * getToLocal();
			auto pts = getPoints(false);

			float d = pts.center.dot(ldir);
			int idx = 0;
			for(int i=0 ; i<4 ; i++) {
				float td = ldir.dot(pts.pt[i]);
				if(d < td) {
					d = td;
					idx = i+1;
				}
			}
			return pts.ar[idx].asVec4(1) * getToWorld();
		}
		Frustum Frustum::operator * (const AMat43& m) const {
			auto ap = spn::DecompAffine(m);
			return Frustum(Pose3D(getOffset() + ap.offset,
							getRot() >> ap.rotation,
							getScale() * ap.scale));
		}

		bool Frustum::hit(const Sphere& sp) const {
			// 球をローカル空間に変換する
			Sphere tsp = sp * getToLocal();
			auto& spc = tsp.center;

			// Frustumは左右対称なのでXの絶対値をとる
			spc.x = std::fabs(spc.x);
			// 同じく上下対称なのでYの絶対値
			spc.y = std::fabs(spc.y);

			const auto& sc = getScale();
			float w = sc.x / sc.z,
				h = sc.y / sc.z;
			float angW = std::atan(w),
				cW = std::cos(angW),
				sW = std::sin(angW);
			float angH = std::atan(h),
				cH = std::cos(angH),
				sH = std::sin(angH);
			Vec3 nml[4] = {
				Vec3(),
				Vec3(),
				Vec3(0, 0, 1),
				Vec3(0, 0, -1)
			};
			std::tie(nml[0],nml[1]) = getUpRightLocalDir();

			int_fast8_t flag = 0;
			for(int i=0 ; i<3 ; i++) {
				float d = nml[i].dot(spc);
				if(d < -tsp.radius)
					return false;
				if(d >= 0)
					flag |= 1<<i;
			}
			// far 面だけは別
			float d = nml[3].dot(spc) + sc.z;
			if(d >= sc.z) {
				// 始点と判定
				return tsp.center.dot(tsp.center) <= spn::Square(tsp.radius);
			} else {
				if(d < -tsp.radius)
					return false;
				if(flag == 0x07)
					return true;

				// 2辺と判定
				Segment sg0(Vec3(0,sc.y,sc.z), Vec3(sc.x,sc.y,sc.z)),
					sg1(Vec3(sc.x,0,sc.z), Vec3(sc.x,sc.y,sc.z));
				return sg0.hit(tsp) || sg1.hit(tsp);
			}
		}
		bool Frustum::hit(const Frustum& fr) const {
			// 四角錘をローカル空間に変換する
			auto m = getToLocal();
			auto pts = fr.getPoints(true);
			for(auto& p : pts.ar)
				p = p.asVec4(1) * m;

			if(inHit_LocalB(pts) || crossHit_LocalB(pts))
				return true;

			m = fr.getToLocal();
			pts = getPoints(true);
			for(auto& p : pts.ar)
				p = p.asVec4(1) * m;
			return fr.inHit_LocalB(pts);
		}
		bool Frustum::hit(const Cone& cone) const {
// 			// Cone始点がVFrusに入っているか
// 			if(cone.hit(bs_getGCenter()))
// 				return true;
//
// 			// 視錐台4隅 + 4辺の線分がConeと交差しているか
// 			auto pts = getPoints(true);
// 			const Vec3* pPts[5] = {
// 				&pts.center//, pts.pt, pts.pt+1, pts.pt+2, pts.pt+3
// 			};
// 			const int pIdx[8][2] = {
// 				{0,1},{0,2},{0,3},{0,4},
// 				{1,2},{1,3},{2,4},{3,4}
// 			};
// 			for(int i=0 ; i<countof(pIdx) ; i++) {
// 				if(Line(*pPts[pIdx[i][0]], *pPts[pIdx[i][1]]).hit(cone))
// 					return true;
// 			}
//
// 			// 視錐台5平面に対しての最深点がFrustumに入っているか
// 			if(cone.nearestPoint(upper) >> mvskip >> mvget > 0) return false;
// 			if(cone.nearestPoint(lower) >> mvskip >> mvget > 0) return false;
// 			if(cone.nearestPoint(left) >> mvskip >> mvget > 0) return false;
// 			if(cone.nearestPoint(right) >> mvskip >> mvget > 0) return false;
// 			if(cone.nearestPoint(back) >> mvskip >> mvget > 0) return false;
//
			return true;
		}
		bool Frustum::hit(const Plane& p) const {
			return getPoints(true).hit(p);
		}
	}
}
