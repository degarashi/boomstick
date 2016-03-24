#include "geom3D.hpp"

namespace boom {
	namespace geo3d {
		void VFrusFloat::init(float nearZ, float farZ, float aspect, spn::RadF fov) {
			nearDistance = nearZ;
			farDistance = farZ;

			float a = aspect,
				e = 1.0f / tan(fov.get()/2),
				f = sqrt(e*e + 1),
				f2 = sqrt(e*e + a*a);
			leftRightX = e / f2;
			leftRightZ = a / f2;

			topBottomY = -e / f;
			topBottomZ = 1.0f / f;
		}
		Frustum::Points::Points(const Points& pts) {
			std::memcpy(this, &pts, sizeof(*this)); }
		Frustum::Points::Points(const Vec3& cen, const Vec3& lt, const Vec3& rt, const Vec3& lb, const Vec3& rb):
			point{cen, lt, rt, lb, rb} {}
		Frustum::Points& Frustum::Points::operator = (const Points& pts) {
			return *(new(this) Points(pts)); }
		Frustum::Points Frustum::Points::operator * (const AMat43& m) const {
			Points pts;
			for(int i=0 ; i<static_cast<int>(countof(point)) ; i++)
				pts.point[i] = point[i].asVec4(1) * m;
			return pts;
		}
		bool Frustum::Points::hit(const Plane& p) const {
			// centerから四隅への線分が交差しているかを調べる
			Segment seg;
			seg.from = point[Center];
			for(int i=1 ; i<NumPos ; i++) {
				seg.to = point[i];
				if(seg.hit(p))
					return true;
			}
			return false;
		}

		Frustum::Planes::Planes(const Planes& ps) {
			std::memcpy(this, &ps, sizeof(Planes)); }
		Frustum::Planes::Planes(const Plane& pL, const Plane& pR, const Plane& pB, const Plane& pT, const Plane& pF, const Plane& pBack):
			plane{pL, pR, pB, pT, pF, pBack} {}
		Frustum::Planes& Frustum::Planes::operator = (const Planes& ps) {
			return *(new(this) Planes(ps)); }
		Frustum::Planes Frustum::Planes::operator * (const AMat43& m) const {
			Planes ret;
			for(int i=0 ; i<static_cast<int>(countof(plane)) ; i++)
				ret.plane[i] = plane[i] * m;
			return ret;
		}
		Frustum::Frustum(const Vec3& ori, const Vec3& dir, const Vec3& up, const spn::RadF fov, const float dist, const float aspect) {
			setOffset(ori);
			setRot(AQuat::FromAxis(up % dir, up, dir));
			const float h = std::tan(fov.get()/2);
			setScale({h*aspect, h, dist});
		}
		Frustum::Frustum(const Vec3& ori, const Vec3& at, const Vec3& up, spn::RadF fov, float aspect):
			Frustum(ori, (at-ori).normalization(), up, fov, (at-ori).length(), aspect)
		{}
		Vec3 Frustum::bs_getGCenter() const {
			auto pts = getPoints(true);
			Vec3 c(0,0,0);
			for(int i=0 ; i<Points::NumPos ; i++)
				c += pts.point[i];
			return c / Points::NumPos;
		}
		Vec3 Frustum::bs_getCenter() const {
			return bs_getGCenter();
		}
		Sphere Frustum::bs_getBVolume() const {
			auto pts = getPoints(true);
			Vec3 c = bs_getCenter();
			float len = c.dist_sq(pts.point[Points::Center]);
			for(int i=1 ; i<Points::NumPos ; i++)
				len = std::max(c.dist_sq(pts.point[i]), len);
			return Sphere(c, spn::Sqrt(len));
		}
		float Frustum::bs_getArea() const { INVOKE_ERROR }
		Mat33 Frustum::bs_getInertia() const { INVOKE_ERROR }

		namespace {
			const float c_len = 1.f / spn::Sqrt(2.f);
		}
		const Frustum::Points Frustum::cs_pointsLocal(
			Vec3(0,0,0),
			Vec3(-1, 1, 1),
			Vec3(1, 1, 1),
			Vec3(-1, -1, 1),
			Vec3(1, -1, 1)
		);
		Frustum::Points Frustum::getPoints(bool bWorld) const {
			if(!bWorld)
				return cs_pointsLocal;
			return cs_pointsLocal * getToWorld();
		}
		Frustum::Planes Frustum::getPlanes(bool bWorld) const {
			const auto& sc = getScale();
			const float lenv = 1.f / Vec2(sc.x, sc.z).length();
			Planes planes(
				// -X: (-sc.x, 0, sc.z) -> (sc.z, 0, sc.x)
				Plane::FromPtDir(Vec3(0), Vec3(sc.z, 0, sc.x)*lenv),
				// +X: (sc.x, 0, sc.z) -> (-sc.z, 0, sc.x)
				Plane::FromPtDir(Vec3(0), Vec3(-sc.z, 0, sc.x)*lenv),
				// -Y: (0, sc.y, sc.z) -> (0, -sc.z, sc.y)
				Plane::FromPtDir(Vec3(0), Vec3(0, -sc.z, sc.y)*lenv),
				// +Y: (0, -sc.y, sc.z) -> (0, sc.z, sc.y)
				Plane::FromPtDir(Vec3(0), Vec3(0, sc.z, sc.y)*lenv),
				Plane::FromPtDir(Vec3(0,0,1), Vec3(0, 0, -1)),
				Plane::FromPtDir(Vec3(0,0,0), Vec3(0,0,1))
			);
			if(!bWorld)
				return planes;
			const auto& tw = getToWorld();
			const auto& rot = getRot();
			for(auto& p : planes.plane) {
				auto ori = p.getOrigin().asVec4(1) * tw;
				auto dir = p.getNormal() * rot;
				p = Plane::FromPtDir(ori, dir);
			}
			return planes;
		}
		//TODO: 遅い実装なので後で何とかする
		Vec3 Frustum::support(const Vec3& dir) const {
			const auto pts = getPoints(true);

			float d = pts.point[Points::Center].dot(dir);
			int idx = 0;
			for(int i=1 ; i<Points::NumPos ; i++) {
				const float td = dir.dot(pts.point[i]);
				if(d < td) {
					d = td;
					idx = i;
				}
			}
			return pts.point[idx];
		}
		Frustum Frustum::operator * (const AMat43& m) const {
			auto ap = spn::AffineParts::Decomp(m);
			return Frustum(Pose3D(getOffset() + ap.offset,
							getRot() >> ap.rotation,
							getScale() * ap.scale));
		}
		Vec3x2 Frustum::getUpRightLocalDir() const {
			const auto& sc = getScale();
			float w = sc.x / sc.z,
				h = sc.y / sc.z;
			float angW = std::atan(w),
				cW = std::cos(angW),
				sW = std::sin(angW);
			float angH = std::atan(h),
				cH = std::cos(angH),
				sH = std::sin(angH);
			return Vec3x2(Vec3(0, -cH, sH), Vec3(-cW, 0, sW));
		}
		//! 片方または両方の点が内部に入っているケースは扱わない
		bool Frustum::crossHit_LocalB(const Points& pts) const {
			Vec3 vU, vR;
			std::tie(vU,vR) = getUpRightLocalDir();
			// Left, Top, Right, Bottom, Back の順
			Plane plane[5] = {Plane(), Plane(vU,0), Plane(vR,0), Plane(), Plane(0,0,-1, getScale().z)};
			vR.x = -vR.x;
			plane[0] = Plane(vR,0);
			vU.y = -vU.y;
			plane[3] = Plane(vU,0);

			// 線分の絶対値化は面倒なので5平面と比較
			for(int i=0 ; i<static_cast<int>(countof(pts.point)) ; i++) {
				auto &p0 = pts.point[i],
					&p1 = pts.point[(i+1)%countof(pts.point)];
				// 線分が必ず2つの平面を横切っている前提 (1つ少なくて済む)
				for(int j=0 ; j<static_cast<int>(countof(plane)-1) ; j++) {
					if(auto cp = Segment(p0,p1).crossPoint(plane[j])) {
						bool bHit = true;
						// cpが前にあるか他の平面とチェック
						for(int k=0 ; k<static_cast<int>(countof(plane)) ; k++) {
							float d = plane[k].dot(*cp);
							if(d < 0) {
								bHit = false;
								break;
							}
						}
						if(bHit) {
							return true;
						}
					}
				}
			}
			return false;
		}
		bool Frustum::inHit_LocalB(const Points& pts) const {
			Vec3 vU, vR;
			std::tie(vU,vR) = getUpRightLocalDir();
			Plane pZ(0,0,-1, getScale().z);

			// check in-hit
			// 効率化のために絶対値の座標で比較する (上, 右, 奥)
			for(auto& p : pts.point) {
				if(p.z >= 0) {
					Vec3 tpt(std::fabs(p.m[0]), std::fabs(p.m[1]), p.m[2]);
					if(vU.dot(tpt)>=0 && vR.dot(tpt)>=0 && pZ.dot(tpt)>=0)
						return true;
				}
			}
			return false;
		}
		int Frustum::checkSide(const Plane& plane, float threshold) const {
			uint32_t flag = 0;
			auto cb = [&flag, threshold](int /*idx*/, float d) {
				if(d >= threshold)
					flag |= 0x02;
				else if(d <= -threshold)
					flag |= 0x01;
			};
			iterateLocalPlane(plane, cb);
			return flag;
		}
		void Frustum::iterateLocalPlane(const Plane& plane, std::function<void (int, float)> cb) const {
			// 平面をローカル空間へ変換
			auto m = getToLocal();
			Vec3 pos(plane.getNormal() * -plane.d),
				dir(plane.getNormal());
			pos = pos.asVec4(1) * m;
			dir = dir.asVec4(0) * m;
			dir.normalize();
			const APlane tplane = APlane::FromPtDir(pos, dir);

			// 5点の符号を比べる
			const auto& pts = getPoints(false);
			const auto& sc = getScale();
			for(int i=0 ; i<static_cast<int>(countof(pts.point)) ; i++) {
				auto& p = pts.point[i];
				cb(i, tplane.dot(Vec3(p.x*sc.x, p.y*sc.y, p.z*sc.z)));
			}
		}
		spn::Rect Frustum::get2DRect(const AMat43& mV, const AMat44& mP) const {
			// カメラZ平面に正対する矩形を算出
			// 四隅と先端のローカル座標
			int nV2 = 0;
			Vec3 v[5] = { Vec3(0,0,0),Vec3(-1,1,1), Vec3(1,1,1), Vec3(1,-1,1), Vec3(-1,-1,1) };
			Vec3 v2[8];
			AMat43 m = getToWorld().convertA44() * mV;
			for(int i=0 ; i<static_cast<int>(countof(v)) ; i++)
				v[i] = v[i].asVec4(1) * m;
			// ビュー座標へ変換し，Z平面0でクリップ
			const int cIdx[8][2] = {
				{0,1},{0,2},{0,3},{0,4},
				{1,2},{2,3},{3,4},{4,1}
			};
			// 裏->表のクリップ座標 + 目標座標をリストに出力
			for(int i=0 ; i<4 ; i++) {
				const auto &v0 = v[cIdx[i][0]],
							&v1 = v[cIdx[i][1]];
				if(v0.z * v1.z < 0) {
					float r = std::fabs(v0.z)/(std::fabs(v0.z)+std::fabs(v1.z));
					Vec3 cp(v0 + (v1-v0)*r);
					cp.z = std::fabs(cp.z);
					v2[nV2++] = cp;
					AssertP(Trap, v2[nV2-1].z >= -1e-3f);
				}
			}
			if(v[0].z > 0)
				v2[nV2++] = v[0];
			for(int i=4 ; i<8 ; i++) {
				const auto &v0 = v[cIdx[i][0]],
							&v1 = v[cIdx[i][1]];
				if(v0.z * v1.z < 0) {
					// クリップ座標を出力
					float r = std::fabs(v0.z)/(std::fabs(v0.z)+std::fabs(v1.z));
					Vec3 cp(v0 + (v1-v0)*r);
					cp.z = std::fabs(cp.z);
					v2[nV2++] = cp;
					AssertP(Trap, v2[nV2-1].z >= -1e-3f)
				}
				if(v1.z > 0) {
					// 目標座標を出力
					v2[nV2++] = v1;
					AssertP(Trap, v2[nV2-1].z >= -1e-3f)
				}
			}
			AssertP(Trap, nV2 <= 8)

			// それぞれプロジェクション変換
			struct MinMax {
				float vmin, vmax;
				MinMax(): vmin(std::numeric_limits<float>::max()),
						vmax(-std::numeric_limits<float>::max()) {}
				void input(float val) {
					vmin = std::min(vmin, val);
					vmax = std::max(vmax, val);
				}
			};
			MinMax mm[2];
			for(int i=0 ; i<nV2 ; i++) {
				Vec4 res = v2[i].asVec4(1) * mP;
				AssertP(Trap, res.w >= -1e-3f)
				res.w = std::max(1e-5f, res.w);
				mm[0].input(res.x/res.w);
				mm[1].input(res.y/res.w);
			}

			return spn::Rect(spn::Saturate(mm[0].vmin, 1.0f),
							 spn::Saturate(mm[0].vmax, 1.0f),
							spn::Saturate(mm[1].vmin, 1.0f),
							spn::Saturate(mm[1].vmax, 1.0f));
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
// 			float w = sc.x / sc.z,
// 				h = sc.y / sc.z;
// 			float angW = std::atan(w),
// 				cW = std::cos(angW),
// 				sW = std::sin(angW);
// 			float angH = std::atan(h),
// 				cH = std::cos(angH),
// 				sH = std::sin(angH);
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
				return tsp.hit(sg0) || tsp.hit(sg1);
			}
		}
		bool Frustum::hit(const Frustum& fr) const {
			// 四角錘をローカル空間に変換する
			auto m = getToLocal();
			auto pts = fr.getPoints(true);
			for(auto& p : pts.point)
				p = p.asVec4(1) * m;

			if(inHit_LocalB(pts) || crossHit_LocalB(pts))
				return true;

			m = fr.getToLocal();
			pts = getPoints(true);
			for(auto& p : pts.point)
				p = p.asVec4(1) * m;
			return fr.inHit_LocalB(pts);
		}
		bool Frustum::hit(const Cone& cone) const {
			// Cone始点がVFrusに入っているか
			if(cone.hit(bs_getGCenter()))
				return true;

			// 視錐台4隅 + 4辺の線分がConeと交差しているか
			auto pts = getPoints(true);
			const Vec3* pPts = pts.point;
			const int pIdx[8][2] = {
				{0,1},{0,2},{0,3},{0,4},
				{1,2},{1,3},{2,4},{3,4}
			};
			for(int i=0 ; i<static_cast<int>(countof(pIdx)) ; i++) {
				if(cone.hit(Segment(pPts[pIdx[i][0]], pPts[pIdx[i][1]])))
					return true;
			}

			auto planes = getPlanes(true);
			// 視錐台5平面に対しての最深点がFrustumに入っているか
			for(auto& pl : planes.plane) {
				if(cone.nearestPoint(pl).second > 0)
					return false;
			}
			return true;
		}
		bool Frustum::hit(const Plane& p) const {
			return getPoints(true).hit(p);
		}
		AABB Frustum::toAABB() const {
			const auto& sc = getScale();
			const auto& rot = getRot();
			const auto& ofs = getOffset();
			const auto xa = rot.getXAxis() * sc.x,
						ya = rot.getYAxis() * sc.y,
						za = rot.getZAxis() * sc.z;
			const AVec3 xA[2] = {xa, -xa},
						yA[2] = {ya, -ya},
						zA[2] = {za, AVec3(0)};
			AABB ret{Vec3(std::numeric_limits<float>::max()),
					Vec3(std::numeric_limits<float>::lowest())};
			for(int i=0 ; i<0b1000 ; i++) {
				const auto p = xA[i&1] +
								yA[(i>>1)&1] +
								zA[(i>>2)&1];
				ret.vmin.selectMin(p);
				ret.vmax.selectMax(p);
			}
			return ret;
		}
		const EdgeList Frustum::Points::cs_edge = {
			{Points::LeftTop, Points::Center},
			{Points::RightTop, Points::Center},
			{Points::LeftBottom, Points::Center},
			{Points::RightBottom, Points::Center},
			{Points::RightTop, Points::LeftTop},
			{Points::RightBottom, Points::LeftBottom},
			{Points::LeftTop, Points::LeftBottom},
			{Points::RightTop, Points::RightBottom}
		};

		bool Frustum::hit(const AABB& a) const {
			const auto pts = getPoints(true);
			Vec3List ptsv(Points::NumPos);
			for(int i=0 ; i<Points::NumPos ; i++)
				ptsv[i] = pts.point[i];
			return HitCheck(ptsv, Points::cs_edge,
							a.getPoints(), AABB::cs_edge, 1e-5f);
		}
		std::ostream& operator << (std::ostream& os, const Frustum& f) {
			return os << "Frustum(3d) [ pose: " << static_cast<const Pose3D&>(f) << ']';
		}
	}
}
