//! 3D形状自体の実装 (hit関数以外)
#include "geom3D.hpp"

namespace boom {
	namespace geo3d {
		// --------------------- IModel ---------------------
		void IModel::notifyChange() {
			if(pParent)
				pParent->notifyChange();
		}
		void IModel::applyChange() {}
		MdlIP IModel::getInner() const { return MdlIP(); }
		void* IModel::getUserData() const { return nullptr; }

		// --------------------- Point ---------------------
		const Vec3& Point::bs_getGCenter() const { return *this; }
		const Vec3& Point::bs_getCenter() const { return *this; }
		Sphere Point::bs_getBSphere() const { return Sphere(*this, NEAR_THRESHOLD); }
		const float& Point::bs_getArea() const { INVOKE_ERROR }
		const Mat33& Point::bs_getInertia() const { INVOKE_ERROR }
		const Vec3& Point::support(const Vec3& dir) const { return *this; }

		// --------------------- Line ---------------------
		Line::Line(const Vec3& p, const Vec3& d): pos(p), dir(d) {}
		Line Line::FromPoints(const Vec3& p0, const Vec3& p1) {
			return Line(p0, (p1-p0).normalization());
		}
		Line Line::operator * (const AMat43& m) const {
			Line ls;
			ls.pos = pos.asVec4(1) * m;
			ls.dir = dir.asVec4(0) * m;
			return ls;
		}
		Vec3x2 Line::nearestPoint(const Line& s) const {
			auto fn = [](float f) { return f; };
			return NearestPoint(*this, s, fn, fn);
		}
		Vec3 Line::nearest(const Vec3& p) const {
			return NearestPoint(*this, p, [](float f){ return f; });
		}
		float Line::dist_sq(const Vec3& p) const {
			return nearest(p).dist_sq(p);
		}

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

		// --------------------- Segment ---------------------
		Segment::Segment(const Vec3& p0, const Vec3& p1): from(p0), to(p1) {}
		Vec3 Segment::bs_getGCenter() const {
			return (from + to)/2;
		}
		Vec3 Segment::bs_getCenter() const {
			return bs_getGCenter();
		}
		Sphere Segment::bs_getBSphere() const {
			auto v = bs_getCenter();
			return Sphere(v, (to-from).length()/2);
		}
		const Vec3& Segment::support(const Vec3& dir) const {
			return (to-from).dot(dir) > 0 ? to : from;
		}
		float Segment::dist_sq(const Vec3& p) const {
			return nearest(p).first.dist_sq(p);
		}
		float Segment::getRatio(const Vec3& pos) const {
			Vec3 tv = pos - from;
			float rlen = spn::Rcp22Bit(from.distance(to));
			Vec3 dir = (to-from) * rlen;
			float r = dir.dot(tv);
			return r * rlen;
		}
		LNear Segment::nearest(const Vec3& p) const {
			LINEPOS lpos;
			float len = from.distance(to);
			Vec3 res = NearestPoint(asLine(), p, [&lpos, len](float f) ->float{
				if(f < 0.0f) {
					lpos = LINEPOS::Begin;
					return 0;
				}
				if(f > len) {
					lpos = LINEPOS::End;
					return len;
				}
				lpos = LINEPOS::OnLine;
				return f;
			});
			return LNear(res, lpos);
		}
		Vec3x2 Segment::nearest(const Segment& l) const {
			auto fn = [](float f){ return spn::Saturate(f, 0.f, 1.f); };
			return NearestPoint(asLine(), l.asLine(), fn, fn);
		}
		float Segment::getLength() const {
			return from.distance(to);
		}
		Vec3 Segment::getDir() const {
			return (to - from).normalization();
		}
		Line Segment::asLine() const {
			return Line(from, getDir());
		}
		Ray Segment::asRay() const {
			return Ray(from, getDir());
		}
		Segment Segment::operator * (const AMat43& m) const {
			Segment s;
			s.from = from.asVec4(1) * m;
			s.to = to.asVec4(1) * m;
			return s;
		}
		std::tuple<bool,float,float> Segment::_crossPoint(const Plane& plane) const {
			float f0 = plane.dot(from),
					f1 = plane.dot(to);
			return std::make_tuple(f0*f1>0, std::fabs(f0), std::fabs(f1));
		}
		spn::Optional<Vec3> Segment::crossPoint(const Plane& plane) const {
			bool b;
			float f0,f1;
			std::tie(b,f0,f1) = _crossPoint(plane);
			if(b)
				return spn::none;
			return getDir() * (f0 / (f0+f1));
		}
		spn::Optional<Vec3> Segment::crossPointFit(const Plane& plane, float threshold) const {
			bool b;
			float f0,f1;
			std::tie(b,f0,f1) = _crossPoint(plane);
			if(b) {
				int flag = 0;
				// 面上に点があるか
				if(f0 <= threshold)
					flag |= 0x01;
				if(f1 <= threshold)
					flag |= 0x02;
				switch(flag) {
					case 0x00: return spn::none;
					case 0x01: return from;
					case 0x02: return to;
					case 0x03: return (from+to)/2;
					default: return spn::none;
				}
			}
			return getDir() * (f0 / (f0+f1));
		}

		// --------------------- Sphere ---------------------
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

		// --------------------- Capsule ---------------------
		Vec3 Capsule::bs_getCenter() const {
			return seg.bs_getCenter();
		}
		Vec3 Capsule::bs_getGCenter() const {
			return seg.bs_getGCenter();
		}
		Sphere Capsule::bs_getBSphere() const {
			auto s = seg.bs_getBSphere();
			s.radius += radius;
			return s;
		}
		float Capsule::bs_getArea() const {
			INVOKE_ERROR
		}
		Mat33 Capsule::bs_getInertia() const {
			INVOKE_ERROR
		}
		Vec3 Capsule::support(const Vec3& dir) const {
			if((seg.to - seg.from).dot(dir) > 0)
				return Sphere(seg.to, radius).support(dir);
			return Sphere(seg.from, radius).support(dir);
		}
		Capsule Capsule::operator * (const AMat43& m) const {
			Capsule c;
			c.seg = seg * m;
			c.radius = radius;
			return c;
		}

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

		Vec3 AABB::bs_getGCenter() const {
			return (vmin + vmax) /2;
		}
		Vec3 AABB::bs_getCenter() const {
			return bs_getGCenter();
		}
		Sphere AABB::bs_getBSphere() const {
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

		// --------------------- Frustum ---------------------
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
// 		bool Frustum::hit(const Cone& cone) const {
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
// 			return true;
// 		}
		bool Frustum::hit(const Plane& p) const {
			return getPoints(true).hit(p);
		}
		Frustum Frustum::operator * (const AMat43& m) const {
			auto ap = spn::DecompAffine(m);
			return Frustum(Pose3D(getOffset() + ap.offset,
							getRot() >> ap.rotation,
							getScale() * ap.scale));
		}
	}
}
