#include "geom3D.hpp"

namespace boom {
	namespace geo3d {
		namespace {
			//! 直線と，原点からZ方向へ伸びる無限円錐との交点を求める (内部関数)
			/*! @param[in] rad_d 円錐の長さ1に対する半径の割合
				@param[in] q 直線オフセット
				@param[in] v 直線の傾き
				@return 交点の0番目のtの値, 交点の1番目のtの値, 交点の数 */
			std::tuple<float,float,int> CrossCone(float rad_d, const Vec3& q, const Vec3& v) {
				using spn::Square;
				float a = Square(v.x) + Square(v.y) + Square(rad_d*v.z),
					b = 2*q.x*v.x + 2*q.y*v.y + 2*Square(rad_d)*q.z*v.z,
					c = Square(q.x) + Square(q.y) - Square(rad_d*q.z);

				float D = Square(b)-4*a*c;
				if(std::fabs(D) <= Point::NEAR_THRESHOLD) {
					float val = -b/(2*a);
					return std::make_tuple(val, val, 1);
				}
				if(D <= 0)
					return std::make_tuple(1e8f, 1e8f, 0);
				D = std::sqrt(D);
				float t0 = (-b + D) / (2*a),
					t1 = (-b - D) / (2*a);
				return std::make_tuple(t0, t1, 2);
			}
			//! 直線と, 原点からcdir方向へ伸びる無限円錐との交点を求める (内部関数)
			/*! @param[in] rad_d 円錐の長さ1に対する半径の割合
				@param[in] cdir 円錐の方向
				@param[in] ofs 直線オフセット
				@param[in] dir 直線の方向
				@return 交点の0番目のtの値, 交点の1番目のtの値, 交点の数 */
			std::tuple<float,float,int> _NearConeLine(float rad_d, const Vec3& cdir, const Vec3& ofs, const Vec3& dir) {
				// ofsとdirを円錐座標系へ変換
				Vec3 yA(0,1,0), xA;
				if(std::fabs(cdir.y - 1.0f) <= Point::NEAR_THRESHOLD)
					yA = Vec3(1,0,0);
				xA = yA % cdir;
				xA.normalize();
				yA = cdir % xA;

				Mat33 m(xA.x, yA.x, cdir.x,
						xA.y, yA.y, cdir.y,
						xA.z, yA.z, cdir.z);
				int nVal;
				float f0,f1;
				std::tie(f0,f1,nVal) = CrossCone(rad_d, ofs * m, dir * m);
				if(nVal == 0)
					return std::make_tuple(0.0f,0.0f,0);
				return std::make_tuple(f0,f1,nVal);
			}
			//! 円錐と直線の交差判定
			/*! @param[in] cone 円錐クラス
				@param[in] st 直線クラス
				@return 交点の0番目の座標, 交点の1番目の座標, 交点の数 */
			std::tuple<Vec3,Vec3,int> CrossConeLineV(const Cone& cone, const Line& st) {
				int nf;
				float f0,f1;
				std::tie(f0,f1,nf) = _NearConeLine(cone.radius / cone.length,
									cone.dir,
									st.pos - cone.center,
									st.dir);
				return std::make_tuple(st.pos+st.dir*f0, st.pos+st.dir*f1, nf);
			}
			std::tuple<float,float,int> CrossConeLineF(const Cone& cone, const Line& st) {
				return _NearConeLine(cone.radius / cone.length,
									cone.dir,
									st.pos - cone.center,
									st.dir);
			}
		}

		float Cone::getAngle() const {
			float r = radius / length;
			return std::atan(r);
		}
		bool Cone::hit(const Plane& p) const {
			// 始点と末端の線分が平面と交差していればtrue
			Vec3 vEnd = center + dir*length;
			if(Segment(center, vEnd).crossPoint(p))
				return true;

			// tmp=平面への最短方向
			Vec3 tmp = dir % p.getNormal();
			float len = tmp.length();
			if(std::fabs(len) <= Point::NEAR_THRESHOLD)
				return false;

			tmp = dir % tmp;
			tmp /= len;
			vEnd += tmp * radius;
			return p.dot(vEnd) <= 0;
		}
		bool Cone::hit(const Segment& s) const {
			int nc;
			float f0,f1;
			std::tie(f0,f1,nc) = CrossConeLineF(*this, s.asLine());
			float len = s.getLength();
			return spn::IsInRange(f0, 0.0f, len) || spn::IsInRange(f1, 0.0f, len);
		}
		bool Cone::hit(const Vec3& p) const {
			Vec3 toP = p - center;
			float d = toP.dot(dir);
			if(d > length)
				return false;
			return (center + dir*d).dist_sq(p) <= radius * d / length;
		}
		Vec3 Cone::support(const Vec3& dir) const {
			// 側面の線
			Vec3 endP = center + this->dir * length;
			Vec3 tv = dir % this->dir;
			tv = this->dir % tv;
			tv.normalize();
			tv = endP + tv*radius;

			if(this->dir.dot(dir) >= 1.0f-1e-5f)
				return endP;

			float f0 = center.dot(dir),
				f1 = tv.dot(dir);
			if(f0 > f1)
				return center;
			return tv;
		}
		std::pair<Vec3,float> Cone::nearestPoint(const Plane& plane) const {
			const auto& nml = plane.getNormal();
			Vec3 bnml = nml % dir;
			bnml %= dir;

			Vec3 ret = center + dir*length + bnml*radius;
			return std::make_pair(ret, plane.dot(ret));
		}
		Cone Cone::operator * (const AMat43& m) const {
			Cone c;
			c.center = center.asVec4(1) * m;
			c.dir = dir.asVec4(0) * m;
			c.dir -= c.center;
			c.dir.normalize();
			c.radius = radius;
			c.length = length;
			return c;
		}
	}
}
