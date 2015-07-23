#include "geom3D.hpp"
#include "convex.hpp"
#include "spinner/ulps.hpp"

namespace boom {
	namespace geo3d {
		Cylinder::Cylinder(const Capsule& c): Capsule(c) {}
		Cylinder::Cylinder(const Vec3& p0, const Vec3& p1, float r): Capsule(p0, p1, r) {}
		Sphere Cylinder::bs_getBVolume() const{
			Sphere sp;
			// 中心座標計算
			sp.center = bs_getCenter();

			// 軸に対して直交なベクトル -> tv
			Vec3 tv, axis = getDir();
			tv = axis % Vec3(1,0,0);
			tv.normalize();

			// (中心座標から)最遠座標 -> tv
			tv *= radius;
			tv += from;

			// 中心座標から最遠座標までのベクトル -> tv
			tv -= sp.center;
			// tv -> 球の半径
			sp.radius = std::sqrt(tv.dot(tv));
			return sp;
		}
		Vec3 Cylinder::bs_getGCenter() const {
			return Capsule::bs_getGCenter();
		}
		Vec3 Cylinder::bs_getCenter() const {
			return Capsule::bs_getCenter();
		}
		float Cylinder::bs_getArea() const { INVOKE_ERROR }
		Mat33 Cylinder::bs_getInertia() const { INVOKE_ERROR }
		Vec3 Cylinder::support(const Vec3& dir) const {
			Vec3 mydir = getDir();
			Vec3 tdir = mydir % dir;
			tdir %= mydir;
			float len = tdir.length();
			float d = mydir.dot(dir);
			if(len < 1e-5f) {
				if(d >= 0)
					return to;
				return from;
			}
			tdir /= len;
			if(d >= 0)
				return to + tdir * radius;
			return from + tdir * radius;
		}
		Cylinder Cylinder::operator * (const AMat43& m) const {
			Cylinder ret;
			static_cast<Capsule&>(ret) = static_cast<const Capsule&>(*this) * m;
			auto vt = spn::GetVerticalVec(getDir());
			vt = from + vt*radius;
			vt = vt.asVec4(1) * m;
			ret.radius = vt.distance(ret.from);
			return ret;
		}
		Vec3x2 Cylinder::getBeginPlaneV() const {
			// 法線, 平面の一点
			return std::make_pair(getDir(), from);
		}
		Vec3x2 Cylinder::getEndPlaneV() const {
			// 法線, 平面の一点
			return std::make_pair(getDir(), to);
		}
		Plane Cylinder::getBeginPlane() const {
			Vec3 pNor, pPos;
			std::tie(pNor,pPos) = getBeginPlaneV();
			return Plane::FromPtDir(pPos, pNor);
		}
		Plane Cylinder::getEndPlane() const {
			Vec3 pNor, pPos;
			std::tie(pNor,pPos) = getEndPlaneV();
			return Plane::FromPtDir(pPos, pNor);
		}
		AMat43 MatrixYLocal(const Vec3& from, const Vec3& to) {
			AMat43 m;
			Vec3 xaxis, yaxis, zaxis;

			// Y軸の計算
			yaxis = to - from;
			yaxis.normalize();

			// X軸の計算
			zaxis = spn::GetVerticalVec(yaxis);
			AssertP(Trap, spn::EqAbs(zaxis.length(), 1.f, 1e-3f))

			// Z軸の計算
			xaxis = yaxis % zaxis;
			AssertP(Trap, spn::EqAbs(xaxis.length(), 1.f, 1e-3f))

			// X軸セット
			m.setColumn(0, xaxis.asVec4(from.dot(-xaxis)));
			// Y軸セット
			m.setColumn(1, yaxis.asVec4(from.dot(-yaxis)));
			// Z軸セット
			m.setColumn(2, zaxis.asVec4(from.dot(-zaxis)));
			return m;
		}
		void Cylinder::translateLocal(ColCv& pol) const{
			// ポリゴンを円柱ローカルに変換
			AMat43 toCyl = MatrixYLocal(from, to);
			pol *= toCyl;

			// 底面クリップ
 			pol.splitThis(Plane::FromPtDir(Vec3(0,0,0), Vec3(0,1,0)));
			// 上面クリップ
			pol.splitThis(Plane::FromPtDir(Vec3(0,getLength(),0), Vec3(0,-1,0)));
		}
		Vec3 Cylinder::translateLocal(const Vec3& vec) const {
			// 頂点を円柱ローカルに変換
			return vec.asVec4(1) * MatrixYLocal(from, to);
		}
		std::ostream& operator << (std::ostream& os, const Cylinder& c) {
			return os << "Cylinder(3d) [ from: " << c.from << std::endl
						<< "to: " << c.to << std::endl
						<< "radius: " << c.radius << ']';
		}
	}
}
