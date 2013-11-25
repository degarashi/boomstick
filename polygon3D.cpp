#include "convex.hpp"

namespace boom {
	namespace geo3d {
		Polygon::Polygon() {}
		Polygon::Polygon(const Vec3* vsrc) {
			init(vsrc);
		}
		Polygon::Polygon(const Vec3& v0, const Vec3& v1, const Vec3& v2) {
			init(v0, v1, v2);
		}
		void Polygon::setUserData(uint32_t dat) {
			_ud = dat;
		}
		uint32_t Polygon::getUserData() const {
			return _ud;
		}
		void Polygon::init(const Vec3& v0, const Vec3& v1, const Vec3& v2) {
			_vtx[0] = v0;
			_vtx[1] = v1;
			_vtx[2] = v2;
			spn::Bit::Set(_rflg, NORMALFLAG | PLANEFLAG | CENTERFLAG);
		}
		void Polygon::init(const Vec3* vsrc) {
			_vtx[0] = vsrc[0];
			_vtx[1] = vsrc[1];
			_vtx[2] = vsrc[2];
			spn::Bit::Set(_rflg, NORMALFLAG | PLANEFLAG | CENTERFLAG);
		}
		const Vec3& Polygon::getNormal() const {
			// NORMAL更新フラグチェック
			if(spn::Bit::ChClear(_rflg, NORMALFLAG)) {
				// 法線更新
				_vNormal = *spn::NormalFromPoints(_vtx[0], _vtx[1], _vtx[2]);
			}
			return _vNormal;
		}
		const Vec3& Polygon::getVtx(int n) const {
			return _vtx[n];
		}
		Polygon Polygon::operator * (const AMat43& m) const {
			Polygon ret;
			// 3頂点を行列変換
			for(int i=0 ; i<3 ; i++)
				ret._vtx[i] = _vtx[i].asVec4(1) * m;
			spn::Bit::Set(ret._rflg, NORMALFLAG | PLANEFLAG | CENTERFLAG);
			// ユーザーデータコピー
			ret.setUserData(getUserData());
			return ret;
		}
		void Polygon::setVtx(int n, const Vec3& src) {
			_vtx[n] = src;
			spn::Bit::Set(_rflg, NORMALFLAG | PLANEFLAG | CENTERFLAG);
		}
		const Plane& Polygon::getPlane() const {
			// PLANE更新フラグチェック
			if(spn::Bit::ChClear(_rflg, PLANEFLAG)) {
				// 平面生成
				_plane = Plane::FromPtDir(_vtx[0], getNormal());
			}
			return _plane;
		}

		const Vec3& Polygon::getGCenter3D() const {
			if(spn::Bit::ChClear(_rflg, CENTERFLAG)) {
				// 中心点計算
				_vCenter *= 0;
				for(int i=0 ; i<3 ; i++)
					_vCenter += _vtx[i];
				_vCenter /= 3;
			}
			return _vCenter;
		}

		float Polygon::calcRadius() const{
			getGCenter3D();
			float dist = 0;
			for(int i=0 ; i<3 ; i++) {
				float td = _vtx[i].distance(_vCenter);
				if(td > dist)
					dist = td;
			}
			return dist;
		}
		float Polygon::calcAcreage() const {
			Vec3 v0 = _vtx[1] - _vtx[0],
				v1 = _vtx[2] - _vtx[0];
			if(v0.len_sq() < std::numeric_limits<float>::epsilon() ||
				v1.len_sq() < std::numeric_limits<float>::epsilon())
				return 0.0f;
			return (v0 % v1).length() / 2;
		}

		Segment Polygon::getEdge(int n) const {
			return Segment(_vtx[n], _vtx[(n+1)%3]);
		}

		Vec3 Polygon::_calcLineHit(const Vec3& pos, const Vec3& dir) const {
			Vec3 toV1 = _vtx[1]-_vtx[0],
				toV2 = _vtx[2]-_vtx[0];
			float det = spn::CramerDet(toV1, toV2, dir);
			return spn::CramersRule(toV1, toV2, dir, pos-_vtx[0], 1.0f/det);
		}
		bool Polygon::_IsValidRange(const Vec3& res) {
			return spn::IsInRange(res.x, 0.0f, 1.0f) &&
				spn::IsInRange(res.y, 0.0f, 1.0f) &&
				spn::IsInRange(res.x+res.y, 0.0f, 1.0f);
		}
		bool Polygon::hit(const Ray& ray) const {
			Vec3 res = _calcLineHit(ray.pos, ray.dir);
			if(_IsValidRange(res)) {
				if(res.z < 0)
					return std::fabs(getPlane().dot(ray.pos)) >= 0;
				return true;
			}
			return false;
		}
		bool Polygon::hit(const Line& ls) const {
			// Rayの始点チェックを省いたバージョン
			Vec3 res = _calcLineHit(ls.pos, ls.dir);
			return _IsValidRange(res);
		}
		std::pair<Vec3,bool> Polygon::nearest(const Vec3& p) const {
			const auto& nml = getNormal();
			Vec3 res = _calcLineHit(p, nml);
			Vec3 cp = _vtx[0] + (_vtx[1]-_vtx[0]) * res.x +
						(_vtx[2]-_vtx[0]) * res.y;
			if(_IsValidRange(res))
				return std::make_pair(cp, true);

			// 各線分との距離を計算
			Vec3 minNP = Segment(_vtx[0], _vtx[1]).nearest(p).first;
			float minD = minNP.dist_sq(p);
			for(int i=1 ; i<3 ; i++) {
				Vec3 np = Segment(_vtx[i], _vtx[(i+1)%3]).nearest(p).first;
				float d = np.dist_sq(p);
				if(d < minD) {
					minD = d;
					minNP = np;
				}
			}
			return std::make_pair(minNP, false);
		}
		Vec3 Polygon::support(const Vec3& dir) const {
			float d[3];
			for(int i=0 ; i<3 ; i++)
				d[i] = dir.dot(_vtx[i]);
			if(d[0] > d[1]) {
				if(d[0] > d[2])
					return _vtx[0];
				return _vtx[2];
			} else {
				if(d[1] > d[2])
					return _vtx[1];
				return _vtx[2];
			}
		}
		Plane Polygon::getEdgePlane(int n) const {
			Vec3 tv = _vtx[(n+1)%3] - _vtx[n];
			tv = getNormal() % tv;
			tv.normalize();

			return Plane::FromPtDir(_vtx[n], tv);
		}
		bool Polygon::isOnTriangleSpace(const Vec3& p) const {
			// エッジを含みポリゴンと垂直な平面を3つ定義
			const float epsilon = -1e-5f;
			for(int i=0 ; i<3 ; i++) {
				if(getEdgePlane(0).dot(p) <= epsilon)
					return false;
			}
			return true;
		}
		std::pair<int,int> Polygon::split(Polygon (&dst)[3], const Plane& plane) const {
			ColCv cv = ColCv::FromVtx(_vtx[0], _vtx[1], _vtx[2]), cvd[2];
			cv.split(plane, cvd[0], cvd[1]);

			int cur = 0,
				pcur[2] = {};

			for(int i=0 ; i<2 ; i++) {
				int nP = pcur[i] = cvd[i].getNPoly();
				for(int j=0 ; j<nP ; j++)
					dst[cur++] = cvd[i].getPolygon(j);
			}
			return std::make_pair(pcur[0], pcur[1]);
		}
	}
}
