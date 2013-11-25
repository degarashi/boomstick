#include "convex.hpp"

#define DEF_TEMP	template <class VType, class UD>
#define DEF(ret)	DEF_TEMP ret Convex<VType,UD>
namespace boom {
	namespace geo3d {
		DEF(void)::init() {
			_vtx.clear();
			_rflg = RFLG_ALL;
		}
		DEF_TEMP Convex<VType,UD>::Convex() { init(); }
		DEF_TEMP Convex<VType,UD>::Convex(int n) {
			init();
			setNVtx(n);
		}
		DEF_TEMP Convex<VType,UD>::Convex(Convex&& c) {
			init();
			swap(c);
		}
		DEF_TEMP Convex<VType,UD>::Convex(const Convex& c) {
			init();
			(*this) = c;
		}

		DEF(void)::setNVtx(int nv) {
			_vtx.resize(nv);
		}
		DEF(void)::degeneration() {
			int nVtx = _vtx.size();
			if(spn::Bit::ChClear(_rflg, RFLG_DEGENERATE)) {
				// 頂点が2つ以下の場合は何もしない
				int nv = getNVtx();
				if(nv < 2)
					return;

				int wcur = 1;
				const Vec3* tv = &getPos(0);
				for(int i=1 ; i<nv ; i++) {
					if(tv->dist_sq(getPos(i)) > spn::Square(Point::NEAR_THRESHOLD)) {
						_vtx[wcur++] = _vtx[i];
						tv = &getPos(i);
					}
				}
				if(nv != wcur)
					spn::Bit::Set(_rflg, RFLG_CENTER|RFLG_CPOINT);

				nVtx = wcur;

				for(int i=0 ; i<wcur ; i++)
					AssertP(Trap, !getPos(i).isNaN(), "Poly::removeDegenerate()\n重複頂点を削除した後に重複頂点が発見されました")

				// 始点と終端チェック
				if(wcur > 1 && getPos(0).dist_sq(getPos(wcur-1)) <= spn::Square(Point::NEAR_THRESHOLD)) {
					nVtx--;
					spn::Bit::Set(_rflg, RFLG_CENTER|RFLG_CPOINT);
				}
				setNVtx(nVtx);
				if(checkLinear())
					init();
			}
		}
		DEF_TEMP
		typename Convex<VType,UD>::FList Convex<VType,UD>::checkPlaneFlags(const Plane& plane) const {
			int nv = getNVtx();
			FList fl(nv);
			for(int i=0 ; i<nv ; i++)
				fl[i] = plane.dot(_vtx[i]);
			return fl;
		}

		DEF(void)::adjustToPlane(const Plane& plane) {
			// 平面に頂点を寄せる
			Vec3 nml(plane.a, plane.b, plane.c);
			int nv = getNVtx();
			for(int i=0 ; i<nv ; i++) {
				Vec3& v = (Vec3&)_vtx[i];
				float d = plane.dot(v);
				v += (nml * -d);
			}
			spn::Bit::Set(_rflg, RFLG_ALL);
		}

		DEF(UD&)::refUserData() {
			return _ud;
		}
		DEF(const UD&)::getUserData() const {
			return _ud;
		}
		DEF(void)::setVertexArray(const VType* src, int nv) {
			setNVtx(nv);
			std::memcpy(&_vtx[0], src, sizeof(VType)*nv);
			spn::Bit::Set(_rflg, RFLG_ALL);
		}

		DEF(void)::splitThis(const Plane& plane, Convex& bDst) {
			Convex fPol;
			split(plane, fPol, bDst);
			*this = fPol;
			_rflg = RFLG_CENTER|RFLG_DEGENERATE|RFLG_CPOINT;
		}
		DEF(void)::splitThis(const Plane& plane) {
			Convex bPol;
			splitThis(plane, bPol);
		}
		DEF(void)::split(const Plane& plane, Convex& fDst, Convex& bDst) const {
			fDst.init();
			bDst.init();

			// 面を許容値を超えて跨いでいる頂点があるか
			// シビアに分割
			// 分割しない
			// すべて分割平面の表か裏として扱う

			// 簡易実装: 面を許容値を超えて跨いでいる頂点があるか判定
			int nv = getNVtx();
			bool ffSide = false,
				ffResult = false;
			for(int i=0 ; i<nv ; i++) {
				float dot = plane.dot(getPos(i));
				if(std::fabs(dot) > DOT_TOLERANCE) {
					ffResult = true;
					break;
				}
				if(i == nv-1) {
					ffResult = false;
					Vec3 tv = *spn::NormalFromPoints(_vtx[0], _vtx[1], _vtx[2]);
					if(Vec3(plane.a, plane.b, plane.c).dot(tv) >= -0.5f) {
						// 全部表
						ffSide = true;
					}
					else {
						// 全部裏
						ffSide = false;
					}
				}
			}

			// ffResult: ポリゴンが面を跨いでいればtrue
			if(!ffResult) {
				if(ffSide)
					fDst.setVertexArray(&_vtx[0], getNVtx());
				else
					bDst.setVertexArray(&_vtx[0], getNVtx());
				return;
			}

			std::vector<uint32_t> flags(getNVtx()+1);
			for(int i=0 ; i<nv ; i++) {
				float dot = plane.dot(getPos(i));
				if(dot > 0.0f)
					flags[i] = 0;
				else
					flags[i] = 1;
			}
			flags[nv] = flags[0];

			// 末尾に最初の頂点をコピー
			Convex tc(*this);
			auto tmp = tc.getVtx(0);
			tc.addVtx(tmp);
			for(int i=0 ; i<tc.getNVtx() ; i++)
				AssertP(Trap, !tc.getVtx(i).isNaN() && tc.getVtx(i).length() <= 1e8f)

			for(int i=0 ; i<nv ; i++) {
				uint32_t flag = flags[i] | (flags[i+1] << 1);
				switch(flag) {
					// 両方表の場合は、2つ目の頂点を表リストに加える
					case 0x00:
						fDst.addVtx(tc.getVtx(i+1));
						break;

					// 裏→表の場合はクリップ頂点と2つ目の頂点をフロントへ、
					// クリップ頂点をバックへ加える
					case 0x01: {
							float pd0 = std::fabs(plane.dot(tc.getPos(i))),
								pd1 = std::fabs(plane.dot(tc.getPos(i+1)));
							auto clipVtx = Segment(tc.getVtx(i), tc.getVtx(i+1)).crossPoint(plane);
							if(!clipVtx) {
								// クリップ頂点座標は平面に近いほうを採用
								if(pd0 < pd1) {
									auto& cv = tc.getVtx(i);

									fDst.addVtx(cv);
									fDst.addVtx(tc.getVtx(i+1));
									bDst.addVtx(cv);
								}
								else {
									auto& cv = tc.getVtx(i+1);

									fDst.addVtx(cv);
									bDst.addVtx(cv);
								}
							}
							else {
								fDst.addVtx(*clipVtx);
								fDst.addVtx(tc.getVtx(i+1));
								bDst.addVtx(*clipVtx);

								AssertP(Trap, std::fabs(plane.dot(*clipVtx)) < DOT_TOLERANCE)
							}
						}
						break;

					// 表→裏の場合はクリップ頂点をフロントへ、
					// クリップ頂点と2つ目の頂点をバックへ加える
					case 0x02: {
							float pd0 = std::fabs(plane.dot(tc.getPos(i))),
								pd1 = std::fabs(plane.dot(tc.getPos(i+1)));
							auto clipVtx = Segment(tc.getVtx(i), tc.getVtx(i+1)).crossPoint(plane);
							if(!clipVtx) {
								// クリップ頂点座標は平面に近いほうを採用
								if(pd0 < pd1) {
									auto& cv = tc.getVtx(i);

									fDst.addVtx(cv);
									bDst.addVtx(cv);
									bDst.addVtx(tc.getVtx(i+1));
								}
								else {
									auto& cv = tc.getVtx(i+1);

									fDst.addVtx(cv);
									bDst.addVtx(cv);
								}
							}
							else {
								fDst.addVtx(*clipVtx);
								bDst.addVtx(*clipVtx);
								bDst.addVtx(tc.getVtx(i+1));

								AssertP(Trap, std::fabs(plane.dot(*clipVtx)) < DOT_TOLERANCE)
							}
						}
						break;

					// 裏→裏の場合は2つ目の頂点をバックへ加える
					case 0x03:
						bDst.addVtx(tc.getVtx(i+1));
						break;

					default:
						AssertP(Trap, false, "予期せぬパス");
				}
			}

			// ユーザーデータはそのままコピー
			fDst.refUserData() = bDst.refUserData() = getUserData();

			// (厳密なポリゴン分割チェック)
			if(!fDst.checkValid())
				fDst.init();
			else {
				fDst.adjustToPlane(getPlane());
				if(!fDst.checkValid())
					fDst.init();
			}

			if(!bDst.checkValid())
				bDst.init();
			else {
				bDst.adjustToPlane(getPlane());
				if(!bDst.checkValid())
					bDst.init();
			}
		}
		DEF(bool)::checkValid() {
			degeneration();
			if(getNVtx() <= 2)
				return false;
			return !checkLinear();
		}
		DEF(bool)::checkLinear() const {
			// ポリゴンが一直線になっていないかチェック
			int nv = getNVtx();
			Vec3 iniDir, dir;
			iniDir = _vtx[1] - _vtx[0];
			iniDir.normalize();
			for(int i=2 ; i<nv ; i++) {
				dir = _vtx[i] - _vtx[0];
				dir.normalize();
				if(std::fabs(dir.dot(iniDir)) < LINEAR_TOLERANCE)
					return false;
			}
			return true;
		}
		DEF(bool)::checkConvex() const {
			int nV = getNVtx();
			if(nV < 3)
				return false;		// 2頂点以下は無効
			else if(nV == 3)
				return true;		// 3頂点なら無条件で凸

			// ポリゴン中心方向ベクトル
			Vec3 tv = getPos(1) - getPos(0),
				tv2 = getPos(2) - getPos(1), cv;
			cv = tv % tv2;
			cv = cv % tv2;
			cv.normalize();
			tv = tv2;

			// 右回りに辺が繋がってるかチェック
			for(int i=3 ; i<nV ; i++) {
				tv2 = getVtx(i) - getVtx(i-1);
				if(tv2.dot(cv) < -1e-2f)
					return false;

				// 中心方向ベクトルの更新
				cv = tv % tv2;
				cv = cv % tv2;
				cv.normalize();
				tv = tv2;
			}

			// 最終頂点が辺0->1を含みポリゴン中心方向を法線とする平面の表側なら凸ポリゴンである
			tv = getNormal() % (getPos(1) - getPos(0));
			tv.normalize();
			Plane plane = Plane::FromPtDir(getPos(0), tv);

			return plane.dot(getPos(nV-1)) >= 0;
		}
		DEF(PsType)::detectPlaneSide(const Plane& plane) const {
			int nv = getNVtx();
			uint32_t tst = 0;
			for(int i=0 ; i<nv ; i++) {
				float dot = plane.dot(getPos(i));
				if(dot > DOT_TOLERANCE)
					tst |= 0x01;
				else if(dot < -DOT_TOLERANCE)
					tst |= 0x02;
			}

			switch(tst) {
				case 0x00:	// すべて頂点が平面上に乗っている
					return PsType::OnPlane;
				case 0x01:	// 面の表側にある
					return PsType::Front;
				case 0x02:	// 面の裏側にある
					return PsType::Back;
				case 0x03:	// 面をまたいでいる
					return PsType::Bridge;
				default:
					AssertP(Trap, false, "Poly::detectPlaneSize()\n想定外のパターン");
					return PsType::Invalid;
			}
		}
		DEF_TEMP Convex<VType,UD> Convex<VType,UD>::FromPlane(const Plane& plane, float dist) {
			Convex ret(4);
			Vec3 tv(plane.a, plane.b, plane.c), tv2, tv3;
			tv.normalize();
			ret._vNormal = tv;
			tv2 = Vec3(1,0,0) % tv;
			tv2.normalize();
			if(tv2.length() < 0.1f)
				tv2 = tv % Vec3(0,0,1);

			tv2.normalize();
			tv3 = ret._vNormal % tv2;
			tv3.normalize();
			tv2 = tv2 * dist;
			tv3 = tv3 * dist;

			Vec3 offset = tv * -plane.d;

			ret.refVtx(0) = offset - tv3 + tv2;
			ret.refVtx(1) = offset + tv3 + tv2;
			ret.refVtx(2) = offset + tv3 - tv2;
			ret.refVtx(3) = offset - tv3 - tv2;

			Vec3 v0 = ret.refVtx(1) - ret.refVtx(0),
				v1 = ret.refVtx(2) - ret.refVtx(0), v2;
			v2 = v0 % v1;
			v2.normalize();
			if(v2.dot(tv) <= 0.0f) {
				std::swap(ret._vtx[0], ret._vtx[3]);
				std::swap(ret._vtx[1], ret._vtx[2]);
			}
			return ret;
		}
		DEF_TEMP Convex<VType,UD> Convex<VType,UD>::FromVtx(const VType& v0, const VType& v1, const VType& v2, const VType& v3) {
			Convex ret(4);
			ret.refVtx(0) = v0;
			ret.refVtx(1) = v1;
			ret.refVtx(2) = v2;
			ret.refVtx(3) = v3;
			return ret;
		}
		DEF_TEMP Convex<VType,UD> Convex<VType,UD>::FromVtx(const VType& v0, const VType& v1, const VType& v2) {
			Convex ret(3);
			ret.refVtx(0) = v0;
			ret.refVtx(1) = v1;
			ret.refVtx(2) = v2;
			return ret;
		}
		DEF(void)::addVtx(const VType& v) {
			AssertP(Trap, !v.isNaN() && v.length() <= 1e8f)
			int nv = getNVtx();
			setNVtx(nv+1);
			refVtx(nv) = v;
			spn::Bit::Set(_rflg, RFLG_ALL);
		}
		DEF(void)::popVtx() {
			int nv = getNVtx();
			if(nv > 0)
				setNVtx(nv-1);
			spn::Bit::Set(_rflg, RFLG_CENTER|RFLG_CPOINT);
		}
		DEF(int)::getNVtx() const {
			return _vtx.size();
		}
		DEF(const Vec3&)::getCPoint() const {
			if(spn::Bit::ChClear(_rflg, RFLG_CPOINT)) {
				// 各軸の最大と最小値を計算
				Vec3 tmin(1e7f, 1e7f, 1e7f), tmax(-1e7f, -1e7f, -1e7f);
				int nv = getNVtx();
				for(int i=0 ; i<nv ; i++) {
					const Vec3& v = _vtx[i];

					if(tmin.x > v.x)
						tmin.x = v.x;
					if(tmin.y > v.y)
						tmin.y = v.y;
					if(tmin.z > v.z)
						tmin.z = v.z;

					if(tmax.x < v.x)
						tmax.x = v.x;
					if(tmax.y < v.y)
						tmax.y = v.y;
					if(tmax.z < v.z)
						tmax.z = v.z;
				}
				_vCPoint = (tmin + tmax) / 2;
			}
			return _vCPoint;
		}
		DEF(const Vec3&)::getCenter() const {
			if(spn::Bit::ChClear(_rflg, RFLG_CENTER)) {
				// 中心座標の計算
				Vec3 c(0,0,0);
				int nv = getNVtx();
				for(int i=0 ; i<nv ; i++)
					c += getPos(i);
				c /= nv;
				_vCenter = c;
			}
			return _vCenter;
		}
		DEF(const Vec3&)::getNormal() const {
			if(spn::Bit::ChClear(_rflg, RFLG_NORMAL)) {
				if(getNVtx() <= 2) {
					AssertP(Trap, false, "Poly::calcNormal()\n頂点数が2個以下");
					_vNormal *= 0;
					return _vNormal;
				}

				// 法線の計算
				_vNormal = Vec3(getVtx(1)-getVtx(0)) % (getVtx(2)-getVtx(0));
				_vNormal.normalize();
			}
			return _vNormal;
		}
		DEF(Vec3)::getSumNormal() const {
			Vec3 sub0 = getVtx(1) - getVtx(0),
				sub1 = getVtx(2) - getVtx(0);
			Vec3 sum(0,0,0);
			int nv = getNVtx();
			for(int i=2 ; i<nv ; i++)
				sum += Vec3(getVtx(i-1) - getVtx(0)) % (getVtx(i) - getVtx(0));

			sum.normalize();
			AssertP(Trap, !(std::fabs(sum.x) <= 0.0f &&
				std::fabs(sum.y) <= 0.0f &&
				std::fabs(sum.z) <= 0.0f), "Poly::calcNormal()\n法線の長さが0");
			return sum;
		}
		DEF(const VType&)::getVtx(int n) const {
			AssertP(Trap, spn::IsInRange(n, 0, getNVtx()-1), "不正なインデックス(%d[%d])", getNVtx(), n);
			return _vtx[n];
		}
		DEF(const Vec3&)::getPos(int n) const {
			return (const Vec3&)getVtx(n);
		}
		DEF(VType&)::refVtx(int n) {
			return const_cast<VType&>(getVtx(n));
		}
		DEF(Vec3&)::refPos(int n) {
			spn::Bit::Set(_rflg, RFLG_ALL);
			return (Vec3&)refVtx(n);
		}
		DEF(int)::getNPoly() const {
			return std::max(0, getNVtx()-2);
		}
		DEF(void)::extractPolygons(VList& vl, IndexList& il) const {
			// 頂点をコピー
			int nV = getNVtx();
			int baseID = vl.size();
			for(int i=0 ; i<nV ; i++) {
				VType vt;
				vt = (const VType&)_vtx[i];
				vl.push_back(vt);

				if(i>=2) {
					il.push_back(baseID);
					il.push_back(baseID+i-1);
					il.push_back(baseID+i);
				}
			}
		}
		DEF_TEMP
		Convex<VType,UD>& Convex<VType,UD>::operator = (const Convex& c) {
			// copy vertex
			_vtx = c._vtx;
			// copy flag
			_rflg = c._rflg;
			// copy cache
			_vCenter = c._vCenter;
			_vNormal = c._vNormal;
			_vCPoint = c._vCPoint;
			_plane = c._plane;

			// copy userdata
			refUserData() = c.getUserData();
			return *this;
		}
		DEF_TEMP
		Vec3x2 Convex<VType,UD>::support(const Vec3& vc) const {
			const Vec3& center = getCenter();
			const Vec3& normal = getNormal();
			const Plane& plane = getPlane();
			Vec3 opvc = plane.placeOnPlane(vc, 0);
			Vec3 dir = vc-center;
			dir.normalize();
			Vec3 vanother = center+dir*100.0f, cp;
			Vec3 dstcpos, dstnml;
			float dist = 1e8f;
			int nv = getNVtx();
			for(int i=0 ; i<nv ; i++) {
				Plane tp = getEdgePlane(i);
				if(auto cp = Segment(center, vanother).crossPoint(tp)) {
					float td = center.dist_sq(*cp);
					if(dist > td) {
						dist = td;
						dstcpos = *cp;
						dstnml = getEdgeNormal(i);
					}
				}
			}
			return Vec3x2(dstcpos, dstnml);
		}
		DEF(Plane)::getEdgePlane(int eid) const {
			Vec3 tv;
			tv = getPos((eid+1)%getNVtx()) - getPos(eid);
			tv = getNormal() % tv;
			tv.normalize();
			return Plane::FromPtDir(getPos(eid), tv);
		}
		DEF(Vec3)::getEdgeNormal(int eid) const {
			Vec3 normal = getVtx((eid+1)%getNVtx()) - getVtx(eid);
			normal = normal % getNormal();
			normal.normalize();
			return normal;
		}
		DEF(bool)::selectSuitableOrigin() {
			int nV = getNVtx();
			AssertP(Trap, nV > 2, "Poly::selectSuitableOrigin()")
			if(nV == 3)
				return false;

			int idx = -1;
			float minDif = 0;
			for(int i=0 ; i<nV ; i++) {
				const Vec3 &v0 = getVtx(i),
							&v1 = getVtx((i+1)%nV);
				Vec3 tv0 = v1-v0;
				tv0.normalize();

				float curAC = 0.0f,
					tMinDif = 1.0f;
				for(int j=(i+2)%nV, count=0 ; count<nV-2 ; j=(j+1)%nV, count++) {
					const Vec3& v = getVtx(j);
					Vec3 tv1 = v-v0;
					tv1.normalize();

					float tac = std::acos(tv0.dot(tv1)) / spn::PI;
					float dif = tac - curAC;
					if(dif < tMinDif)
						tMinDif = dif;
					curAC = tac;
				}

				if(tMinDif > minDif) {
					minDif = tMinDif;
					idx = i;
				}
			}
			if(idx >= 0 && idx < nV) {
				selectOrigin(idx);
				return true;
			}
			return false;
		}
		DEF(void)::selectOrigin(int idx) {
			Convex tp;
			int nV = getNVtx();
			int rcur = idx,
				wcur = 0;
			do {
				tp._vtx[wcur++] = _vtx[rcur];
				rcur = (rcur+1) % nV;
			} while(rcur != idx);
			(*this) = tp;
		}
		DEF(void)::invert() {
			spn::Bit::Set(_rflg, RFLG_NORMAL|RFLG_PLANE);
			int cur0=0, cur1=getNVtx()-1;
			while(cur0 < cur1) {
				std::swap(_vtx[cur0], _vtx[cur1]);
				cur0 ++;
				cur1 --;
			}
		}
		DEF_TEMP
		Convex<VType,UD> Convex<VType,UD>::loopExtract(int vtxI) const {
			int nV = getNVtx();
			Convex ret(nV);
			AssertP(Trap, vtxI >= 0 && vtxI < nV, "Poly::loopExtract");
			for(int i=vtxI, count=0 ; count<nV-1 ; i=(i+1)%nV, count++)
				ret.refVtx(i) = getVtx(i);
			return ret;
		}
		DEF_TEMP
		std::pair<typename Convex<VType,UD>::VList, IndexList> Convex<VType,UD>::subdivide(float dot_area) const {
			VList dstV;
			IndexList dstI;

			// 軸探し
			Vec3 axis[3];
			axis[2] = getNormal();
			axis[0] = Vec3(0,1,0) % axis[2];
			if(axis[0].length() < 5e-5f)
				axis[0] = Vec3(1,0,0) % axis[2];
			axis[0].normalize();
			axis[1] = axis[2] % axis[0];
			axis[1].normalize();

			// ワールドtoローカル行列計算
			Vec3 tmin, tmax;
			AMat43 mLoc;
			mLoc.setColumn(0, axis[0].asVec4(0));
			mLoc.setColumn(1, axis[1].asVec4(0));
			mLoc.setColumn(2, axis[2].asVec4(0));
			std::tie(tmin, tmax) = calcMinMax(axis[0], axis[1], axis[2]);
			reinterpret_cast<AMat44&>(mLoc).transpose();

			float width = tmax.x - tmin.x,
				height = tmax.y - tmin.y;
			int tw = std::max(1.0f, width / dot_area),
				th = std::max(1.0f, height / dot_area);
			Plane divpX, divpY;
			tmin = tmin.asVec4(1) * mLoc;
			divpX = Plane::FromPtDir(tmin, axis[0]);
			Convex tpX = *this, tpbX,
					tpY, tpbY;
			divpX.d += 1e-3f;
			// X軸方向千切り
			for(int i=0 ; i<tw ; i++) {
				// 平面をX軸方向にdot_area分ずらす
				divpX.d -= dot_area;
				if(i < tw-1)
					tpX.splitThis(divpX, tpbX);
				else
					tpbX = tpX;

				tpY = tpbX;
				divpY = Plane::FromPtDir(tmin, axis[1]);
				// Y軸方向千切り
				divpY.d += 1e-3f;
				for(int j=0 ; j<th ; j++) {
					// 平面をY軸方向にdot_area分ずらす
					divpY.d -= dot_area;
					if(j==th-1) {
						// tpbをそのまま出力
						tpbY = tpY;
					}
					else {
						// tpbを平面で分割して出力
						tpY.splitThis(divpY, tpbY);
					}
					if(tpbY.getNPoly() == 0)
						continue;

					tpbY.extractPolygons(dstV, dstI);
					if(tpY.getNPoly() == 0)
						break;
				}
			}
			return std::make_pair(std::move(dstV), std::move(dstI));
		}
		DEF(float)::calcRadius() const {
			// 中心座標計算
			const Vec3& center = getCenter();
			float mDist = 0.0f;
			int nv = getNVtx();
			for(int i=0 ; i<nv ; i++) {
				Vec3 tv = getPos(i) - center;
				float tDist = tv.dot(tv);
				mDist = std::max(mDist, tDist);
			}
			return sqrt(mDist);
		}
		DEF(Sphere)::calcSphere() const {
			Sphere s;
			s.radius = calcRadius();
			s.center = getCenter();
			return s;
		}
		DEF_TEMP
		Vec3x2 Convex<VType,UD>::calcMinMax() const {
			Vec3 vmin(1e32f, 1e32f, 1e32f),
				vmax(-1e32f, -1e32f, -1e32f);
			int nv = getNVtx();
			for(int i=0 ; i<nv ; i++) {
				const Vec3& vtx = getPos(i);
				if(vmin.x > vtx.x)
					vmin.x = vtx.x;
				if(vmax.x < vtx.x)
					vmax.x = vtx.x;

				if(vmin.y > vtx.y)
					vmin.y = vtx.y;
				if(vmax.y < vtx.y)
					vmax.y = vtx.y;

				if(vmin.z > vtx.z)
					vmin.z = vtx.z;
				if(vmax.z < vtx.z)
					vmax.z = vtx.z;
			}

			return Vec3x2(vmin,vmax);
		}
		DEF(void)::addOfsVec(const Vec3& ofs) {
			int nv = getNVtx();
			for(int i=0 ; i<nv ; i++)
				refPos(i) += ofs;
			spn::Bit::Set(_rflg, RFLG_CENTER|RFLG_CPOINT);
		}
		DEF_TEMP
		Convex<VType,UD> Convex<VType,UD>::cloneReverse() const {
			int nv = getNVtx();
			Convex c(nv);
			for(int i=0 ; i<nv ; i++)
				c.refVtx(i) = getVtx(nv-1-i);
			c.refUserData() = getUserData();
			return c;
		}
		DEF_TEMP
		spn::Optional<Vec3> Convex<VType,UD>::lineCollision(const Vec3& begin, const Vec3& end, float offset, bool planeThres) const {
			Vec3 cp;
			// ポリゴンを含む平面と交わるか判定
			const Vec3& normal = getNormal();
			Plane plane = getPlane();
			plane.d -= offset;
			float bd = plane.dot(begin);
			float ed = plane.dot(end);
			if(planeThres &&
				(std::fabs(bd) < 1e-4f ||
				std::fabs(ed) < 1e-4f))
				return spn::none;
			if(bd * ed >= 0.0f)
				return spn::none;

			// 平面との交点を求める
			bd = fabs(bd);
			ed = fabs(ed);

			cp = end - begin;
			cp = cp * (bd / (bd + ed));
			cp = cp + begin;

			AssertP(Trap, std::fabs(plane.dot(cp)) <= 0.001f)

			// 三角ポリゴンごとに1つ1つ判定
			Vec3 tv0, tv1;
			int nV = getNVtx();
			for(int i=0 ; i<nV-2 ; i++) {
				Vec3 tmp0 = getVtx(0) + normal*offset;
				Vec3 tmp1 = getVtx(i+1) + normal*offset;
				Vec3 tmp2 = getVtx(i+2) + normal*offset;

				tv0 = tmp1 - tmp0;
				tv0.normalize();
				tv1 = tmp2 - tmp0;
				tv1.normalize();
				if(spn::IsInTriangle(tmp0, tmp1, tmp2, cp))
					return cp;
			}
			return spn::none;
		}
		DEF(bool)::hit(const Sphere& s) const {
			// Convex平面と球の距離が半径より大きかったらfalse
			const Plane& p = getPlane();
			if(p.dot(s.center) > s.radius)
				return false;

			// 各ポリゴンと判定
			int np = getNPoly();
			for(int i=0 ; i<np ; i++) {
				getPos(0);
				getPos(i+1);
				getPos(i+2);
			}
			return false;
		}
		DEF(Polygon)::getPolygon(int n) const {
			return Polygon(getPos(0), getPos(n+1), getPos(n+2));
		}
		DEF(const Plane&)::getPlane() const {
			if(spn::Bit::ChClear(_rflg, RFLG_PLANE)) {
				const Vec3& nml = getNormal();

				// 平面の計算
				int nv = getNVtx();
				Vec3 sum(0,0,0);
				for(int i=0 ; i<nv ; i++)
					sum += getPos(i);
				sum /= nv;

				_plane = Plane::FromPtDir(sum, nml);
				AssertP(Trap, !(std::fabs(_plane.a) <= 0.0f &&
					std::fabs(_plane.b) <= 0.0f &&
					std::fabs(_plane.c) <= 0.0f) &&
					isOnPlane(_plane, 1e-2f), "Poly::calcPlane()\n不正な平面が生成されました")
			}
			return _plane;
		}
		DEF(bool)::isOnPlane(const Plane& p, float th) const {
			FList fl = checkPlaneFlags(_plane);
			int nv = getNVtx();
			for(int i=0 ; i<nv ; i++) {
				if(fl[i] > th)
					return false;
			}
			return true;
		}
		DEF(void)::swap(Convex& c) throw() {
			std::swap(_vtx, c._vtx);
			std::swap(_ud, c._ud);
			std::swap(_vCenter, c._vCenter);
			std::swap(_vNormal, c._vNormal);
			std::swap(_vCPoint, c._vCPoint);
			std::swap(_plane, c._plane);
			std::swap(_rflg, c._rflg);
		}
		// 適時実体化
		template class Convex<spn::Vec3, ConvexUD_Col>;
	}
}
