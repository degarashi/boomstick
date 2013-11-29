#include "geom3D.hpp"
#include "convex.hpp"

namespace boom {
	namespace geo3d {
		Vec3 MinkowskiSub(const IModel& m0, const IModel& m1, const Vec3& dir) {
			return m0.im_support(dir) - m1.im_support(-dir);
		}

		// ---------------------- GSimplex ----------------------
		GSimplex::GSimplex(const IModel& m0, const IModel& m1): _m0(m0), _m1(m1) {
			_clear();
		}
		void GSimplex::_clear() {
			_nVtx = 0;
			_retryCount = 0;
			_hr = Result::Default;
		}
		bool GSimplex::_checkDuplication(const Vec3& v) const {
			int nV = _nVtx;
			for(int i=0 ; i<nV ; i++) {
				if(_vtx[i].dist_sq(v) < 1e-4f)
					return true;
			}
			return false;
		}
		bool GSimplex::_addVtx(const Vec3& vA, const Vec3& vB, int id) {
			AssertP(Trap, _nVtx <= 4)
			Vec3 v = vA - vB;
			if(_checkDuplication(v))
				return true;

			if(_nVtx == 4) {
				// 不要な頂点を省く
				AssertP(Trap, id >= 0)
				_vtx[id] = v;
				_posB[id] = vB;
			} else {
				_vtx[_nVtx] = v;
				_posB[_nVtx++] = vB;
			}
			return false;
		}

		const Vec3 GSimplex::cs_randVec[6] = {Vec3(1,0,0), Vec3(0,1,0), Vec3(0,0,1),
												Vec3(-1,0,0), Vec3(0,-1,0), Vec3(0,0,-1)};
		std::tuple<Vec3,GSimplex::Result,int> GSimplex::getNearestDir() {
			AssertP(Trap, _nVtx > 1)
			const Vec3 ori(0,0,0);
			if(_nVtx == 2) {
				if(_retryCount == 0) {
					// 線分と点の最短距離
					LNear lnear = Segment(_vtx[0], _vtx[1]).nearest(ori);
					_nearP = lnear.first;
					if(lnear.second == LinePos::OnLine) {
						float len = _nearP.length();
						if(len >= 1e-5f)
							return std::make_tuple(Vec3(_nearP/len), Result::Default, -1);
					}
				}
				while(_retryCount != countof(cs_randVec)) {
					Vec3 toV1 = _vtx[1]-_vtx[0];
					Vec3 dir = cs_randVec[_retryCount].cross(toV1);
					float len = dir.length();
					if(len > 1e-5f)
						return std::make_tuple(Vec3(dir/len), Result::Default, -1);
					++_retryCount;
				}
				return std::make_tuple(Vec3(), Result::NoHit, -1);
			}
			if(_nVtx == 3) {
				Polygon poly(_vtx[0], _vtx[1], _vtx[2]);
				if(_retryCount == 0) {
					// ポリゴンと点の最短距離
					auto res = poly.nearest(ori);
					_nearP = res.first;
					if(res.second)
						return std::make_tuple(Vec3(_nearP.normalization()), Result::Default, -1);
				}
				int rc = _retryCount;
				if(rc == 0)
					return std::make_tuple(Vec3(-poly.getNormal()), Result::Default, -1);
				else if(rc == 1)
					return std::make_tuple(poly.getNormal(), Result::Default, -1);
				return std::make_tuple(Vec3(), Result::NoHit, -1);
			}
			// 4面体と点の最短ベクトル
			float minD = std::numeric_limits<float>::max();
			bool bInPoly = true;
			if(std::fabs(Plane::FromPts(_vtx[0],_vtx[1],_vtx[2]).dot(_vtx[3])) > 1e-4f) {
				Vec3 center = (_vtx[0]+_vtx[1]+_vtx[2]+_vtx[3])/4,
					cp;
				for(int i=0 ; i<countof(cs_index) ; i++) {
					auto& idx = cs_index[i];
					Polygon poly(_vtx[idx[0]],_vtx[idx[1]],_vtx[idx[2]]);
					const auto& plane = poly.getPlane();

					if(plane.dot(center) * plane.d < -std::numeric_limits<float>::min())
						bInPoly = false;
					cp = poly.nearest(ori).first;

					float d = cp.len_sq();
					if(minD > d) {
						minD = d;
						_nearP = cp;
						_nearIDX = i;
					}
				}
			} else {
				// 4面体がつぶれている場合は別の判定方法
				bInPoly = false;
				for(int i=0 ; i<countof(cs_index) ; i++) {
					auto& idx = cs_index[i];
					Vec3 toV1 = _vtx[idx[1]]-_vtx[idx[0]],
						toV2 = _vtx[idx[2]]-_vtx[idx[0]],
						vC = (toV1%toV2).normalization();
					float det = spn::CramerDet(toV1, toV2, vC);
					Vec3 res = spn::CramersRule(toV1, toV2, vC, -_vtx[idx[0]], 1.0f/det);
					if(std::fabs(res.z) > 1e-4f ||
						!spn::IsInRange(res.x, 0.0f, 1.0f) ||
						!spn::IsInRange(res.y, 0.0f, 1.0f) ||
						!spn::IsInRange(res.x+res.y, 0.0f, 1.0f))
					{}
					else {
						bInPoly = true;
						_nearIDX = i;
						_nearP *= 0;
						break;
					}
				}
				if(!bInPoly) {
					// 線分の最短距離
					minD = std::numeric_limits<float>::max();
					for(int i=0 ; i<4 ; i++) {
						for(int j=i+1 ; j<4 ; j++) {
							Segment ls(_vtx[i], _vtx[j]);
							float r = ls.getRatio(ori);
							if(spn::IsInRange(r, 0.0f, 1.0f)) {
								float d = ls.dist_sq(ori);
								if(d < minD) {
									minD = d;
									_nearPosB = _posB[i].l_intp(_posB[j], r);
								}
							}
						}
					}
					_nearIDX = -1;
					return std::make_tuple(Vec3(_nearPosB.normalization()), Result::NoHit, -1);
				}
			}
			return std::make_tuple(Vec3(_nearP.normalization()), bInPoly ? Result::Hit : Result::Default, cs_index[_nearIDX][3]);
		}

		const int GSimplex::cs_index[4][4] = {{0,1,2, 3}, {0,3,2, 1}, {0,1,3, 2}, {2,3,1, 0}};
		//! GJKを用いて凸物体の衝突判定
		//! @return 算出結果の座標, 衝突したかのフラグ
		bool GSimplex::_gjkMethod() {
			_clear();

			// GJKにおいて support方向が前回と同じ又は既存の頂点と同じ結果を得た: false
			// 点や線分やポリゴンが原点を含んだ: true
			// 4面体が原点を内包した: true

			// 1個目の頂点
			Vec3 dir = (_m1.im_getCenter() - _m0.im_getCenter()).normalization();
			Vec3 pos0 = _m0.im_support(-dir),
				pos1 = _m1.im_support(dir);
			_vtx[0] = pos0 - pos1;
			_posB[0] = pos1;
			_nVtx = 1;
			if(_vtx[0].len_sq() < 1e-4f) {
				_hr = Result::Hit;
				return true;
			}
			dir = _vtx[0].normalization();

			int delID = -1;
			for(;;) {
				// ワールドへ変換
				pos0 = _m0.im_support(-dir);
				pos1 = _m1.im_support(dir);
				//if(_nVtx == 4 && (pos0-pos1).dot(dir) < 0) {
				//	_hr = Result::NoHit;
				//	return false;
				//}
				// 同時に不要な点を削除 (supportベクトルの算出にかかわらなかった頂点)
				if(_addVtx(pos0, pos1, delID)) {
					++_retryCount;
					if(_nVtx == 4) {
						_hr = Result::NoHit;
						return false;
					}
				} else
					_retryCount = 0;

				std::tie(dir, _hr, delID) = getNearestDir();
				if(_hr != Result::Default)
					return _hr == Result::Hit;
			}
		}
		namespace {
			bool IsValidValue(float v) { return v>=0 || v<=0; }
			bool IsValidValue(const Vec3& v) {
				return IsValidValue(v.x+v.y+v.z);
			}
		}
		Vec3 GSimplex::getInterPoint() const {
			AssertP(Trap, _hr == Result::Hit)
			if(_nVtx == 1)
				return _posB[0];
			AssertP(Trap, _nVtx == 4)

			// 直前に消された側の面からposBを補間する
			// 最寄り点ベクトルの半分を移動すれば中点が求まる
			Vec3 ori(0,0,0);
			auto& idx = cs_index[_nearIDX];
			Polygon poly(_vtx[idx[0]], _vtx[idx[1]], _vtx[idx[2]]);
			auto& nml = poly.getNormal();
			if(IsValidValue(nml)) {
				Vec3 cp = poly.nearest(ori).first;
				float d = cp.len_sq();

				Vec3 nml = poly.getNormal();
				Vec3 toV1 = poly.getVtx(1)-poly.getVtx(0),
					toV2 = poly.getVtx(2)-poly.getVtx(0);
				float detInv = 1.0f / spn::CramerDet(toV1, toV2, nml);
				Vec3 res = spn::CramersRule(toV1,toV2, nml, cp-poly.getVtx(0), detInv);
				Vec3 toP1 = _posB[idx[1]] - _posB[idx[0]],
						toP2 = _posB[idx[2]] - _posB[idx[0]];
				Vec3 ret = toP1 * res.x + toP2 * res.y + _posB[idx[0]];
				return ret;
			}
			return _posB[0];
		}
		namespace {
			Vec3 CalcPolyRatio(const Vec3& dir0, const Vec3& dir1, const Vec3& cp,
				const Vec3& od0, const Vec3& od1, const Vec3& odOri)
			{
				Vec3 nml = dir0.cross(dir1);
				float detInv = 1.0f / spn::CramerDet(dir0,dir1,nml);
				Vec3 res = spn::CramersRule(dir0, dir1, nml, cp, detInv);
				AssertP(Trap, (dir0*res.x + dir1*res.y).distance(cp) < 1e-3f)
				return od0*res.x + od1*res.y + odOri;
			}
			Vec3 CalcPolyRatio(const Vec3& src0, const Vec3& src1, const Vec3& src2, const Vec3& cp,
							const Vec3& dst0, const Vec3& dst1, const Vec3& dst2)
			{
				return CalcPolyRatio(src1 - src0,
									src2 - src0,
									cp - src0,
									dst1 - dst0,
									dst2 - dst0,
									dst0);
			}
		}
		Vec3x2 GSimplex::getNearestPair() const {
			AssertP(Trap, _hr == Result::NoHit)
			if(_nVtx == 1)
				return Vec3x2(_posB[0], Vec3(_posB[0]+_vtx[0]));
			if(_nVtx == 2) {
				// 2回分の履歴を線分の近傍点で割合をかけたもの
				// (nearPは線分上にある)
				float r = (_vtx[1] - _vtx[0]).dot(_nearP - _vtx[0]);
				AssertP(Trap, spn::IsInRange(r, 0.0f-1e-5f, 1.0f+1e-5f))
				Vec3 v0 = _posB[0].l_intp(_posB[1], r),
					v1 = v0 + _nearP;
				return Vec3x2(v0, v1);
			}
			// 4面体は最寄りの3角ポリゴンと同義
			if(_nearIDX >= 0) {
				// (nearPはポリゴン上にある)
				const int c_defaultID[3] = {0,1,2};
				const int* use = cs_index[(_nVtx==3) ? 0 : _nearIDX];
				Vec3 v0 = CalcPolyRatio(_vtx[use[0]], _vtx[use[1]], _vtx[use[2]], _nearP,
										_posB[use[0]], _posB[use[1]], _posB[use[2]]),
					v1 = v0 + _nearP;
				return Vec3x2(v0, v1);
			}
			return Vec3x2(_nearPosB, Vec3(_nearPosB+_nearP));
		}
		int GSimplex::getNVtx() const { return _nVtx; }
		const Vec3& GSimplex::getVtx(int n) const { return _vtx[n]; }
		bool GSimplex::getResult() const {
			AssertP(Trap, _hr!=Result::Default)
			return _hr == Result::Hit;
		}

		// ---------------------- EConvex ----------------------
		uint32_t EConvex::_ComposeIDX(int id0, int id1, int id2) {
			if(id2 > id1) std::swap(id2,id1);
			if(id1 > id0) std::swap(id1,id0);
			if(id2 > id1) std::swap(id2,id1);
			AssertP(Trap, id0 > id1 && id1 > id2)
			return (id0<<20) | (id1<<10) | id2;
		}
		std::tuple<int,int,int> EConvex::_DecompIDX(uint32_t id) {
			int id2 = (id & 0x3f),
				id1 = (id >> 10) & 0x3f,
				id0 = (id >> 20) & 0x3f;
			return std::make_tuple(id0, id1, id2);
		}
		EConvex::EConvex(const GSimplex& gs, const IModel& m0, const IModel& m1):
			ConvexP(&gs.getVtx(0), gs.getNVtx()), _m0(m0), _m1(m1) {}
		std::pair<Vec3,float> EConvex::_epaMethod() {
			// 点が4個に満たない場合の処理
			int nV = getNVtx();
			AssertP(Trap, spn::IsInRange(nV, 1, 4))
			if(nV == 1) {
				// 1個なら点が原点と重なっているハズ: c1からc0へのsupport写像を試す => 頂点2個のケースと同様に処理
				Vec3 p = MinkowskiSub(_m0, _m1, (_m1.im_getGCenter() - _m0.im_getGCenter()).normalization());
				if(getVtx(0).dist_sq(p) < 1e-6f)
					return std::make_pair(Vec3(0,0,0), 0.0f);
				addVtx(p);
			}
			// 既存の頂点とダブってない場合にのみMinkowski差を追加
			auto addv = [this](const Vec3& dir) -> bool{
				Vec3 p = MinkowskiSub(_m0, _m1, dir);
				int cur=0,
					nV = getNVtx();
				while(cur != nV) {
					if(getVtx(cur).dist_sq(p) < 1e-6f)
						break;
					++cur;
				}
				if(cur == nV) {
					addVtx(p);
					return true;
				}
				return false;
			};

			if(nV == 2) {
				// 2個なら辺と重なっているハズ: 辺に垂直な4方向に対してsupport写像後にquickhull
				Vec3 toV1(getVtx(1) - getVtx(0));
				Vec3 nml0 = toV1.cross(Vec3(0,1,0));
				if(nml0.len_sq() < 1e-6f)
					nml0 = toV1.cross(Vec3(1,0,0));
				nml0.normalize();
				Vec3 nml1 = toV1.cross(nml0);
				nml1.normalize();

				addv(nml0);
				addv(-nml0);
				addv(nml1);
				addv(-nml1);

				Idx3List polL = getPolyFace();
				for(auto& p : polL)
					p.flip();
				_addPoly(polL);
			} else if(nV == 3) {
				// 3個なら面と重なっているハズ: 面の法線両方に対してsupport写像後にquickhull
				Vec3 nml = Polygon(getVtx(0), getVtx(1), getVtx(2)).getNormal();

				addv(nml);
				addv(-nml);

				Idx3List polL = getPolyFace();
				for(auto& p : polL)
					p.flip();

				_addPoly(polL);
			} else {
				// DistanceMapに登録
				const int c_idx[4][3] = {{0,1,2}, {2,1,3}, {3,1,0}, {3,0,2}};
				for(auto& idx : c_idx)
					_addPoly(idx[0], idx[1], idx[2]);
			}
			AssertP(Trap, getNVtx() >= 4)

			// 最短距離の面を1つ取り出す
			auto itr = _d_id.begin();
			bool bLoop = true;
			do {
				// 面と重なっている場合は凸包の内部点と反対側の法線を使う
				const auto& poly = itr->second;
				const auto& nml = poly.getPlane().getNormal();
				bool bOnPlane = itr->first < 1e-5f;

				// ポリゴンの上空にあって距離が負数の物を採用
				if(!_dividePoly(_m0, _m1, poly, -nml)) {
					if(bOnPlane) {
						// もう片方を試してダメなら終了
						if(!_dividePoly(_m0, _m1, poly, nml))
							_d_id.erase(itr);
					} else
						bLoop = false;
				}
				AssertP(Trap, !_d_id.empty())
				itr = _d_id.begin();
			} while(bLoop);
			const auto& plane = itr->second.getPlane();
			return std::make_pair(Vec3(-plane.getNormal()), plane.d);
		}
		void EConvex::_addPoly(const Idx3List& src) {
			for(auto& p : src)
				_addPoly(p.getID(0), p.getID(1), p.getID(2));
		}
		void EConvex::_addPoly(int idx0, int idx1, int idx2) {
			uint32_t id = _ComposeIDX(idx0, idx1, idx2);
			AssertP(Trap, _id_d.count(id) == 0)

			// もし面が原点方向を向いてなかったら反転
			Idx3 poly(idx0,idx1,idx2, getVtxArray(), 0);
			auto& v = getVtxArray();
			Vec3 cp = Polygon(v[idx0],v[idx1],v[idx2]).nearest(Vec3(0,0,0)).first;
			float d = cp.length();//poly.getPlane().d;
			AssertP(Trap, d >= 0)
			d = std::fabs(d);
			_d_id.insert(std::make_pair(d, poly));
			_id_d.insert(std::make_pair(id, d));
		}
		void EConvex::_delPoly(uint32_t id) {
			auto itr = _id_d.find(id);
			// 同じ距離の面は他にあるかもしれないので線形探索
			auto itrB = _d_id.lower_bound(itr->second);
			auto itrE = _d_id.upper_bound(itr->second);
			while(itrB != itrE) {
				const auto& p = itrB->second;
				uint32_t pi_id = _ComposeIDX(p.getID(0), p.getID(1), p.getID(2));
				if(pi_id == id) {
					_d_id.erase(itrB);
					_id_d.erase(itr);
					return;
				}
				++itrB;
			}
			AssertP(Trap, false)
		}
		bool EConvex::_dividePoly(const IModel& m0, const IModel& m1, const Idx3& poly, const Vec3& dir) {
			const auto& vtx = getVtxArray();
			const auto& plane = poly.getPlane();
			const int id[3] = {poly.getID(0), poly.getID(1), poly.getID(2)};
			Polygon ppoly(vtx[id[0]], vtx[id[1]], vtx[id[2]]);
			Vec3 p = MinkowskiSub(m0, m1, dir);
			if(plane.dot(p) < 0 &&
				ppoly.isOnTriangleSpace(p))
			{
				// 既存の頂点とダブってる場合(恐らく対象となった面の構成頂点)は省く
				int nV = vtx.size();
				for(int i=0 ; i<nV ; i++) {
					if(vtx[i].dist_sq(p) < 1e-6f)
						return false;
				}

				// 古い面の削除
				_delPoly(_ComposeIDX(poly.getID(0), poly.getID(1), poly.getID(2)));
				// 新しく生成された面の追加
				// [old0,old1,new] [old1,old2,new] [old2,old0,new]
				const int new_id = nV;
				addVtx(p);
				_addPoly(id[0], id[1], new_id);
				_addPoly(id[1], id[2], new_id);
				_addPoly(id[2], id[0], new_id);
				return true;
			}
			return false;
		}
	}
}
