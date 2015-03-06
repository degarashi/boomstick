#include "geom2D.hpp"

namespace boom {
	namespace geo2d {
		Vec2 MinkowskiSub(const IModel& m0, const IModel& m1, const Vec2& dir) {
			return m0.im_support(dir) - m1.im_support(-dir);
		}
		// ------------------ GSimplex ------------------
		namespace {
			constexpr int ToSearch[3][3] = {
				{2,0,1},
				{0,1,2},
				{1,2,0}
			};
		}
		void GSimplex::_minkowskiSub(const Vec2& dir, int n) {
			_posB[n] = _m1.im_support(-dir);
			_poly.point[n] = (_m0.im_support(dir) - _posB[n]);
		}
		bool GSimplex::_gjkMethod() {
			// とりあえず物体の中心位置でサポートベクトルを決める
			Vec2 dir(_m1.im_getCenter() - _m0.im_getCenter());
			float lens = dir.len_sq();
			if(lens < NEAR_THRESHOLD_SQ)
				dir = Vec2(1,0);
			else
				dir *= spn::RSqrt(lens);

			_minkowskiSub(dir, 0);
			if(dir.dot(_poly.point[0]) < -DOT_THRESHOLD) {
				_nVtx = 1;
				return false;
			}
			// 原点と重なっていたら終了 = 内部点
			lens = _poly.point[0].len_sq();
			if(lens < NEAR_THRESHOLD) {
				_nVtx = 1;
				_inner = _posB[0];
				return true;
			}
			dir = _poly.point[0] * spn::RSqrt(lens) * -1.f;
			_minkowskiSub(dir, 1);

			const Vec2 ori(0);
			Vec2 tmp = _poly.point[1] - _poly.point[0];
			float tmpLen = tmp.length();
			tmp *= spn::Rcp22Bit(tmpLen);
			float r = tmp.dot(-_poly.point[0]);
			if(r > tmpLen)
				return false;
			tmp = _poly.point[0] + tmp*r;

			lens = tmp.len_sq();
			if(lens < NEAR_THRESHOLD) {
				// ライン上に原点がある
				_inner =_posB[0].l_intp(_posB[1], r * spn::Rcp22Bit(tmpLen));
				_nVtx = 2;
				return true;
			}

			_nVtx = 3;
			LNear res(tmp, LinePos::OnLine);
			int idx = 2;
			for(;;) {
				auto& rV = _poly.point[idx];
				// 新たな頂点を追加
				dir = res.first.normalization() * -1.f;
				_minkowskiSub(dir, idx);
				if(dir.dot(rV) < -DOT_THRESHOLD)
					return false;

				auto& ts = ToSearch[idx];
				if(rV.distance(_poly.point[ts[0]]) < 1e-7f ||
					rV.distance(_poly.point[ts[2]]) < 1e-7f)
				{
					std::cout << _m0 << std::endl << _m1 << std::endl;
					return false;
				}

				// 現時点で三角形が原点を含んでいるか
				bool bIn;
				if((_poly.point[1]-_poly.point[0]).cw(_poly.point[2]-_poly.point[0]) > 0)
					bIn = _poly.hit({0,0});
				else {
					_poly.invert();
					bIn = _poly.hit({0,0});
					_poly.invert();
				}
// 				float cm = (_poly.point[ts[1]] - _poly.point[ts[0]]).cw(-_poly.point[ts[0]]) *
// 							(_poly.point[ts[2]] - _poly.point[ts[1]]).cw(-_poly.point[ts[1]]);
// 				if(cm >= -1e-5f) {
				if(bIn) {
					float cf[3];
					spn::BarycentricCoord(cf, _poly.point[0], _poly.point[1], _poly.point[2], ori);
					_inner = spn::MixVector(cf, _posB[0], _posB[1], _posB[2]);
					_inner = TriangleLerp(_poly.point[0], _poly.point[1], _poly.point[2], ori,
										_posB[0], _posB[1], _posB[2]);
					Assert(Trap, _m0.im_hitPoint(_inner) && _m1.im_hitPoint(_inner));
					return true;
				}
				float dist = std::numeric_limits<float>::max();
				idx = -1;
				res = Segment(_poly.point[ts[0]], _poly.point[ts[1]]).nearest(ori);
				if(res.second == LinePos::OnLine) {
					float td = res.first.len_sq();
					if(dist > td) {
						dist = td;
						// この辺について調べる
						idx = ts[2];
					}
				}
				res = Segment(_poly.point[ts[1]], _poly.point[ts[2]]).nearest(ori);
				if(res.second == LinePos::OnLine) {
					float td = res.first.len_sq();
					if(dist > td) {
						dist = td;
						// この辺について調べる
						idx = ts[0];
					}
				}
				if(dist < NEAR_THRESHOLD_SQ) {
					_inner = TriangleLerp(_poly.point[0], _poly.point[1], _poly.point[2], ori,
										  _posB[0], _posB[1], _posB[2]);
					return true;
				}
				if(idx < 0)
					return false;
			}
		}
		GSimplex::GSimplex(const IModel& m0, const IModel& m1): _m0(m0), _m1(m1), _bHit(_gjkMethod()) {}
		bool GSimplex::getResult() const {
			return _bHit;
		}
		const Vec2& GSimplex::getInner() const {
			return _inner;
		}
		namespace {
			bool IsNaN(float f) {
				return !(f>0) && !(f<=0);
			}
			bool IsNaN(const Vec2& v) {
				return IsNaN(v.x) || IsNaN(v.y);
			}
		}
		// ---------------------- GEpa ----------------------
		GEpa::~GEpa() {
			_clear();
		}
		GEpa::VPool thread_local GEpa::tls_vPool(GEpa::MAX_VERT);

		Vec2x2* GEpa::_allocVert(int n) {
			Vec2x2* p = tls_vPool.malloc();
			if(n >= 0) {
				std::memmove(&_vl[n+1], &_vl[n], sizeof(Vec2x2*)*(_szVl-n));
				_vl[n] = p;
			} else
				_vl[_szVl] = p;
			++_szVl;
			AssertP(Trap, _szVl <= _vl.size());
			return p;
		}
		const Vec2x2& GEpa::_minkowskiSub(const Vec2& dir, int n) {
			Vec2x2* vp = _allocVert(n);
			vp->second = _m1.im_support(-dir);
			vp->first = _m0.im_support(dir) - vp->second;
			return *vp;
		}
		// デバッグ用Print関数
		namespace {
			void PrintV(const Vec2& v) {
				std::cout << '[' << v.x << ',' << v.y << ']';
			}
			void PrintV(const Vec2x2& v) {
				PrintV(v.first);
				PrintV(v.second);
			}
		}
		void GEpa::_printASV() const {
			int nV = _szVl;
			for(int i=0 ; i<nV ; i++) {
				PrintV(*_vl[i]);
				std::cout << std::endl;
			}
			for(auto itr=_asv.cbegin() ; itr!=_asv.cend() ; itr++) {
				std::cout << itr->dist << ':';
				PrintV(itr->dir);
				std::cout << std::hex << itr->vtx[0] << ',' << itr->vtx[1] << std::endl;
			}
		}

		int GEpa::_getIndex(const Vec2x2 *vp) const {
			auto itr = std::find(_vl.begin(), _vl.end(), vp);
			if(itr == _vl.end())
				return 0xffff;
			return itr-_vl.begin();
		}
		namespace {
			const Vec2 c_origin(0);
		}
		bool GEpa::_addAsv(const Vec2x2& v0, const Vec2x2& v1) {
			auto res = Segment(v0.first, v1.first).nearest(c_origin);
			if(res.second == LinePos::OnLine) {
				float len = res.first.length();
				res.first *= spn::Rcp22Bit(len);
				_asv.insert(LmLen{len, res.first, {&v0,&v1}});
				return true;
			}
			return res.second != LinePos::NotHit;
		}
		void GEpa::_epaMethodOnHit() {
			// 衝突時: 脱出ベクトルを求める
			float minLen = std::numeric_limits<float>::max();
			while(!_asv.empty()) {
				auto fr = _asv.pop_frontR();

				const Vec2x2 &v1 = _minkowskiSub(fr.dir, _getIndex(fr.vtx[0])),
							&v0 = *fr.vtx[0],
							&v2 = *fr.vtx[1];
				float d1 = fr.dir.dot(v1.first);
				if(spn::EqAbs(d1, fr.dist, 5e-4f)) {
					if(minLen > fr.dist) {
						minLen = fr.dist;
						_pvec = fr.dir * fr.dist;
					}
					if(_asv.empty() || _asv.front().dist > fr.dist)
						break;
				}

				if(_addAsv(v0, v1) || _addAsv(v1, v2)) {
					_pvec = fr.dir * fr.dist;
					break;
				}
			}
		}
		void GEpa::_epaMethodNoHit() {
			// 無衝突時: 最近傍対を求める
			while(!_asv.empty()) {
				auto fr = _asv.pop_frontR();
				if(!fr.vtx[1]) {
					auto& p = *fr.vtx[0];
					_nvec.second = p.second;
					_nvec.first = p.second + fr.dir * fr.dist;
					return;
				}
				if(fr.dist < 1e-6f || IsNaN(fr.dir)) {
					_nvec.first = _nvec.second = fr.vtx[0]->second;
					return;
				}
				const Vec2x2& v1 = _minkowskiSub(-fr.dir, _getIndex(fr.vtx[0])),
							&v0 = *fr.vtx[0],
							&v2 = *fr.vtx[1];
				float d1 = fr.dir.dot(v1.first);
				if(spn::EqAbs(std::fabs(d1), fr.dist, 5e-4f)) {
					_nvec.second = v1.second;
					_nvec.first = v1.second + fr.dir * fr.dist;
					return;
				}
				if(_addAsv(v0, v1) || _addAsv(v1, v2)) {
					_nvec.second = v1.second;
					_nvec.first = v1.second + fr.dir * fr.dist;
					return;
				}
			}
			AssertP(Trap, false)
		}
		void GEpa::_adjustLoop3() {
			AssertP(Trap, _szVl == 3)
			Vec2 v01(_vl[1]->first - _vl[0]->first),
				v02(_vl[2]->first - _vl[0]->first);
			if(v01.ccw(v02) < 0)
				std::swap(_vl[0], _vl[1]);
		}
		void GEpa::_clear() {
			size_t nV = _szVl;
			for(size_t i=0 ; i<nV ; i++)
				tls_vPool.destroy(_vl[i]);
			_szVl = 0;
			_asv.clear();
		}
		void GEpa::_geneASV() {
			int nV = _szVl;
			assert(nV >= 3 || !getResult());
			_asv.clear();

			for(int i=0 ; i<nV ; i++) {
				auto &p0 = *_vl[i],
					&p1 = *_vl[(i+1)%nV];
				auto res = Segment(p1.first, p0.first).nearest(c_origin);
				float len = res.first.length();
				if(len < 1e-5f) {
					// 線分上に原点があるので少し小細工
					// 線分方向90度回転かつ線分に関わっていない頂点側
					Vec2 dir(p1.first - p0.first);
					dir.normalize();
					dir *= cs_mRot90[0];
					if(dir.dot(_vl[(i+2)%nV]->first - p0.first) > 0.f)
						dir *= -1.f;
					res.first = dir;
				} else
					res.first *= spn::Rcp22Bit(len);
				_asv.insert(LmLen{len, res.first, {_vl[i], _vl[(i+1)%nV]}});
			}
		}
		void GEpa::_recover2NoHit() {
			auto res = Segment(_poly.point[0], _poly.point[1]).nearest(Vec2(0));
			float len = res.first.length();
			LmLen lmlen;
			lmlen.dist = len;
			lmlen.dir *= spn::Rcp22Bit(len);
			if(res.second != LinePos::OnLine) {
				lmlen.vtx[0] = lmlen.vtx[1] = _vl[(res.second == LinePos::Begin) ? 0 : 1];
			} else {
				lmlen.vtx[0] = _vl[0];
				lmlen.vtx[1] = _vl[1];
			}
			_asv.insert(lmlen);
		}
		void GEpa::_recover2OnHit() {
			LmLen lm;
			lm.dist = 0;
			lm.dir = _poly.point[0]-_poly.point[1];
			lm.dir.normalize();
			lm.dir *= cs_mRot90[0];
			lm.vtx[0] = _vl[0];
			lm.vtx[1] = _vl[1];
			_asv.insert(lm);

			std::swap(lm.vtx[0], lm.vtx[1]);
			lm.dir *= -1.f;
			_asv.insert(lm);
		}

		GEpa::GEpa(const IModel& m0, const IModel& m1): GSimplex(m0,m1), _szVl(0) {
			int nV = _nVtx;
			_szVl = nV;
			for(int i=0 ; i<nV ; i++) {
				auto res = _vl[i] = tls_vPool.malloc();
				res->first = _poly.point[i];
				res->second = _posB[i];
			}

			do {
				if(getResult()) {
					if(nV == 1) {
						_pvec *= 0;
						break;
					} else if(nV == 2)
						_recover2OnHit();
					else {
						_adjustLoop3();
						_geneASV();
					}
					_epaMethodOnHit();
				} else {
					if(nV == 1) {
						const Vec2& v = _minkowskiSub(-_poly.point[0].normalization()).first;
						auto res = Segment(_poly.point[0], v).nearest(Vec2(0));
						if(res.second == LinePos::OnLine) {
							float len = res.first.length();
							_asv.insert(LmLen{len, res.first*spn::Rcp22Bit(len), {_vl[0],_vl[1]}});
						}
						else {
							const auto& p = *_vl[0];
							_nvec.second = p.second;
							_nvec.first = _nvec.second + p.first;
							break;
						}
					} else if(nV == 2)
						_recover2NoHit();
					else {
						_adjustLoop3();
						_geneASV();
					}
					_epaMethodNoHit();
				}
			} while(false);
		}
		Vec2x2 GEpa::getNearestPair() const {
			AssertP(Trap, !getResult())
			AssertP(Trap, !IsNaN(_nvec.first))
			AssertP(Trap, !IsNaN(_nvec.second))
			return _nvec;
		}
		const Vec2& GEpa::getPVector() const {
			AssertP(Trap, getResult())
			AssertP(Trap, !IsNaN(_pvec))
			return _pvec;
		}
		/*! 算出不能なケースは考えない */
		Float2 TriangleRatio2(const Vec2& v0, const Vec2& v1, const Vec2& v2, const Vec2& vt) {
			Vec2 toV1(v1-v0),
				toV2(v2-v0),
				toVT(vt-v0);

			const float ar[2][3] = {
				{toV1.x, toV2.x, toVT.x},
				{toV1.y, toV2.y, toVT.y}
			};
			reg128 xm[2] = {LOADPS_Z3(ar[0]), LOADPS_Z3(ar[1])};
			xm[0] = spn::_mmDivPs(xm[0], reg_load1_ps(&toV1.x));
			xm[1] = reg_add_ps(reg_mul_ps(xm[0], reg_load1_ps(&toV1.y)), xm[1]);
			xm[1] = spn::_mmDivPs(xm[1], reg_shuffle_ps(xm[1], xm[1], _REG_SHUFFLE(1,1,1,1)));
			xm[0] = reg_sub_ps(xm[0], reg_mul_ps(xm[1], reg_shuffle_ps(xm[0], xm[0], _REG_SHUFFLE(1,1,1,1))));

			Float2 ret;
			reg_store_ss(&ret.first, xm[0]);
			reg_store_ss(&ret.second, xm[1]);
			return ret;
		}
	}
}
