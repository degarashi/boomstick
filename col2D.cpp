#include "geom2D.hpp"
#include "spinner/bits.hpp"
#include "spinner/common.hpp"
#include "spinner/spn_math.hpp"
#include "spinner/misc.hpp"
#include "spinner/assoc.hpp"
#include <limits>
#include <iostream>

namespace boom {
	const spn::AMat22 cs_mRot90[2] = {
		{spn::COS90, -spn::SIN90,
		spn::SIN90, spn::COS90},
		{spn::COS90, spn::SIN90,
		-spn::SIN90, spn::COS90}
	};
	float Area_x2(const Vec2& v0, const Vec2& v1) {
		return std::fabs(v0.cw(v1));
	}
	float Area_x2(const Vec3& v0, const Vec3& v1) {
		return v0.cross(v1).length();
	}

	namespace geo2d {
		// ---------------------- GSimplex ----------------------
		Vec2 MinkowskiSub(const IModel& m0, const IModel& m1, const Vec2& dir) {
			return m0.support(dir) - m1.support(-dir);
		}

		namespace {
			constexpr float DOT_THRESHOLD = 1e-3f,
							DIST_THRESHOLD = 1e-6f;
			constexpr int ToSearch[3][3] = {
				{2,0,1},
				{0,1,2},
				{1,2,0}
			};
		}
		void GSimplex::_minkowskiSub(const Vec2& dir, int n) {
			_posB[n] = _m1.support(-dir);
			_vtx[n] = (_m0.support(dir) - _posB[n]);
		}
		bool GSimplex::_gjkMethod() {
			Vec2 dir(_m1.getCenter() - _m0.getCenter());
			float lens = dir.len_sq();
			if(lens < DIST_THRESHOLD)
				dir = Vec2(1,0);
			else
				dir *= spn::RSqrt(lens);

			_minkowskiSub(dir, 0);
			if(dir.dot(_vtx[0]) < -DOT_THRESHOLD) {
				_nVtx = 1;
				return false;
			}
			// 原点と重なっていたら終了 = 内部点
			lens = _vtx[0].len_sq();
			if(lens < DIST_THRESHOLD) {
				_nVtx = 1;
				_inner = _posB[0];
				return true;
			}
			dir = _vtx[0] * spn::RSqrt(lens) * -1.f;
			_minkowskiSub(dir, 1);

			const Vec2 ori(0);
			Vec2 tmp = _vtx[1] - _vtx[0];
			float tmpLen = tmp.length();
			tmp *= Rcp22Bit(tmpLen);
			float r = tmp.dot(-_vtx[0]);
			if(r > tmpLen)
				return false;
			tmp = _vtx[0] + tmp*r;

			lens = tmp.len_sq();
			if(lens < DIST_THRESHOLD) {
				// ライン上に原点がある
				_inner =_posB[0].l_intp(_posB[1], r * spn::Rcp22Bit(tmpLen));
				_nVtx = 2;
				return true;
			}

			_nVtx = 3;
			LNear res(tmp, LINEPOS::ONLINE);
			int idx = 2;
			for(;;) {
				auto& rV = _vtx[idx];
				// 新たな頂点を追加
				dir = res.first.normalization() * -1.f;
				_minkowskiSub(dir, idx);
				if(dir.dot(rV) < -DOT_THRESHOLD)
					return false;

				auto& ts = ToSearch[idx];
				if(rV.distance(_vtx[ts[0]]) < 1e-7f ||
					rV.distance(_vtx[ts[2]]) < 1e-7f)
				{
					std::cout << _m0 << std::endl << _m1 << std::endl;
					return false;
				}

				// 現時点で三角形が原点を含んでいるか
				bool bIn;
				if((_vtx[1]-_vtx[0]).cw(_vtx[2]-_vtx[0]) > 0) {
					bIn = (_vtx[1]-_vtx[0]).cw(-_vtx[0]) >= 0 &&
					(_vtx[2]-_vtx[1]).cw(-_vtx[1]) >= 0 &&
					(_vtx[0]-_vtx[2]).cw(-_vtx[2]) >= 0;
				} else {
					bIn = (_vtx[1]-_vtx[0]).ccw(-_vtx[0]) >= 0 &&
					(_vtx[2]-_vtx[1]).ccw(-_vtx[1]) >= 0 &&
					(_vtx[0]-_vtx[2]).ccw(-_vtx[2]) >= 0;
				}
// 				float cm = (_vtx[ts[1]] - _vtx[ts[0]]).cw(-_vtx[ts[0]]) *
// 							(_vtx[ts[2]] - _vtx[ts[1]]).cw(-_vtx[ts[1]]);
// 				if(cm >= -1e-5f) {
				if(bIn) {
					_inner = TriangleLerp(_vtx[0], _vtx[1], _vtx[2], ori,
										_posB[0], _posB[1], _posB[2]);
					assert(_m0.isInner(_inner) && _m1.isInner(_inner));
					return true;
				}
				float dist = std::numeric_limits<float>::max();
				idx = -1;
				res = LineCore(_vtx[ts[0]], _vtx[ts[1]]).nearest(ori);
				if(res.second == LINEPOS::ONLINE) {
					float td = res.first.len_sq();
					if(dist > td) {
						dist = td;
						// この辺について調べる
						idx = ts[2];
					}
				}
				res = LineCore(_vtx[ts[1]], _vtx[ts[2]]).nearest(ori);
				if(res.second == LINEPOS::ONLINE) {
					float td = res.first.len_sq();
					if(dist > td) {
						dist = td;
						// この辺について調べる
						idx = ts[0];
					}
				}
				if(dist < 1e-7f) {
					_inner = TriangleLerp(_vtx[0], _vtx[1], _vtx[2], ori,
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
			assert(_szVl <= _vl.size());
			return p;
		}
		const Vec2x2& GEpa::_minkowskiSub(const Vec2& dir, int n) {
			Vec2x2* vp = _allocVert(n);
			vp->second = _m1.support(-dir);
			vp->first = _m0.support(dir) - vp->second;
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
		void GEpa::_printASV() {
			int nV = _szVl;
			for(int i=0 ; i<nV ; i++) {
				PrintV(*_vl[i]);
				std::cout << std::endl;
			}
			for(auto& p : _asv) {
				std::cout << p.dist << ':';
				PrintV(p.dir);
				std::cout << std::hex << p.vtx[0] << ',' << p.vtx[1] << std::endl;
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
			auto res = LineCore(v0.first, v1.first).nearest(c_origin);
			if(res.second == LINEPOS::ONLINE) {
				float len = res.first.length();
				res.first *= Rcp22Bit(len);
				_asv.insert(LmLen{len, res.first, {&v0,&v1}});
				return true;
			}
			return res.second != LINEPOS::NOTHIT;
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
				if(spn::IsNear(d1, fr.dist, 5e-4f)) {
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
				if(spn::IsNear(std::fabs(d1), fr.dist, 5e-4f)) {
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
			assert(false);
		}
		void GEpa::_adjustLoop3() {
			assert(_szVl == 3);
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
				auto res = LineCore(p1.first, p0.first).nearest(c_origin);
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
					res.first *= Rcp22Bit(len);
				_asv.insert(LmLen{len, res.first, {_vl[i], _vl[(i+1)%nV]}});
			}
		}
		void GEpa::_recover2NoHit() {
			auto res = LineCore(_vtx[0], _vtx[1]).nearest(Vec2(0));
			float len = res.first.length();
			LmLen lmlen;
			lmlen.dist = len;
			lmlen.dir *= Rcp22Bit(len);
			if(res.second != LINEPOS::ONLINE) {
				lmlen.vtx[0] = lmlen.vtx[1] = _vl[(res.second == LINEPOS::BEGIN) ? 0 : 1];
			} else {
				lmlen.vtx[0] = _vl[0];
				lmlen.vtx[1] = _vl[1];
			}
			_asv.insert(lmlen);
		}
		void GEpa::_recover2OnHit() {
			LmLen lm;
			lm.dist = 0;
			lm.dir = _vtx[0]-_vtx[1];
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
				res->first = _vtx[i];
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
						const Vec2& v = _minkowskiSub(-_vtx[0].normalization()).first;
						auto res = LineCore(_vtx[0], v).nearest(Vec2(0));
						if(res.second == LINEPOS::ONLINE) {
							float len = res.first.length();
							_asv.insert(LmLen{len, res.first*Rcp22Bit(len), {_vl[0],_vl[1]}});
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
			assert(!getResult());
			assert(!IsNaN(_nvec.first));
			assert(!IsNaN(_nvec.second));
			return _nvec;
		}
		const Vec2& GEpa::getPVector() const {
			assert(getResult());
			assert(!IsNaN(_pvec));
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

		//! 固有のアルゴリズムによる衝突判定
		bool HitCheck(const IModel& mdl0, const IModel& mdl1) {
//			uint32_t id = (mdl0.getCID() << spn::NBits<CTGeo::size>::result) | mdl1.getCID();
			AssertT(Trap, false, (std::domain_error)(const char*), "not implemented yet")
		}
		//! 固有のアルゴリズムでmdlFromのmdlToに対する最深点を算出
		Vec2 HitPos(const IModel& mdlFrom, const IModel& mdlTo) {
			AssertT(Trap, false, (std::domain_error)(const char*), "not implemented yet")
		}
	}
}