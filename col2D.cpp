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
		{spn::COS90, spn::SIN90,
		-spn::SIN90, spn::COS90},
		{spn::COS90, -spn::SIN90,
		spn::SIN90, spn::COS90}
	};
	float Area_x2(const Vec2& v0, const Vec2& v1) {
		return std::fabs(v0.ccw(v1));
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
			Vec2 dir(_m1.center() - _m0.center());
			dir.normalize();

			_minkowskiSub(dir, 0);
			if(dir.dot(_vtx[0]) < -DOT_THRESHOLD) {
				_nVtx = 1;
				return false;
			}
			dir = _vtx[0].normalization() * -1.f;
			_minkowskiSub(dir, 1);

			const Vec2 ori(0);
			Vec2 tmp = _vtx[1] - _vtx[0];
			float tmpLen = tmp.length();
			tmp *= _sseRcp22Bit(tmpLen);
			float r = tmp.dot(-_vtx[0]);
			if(r > tmpLen)
				return false;
			tmp = _vtx[0] + tmp*r;

			float lens = tmp.len_sq();
			if(lens < DIST_THRESHOLD) {
				// ライン上に原点がある
				_inner =_posB[0].l_intp(_posB[1], r * spn::_sseRcp22Bit(tmpLen));
				_nVtx = 2;
				return true;
			}

			_nVtx = 3;
			LineCore::LNear res(tmp, LineCore::ONLINE);
			int idx = 2;
			for(;;) {
				auto& rV = _vtx[idx];
				// 新たな頂点を追加
				dir = res.first.normalization() * -1.f;
				_minkowskiSub(dir, idx);
				if(dir.dot(rV) < -DOT_THRESHOLD)
					return false;

				auto& ts = ToSearch[idx];
				// 現時点で三角形が原点を含んでいるか
				float cm = (_vtx[ts[1]] - _vtx[ts[0]]).ccw(-_vtx[ts[0]]) *
							(_vtx[ts[2]] - _vtx[ts[1]]).ccw(-_vtx[ts[1]]);
				if(cm >= -1e-5f) {
					_inner = TriangleLerp(_vtx[0], _vtx[1], _vtx[2], ori,
										_posB[0], _posB[1], _posB[2]);
					return true;
				}
				res = LineCore(_vtx[ts[0]], _vtx[ts[1]]).nearest(rV);
				if(res.second == LineCore::ONLINE) {
					// この辺について調べる
					idx = ts[2];
					continue;
				}
				res = LineCore(_vtx[ts[1]], _vtx[ts[2]]).nearest(rV);
				if(res.second == LineCore::ONLINE) {
					// この辺について調べる
					idx = ts[0];
					continue;
				}
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

		// ---------------------- GEpa ----------------------
		const Vec2x2& GEpa::_minkowskiSub(const Vec2& dir, int n) {
			int nV = _vl.size();
			auto* vp = new Vec2x2;
			vp->second = _m1.support(-dir);
			vp->first = _m0.support(dir) - vp->second;
			if(n >= 0)
				_vl.insert(_vl.begin()+n, vp);
			else
				_vl.push_back(vp);

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
			for(auto& p : _vl) {
				PrintV(*p);
				std::cout << std::endl;
			}
			for(auto& p : _asv) {
				std::cout << p.dist << ':';
				PrintV(p.dir);
				std::cout << std::hex << p.vtx[0] << ',' << p.vtx[1] << std::endl;
			}
		}

		int GEpa::_getIndex(const Vec2x2 *vp) const {
			int nV = _vl.size();
			for(int i=0 ; i<nV ; i++) {
				if(_vl[i] == vp)
					return i;
			}
			return 0xffff;
		}
		namespace {
			const Vec2 c_origin(0);
		}
		void GEpa::_addAsv(const Vec2& v0, const Vec2& v1, const Vec2x2* (&vtx)[2]) {
			auto res = LineCore(v0,v1).nearest(c_origin);
			float len = res.first.length();
			if(res.second == LineCore::ONLINE) {
				res.first *= _sseRcp22Bit(len);
				_asv.insert(LmLen{len, res.first, {vtx[0], vtx[1]}});
			}
		}
		void GEpa::_epaMethodOnHit() {
			_adjustLoop();
			_geneASV();

			// 衝突時: 脱出ベクトルを求める
			float minLen = std::numeric_limits<float>::max();
			while(!_asv.empty()) {
				auto fr = _asv.pop_frontR();

				_printASV();
				const Vec2x2 &v1 = _minkowskiSub(fr.dir, _getIndex(fr.vtx[0])),
							&v0 = *fr.vtx[0],
							&v2 = *fr.vtx[1];
				_printASV();
				float d1 = fr.dir.dot(v1.first);
				if(spn::IsNear(d1, fr.dist, 5e-4f)) {
					if(minLen > fr.dist) {
						minLen = fr.dist;
						_pvec = fr.dir * fr.dist;
					}
					if(_asv.empty() || _asv.front().dist > fr.dist)
						break;
				}

				_addAsv(v0.first, v1.first, fr.vtx);
				_addAsv(v1.first, v2.first, fr.vtx);
			}
		}
		void GEpa::_epaMethodNoHit() {
			_adjustLoop();
			_geneASV();
			// 無衝突時: 最近傍対を求める
			while(!_asv.empty()) {
				auto fr = _asv.pop_frontR();
				if(!fr.vtx[1]) {
					auto& p = *fr.vtx[0];
					_nvec.second = p.second;
					_nvec.first = p.second + fr.dir * fr.dist;
					return;
				}
				_printASV();
				const Vec2x2& v1 = _minkowskiSub(-fr.dir, _getIndex(fr.vtx[0])),
							&v0 = *fr.vtx[0],
							&v2 = *fr.vtx[1];
				_printASV();
				float d1 = fr.dir.dot(v1.first);
				if(spn::IsNear(d1, fr.dist, 5e-4f)) {
					_nvec.second = v1.second;
					_nvec.first = v1.second + fr.dir * fr.dist;
					return;
				}
				_addAsv(v0.first, v1.first, fr.vtx);
				_addAsv(v1.first, v2.first, fr.vtx);
			}
			assert(false);
		}
		void GEpa::_adjustLoop() {
			int nV = _vl.size();
			if(nV < 3)
				return;
			// 頂点0を基準にCCWの値でソート
			struct Ccw {
				float	val;
				int		idx;

				bool operator > (const Ccw& c) const {
					return val > c.val;
				}
			};
			spn::AssocVec<Ccw, std::greater<Ccw>> asv;

			const Vec2& v0 = _vl[0]->first;
			for(int i=1 ; i<nV ; i++)
				asv.insert(Ccw{v0.ccw(_vl[i]->first), i});

			decltype(_vl) nvl(nV);
			auto* pDst = &nvl[0];
			*pDst++ = _vl[0];
			while(!asv.empty()) {
				auto c = asv.pop_frontR();
				*pDst++ = _vl[c.idx];
			}
			std::swap(_vl, nvl);
		}
		void GEpa::_geneASV() {
			int nV = _vl.size();
			assert(nV >= 3 || !getResult());
			_asv.clear();

			Vec2 ori(0);
			for(int i=0 ; i<nV ; i++) {
				auto &p0 = *_vl[i],
					&p1 = *_vl[(i+1)%nV];
				auto res = LineCore(p1.first, p0.first).nearest(ori);
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
					res.first *= _sseRcp22Bit(len);
				_asv.insert(LmLen{len, res.first, {_vl[i], _vl[(i+1)%nV]}});
			}
		}

		void GEpa::_recover2NoHit() {
			auto res = LineCore(_vtx[0], _vtx[1]).nearest(Vec2(0));
			float len = res.first.length();
			LmLen lmlen;
			lmlen.dist = len;
			lmlen.dir *= _sseRcp22Bit(len);
			if(res.second != LineCore::ONLINE) {
				lmlen.vtx[0] = lmlen.vtx[1] = _vl[(res.second == LineCore::BEGIN) ? 0 : 1];
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

		GEpa::GEpa(const IModel& m0, const IModel& m1): GSimplex(m0,m1) {
			int nV = _nVtx;
			_vl.resize(nV);
			for(int i=0 ; i<nV ; i++)
				_vl[i] = new Vec2x2(_vtx[i], _posB[i]);

			do {
				if(getResult()) {
					if(nV == 1) {
						_pvec *= 0;
						break;
					} else if(nV == 2)
						_recover2OnHit();
					_epaMethodOnHit();
				} else {
					if(nV == 1) {
						const Vec2& v = _minkowskiSub(-_vtx[0].normalization()).first;
						auto res = LineCore(_vtx[0], v).nearest(Vec2(0));
						if(res.second == LineCore::ONLINE) {
							float len = res.first.length();
							_asv.insert(LmLen{len, res.first*_sseRcp22Bit(len), {_vl[0],_vl[1]}});
						}
						else {
							const auto& p = *_vl[0];
							_nvec.second = p.second;
							_nvec.first = _nvec.second + p.first;
							break;
						}
					} else if(nV == 2)
						_recover2NoHit();
					_epaMethodNoHit();
				}
			} while(false);
		}
		Vec2x2 GEpa::getNearestPair() const {
			return _nvec;
		}
		const Vec2& GEpa::getPVector() const {
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
			__m128 xm[2] = {LOADPS_Z3(ar[0]), LOADPS_Z3(ar[1])};
			xm[0] = spn::_mmDivPs(xm[0], _mm_load1_ps(&toV1.x));
			xm[1] = _mm_add_ps(_mm_mul_ps(xm[0], _mm_load1_ps(&toV1.y)), xm[1]);
			xm[1] = spn::_mmDivPs(xm[1], _mm_shuffle_ps(xm[1], xm[1], _MM_SHUFFLE(1,1,1,1)));
			xm[0] = _mm_sub_ps(xm[0], _mm_mul_ps(xm[1], _mm_shuffle_ps(xm[0], xm[0], _MM_SHUFFLE(1,1,1,1))));

			Float2 ret;
			_mm_store_ss(&ret.first, xm[0]);
			_mm_store_ss(&ret.second, xm[1]);
			return ret;
		}

		//! 固有のアルゴリズムによる衝突判定
		bool HitCheck(const IModel& mdl0, const IModel& mdl1) {
			uint32_t id = (mdl0.getCID() << spn::NBits<CTGeo::size>::result) | mdl1.getCID();
			throw std::runtime_error("not implemented yet");
		}
		//! 固有のアルゴリズムでmdlFromのmdlToに対する最深点を算出
		Vec2 HitPos(const IModel& mdlFrom, const IModel& mdlTo) {
			throw std::runtime_error("not implemented yet");
		}
	}
}