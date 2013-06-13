#include "geom2D.hpp"
#include "spinner/bits.hpp"
#include "spinner/common.hpp"
#include "spinner/spn_math.hpp"
#include "spinner/misc.hpp"
#include "spinner/assoc.hpp"
#include <limits>

namespace boom {
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
			if(dir.dot(_vtx[0]) < -DOT_THRESHOLD)
				return false;
			dir = _vtx[0].normalization() * -1.f;
			_minkowskiSub(dir, 1);

			const Vec2 ori(0);
			Vec2 tmp = _vtx[1] - _vtx[0];
			float tmpLen = tmp.length();
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
				float cm = (_vtx[ts[1]] - _vtx[ts[0]]).ccw(-_vtx[ts[0]]) * (_vtx[ts[1]] - _vtx[ts[2]]).ccw(-_vtx[ts[2]]);
				if(cm < 0) {
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
		const Vec2& GEpa::_minkowskiSub(const Vec2& dir, int n) {
			Vec2x2 tmp;
			tmp.second = _m1.support(-dir);
			tmp.first = _m0.support(dir) - tmp.second;

			if(n < 0) {
				_vl.push_back(tmp);
				return _vl.back().first;
			}
			_vl.insert(_vl.begin()+n, tmp);
			return _vl[n].first;
		}
		void GEpa::_epaMethod() {
			if(_vl.size() == 2) {
				// 補助法線でちゃんとした凸包を形成
				Vec2 toV1(_vl[1].first - _vl[0].first);
				Vec2 dir = toV1 * cs_mRot[0];
				dir.normalize();

				const Vec2& v = _minkowskiSub(dir);
				LineCore lc(_vl[0].first, _vl[1].first);
				float r = lc.ratio(v);
				if(spn::IsInRange(r, 0.f, 1.f)) {
					_pvec = Vec2(0);
					return;
				}

				dir = toV1 * cs_mRot[1];
				const Vec2& v2 = _minkowskiSub(dir, 1);
				if(spn::IsInRange(lc.ratio(v2), 0.f, 1.f)) {
					_pvec = Vec2(0);
					return;
				}
			}

			struct LmLen {
				float	dist;
				Vec2	dir;
				int		idx;

				bool operator < (const LmLen& len) const {
					return dist < len.dist;
				}
			};
			spn::AssocVec<LmLen> asv;
			const Vec2 ori(0);
			float minDist = std::numeric_limits<float>::max();
			Vec2 minVec, minPosB;
			// 辺の最短リストを構築
			auto addLineDist = [ori, &asv, &minDist, &minVec, &minPosB](const Vec2& vt0, const Vec2& vt1, int idx) {
				auto ls = LineCore(vt0,vt1).nearest(ori);
				float len = ls.first.length();
				if(minDist > len) {
					minDist = len;
					minVec = ls.first;
				}
				if(ls.second == LineCore::ONLINE)
					asv.insert(LmLen{len, ls.first * _sseRcp22Bit(len), idx});
			};

			// 現在の頂点情報で距離リストを構築
			int nV = _vl.size();
			for(int i=0 ; i<nV ; i++)
				addLineDist(_vl[i].first, _vl[(i+1)%nV].first, i);

			if(!getResult()) {
				// 無衝突時: 最近傍対を求める
				while(!asv.empty()) {
					auto fr = asv.pop_frontR();
					_minkowskiSub(-fr.dir, fr.idx);
					const auto &v0 = _vl[fr.idx].first,
								&v1 = _vl[fr.idx+1].first,
								&v2 = _vl[(fr.idx+2)%nV].first;
					float ccw = (v1 - v0).ccw(v2 - v0);
					if(ccw > -1e-4f)
						break;
					addLineDist(v0,v1, fr.idx);
					addLineDist(v1,v2, fr.idx+1);
				}
				_nvec[0] = minPosB + minVec;
				_nvec[1] = minPosB;
				return;
			} else {
				// 衝突時: 脱出ベクトルを求める
				while(!asv.empty()) {
					auto fr = asv.pop_frontR();
					_minkowskiSub(fr.dir, fr.idx);
					// ライン上なら終了
					const auto &v0 = _vl[fr.idx].first,
								&v1 = _vl[fr.idx+1].first,
								&v2 = _vl[(fr.idx+2)%nV].first;
					Vec2 v01(v1-v0),
						v02(v2-v0);
					v01.normalize();
					v02.normalize();
					if(v01.dot(v02) > 1.f-1e-3f) {
						_pvec = fr.dir * fr.dist;
						return;
					}
					// 距離リストに登録
					addLineDist(v0,v1, fr.idx);
					addLineDist(v1,v2, fr.idx+1);
				}
				// リストが空になる事はない筈
				throw std::runtime_error("something wrong");
			}
		}
		void GEpa::_adjustLoop() {
			int nV = _vl.size();
			if(nV > 2) {
				// 頂点0を基準にCCWの値でソート
				struct Ccw {
					float	val;
					int		idx;

					bool operator < (const Ccw& c) const {
						return val < c.val;
					}
				};
				spn::AssocVec<Ccw> asv;

				const Vec2& v0 = _vl[0].first;
				for(int i=1 ; i<nV ; i++)
					asv.insert(Ccw{v0.ccw(_vl[i].first), i});

				decltype(_vl) nvl(nV);
				auto* pDst = &nvl[0];
				*pDst++ = _vl[0];
				while(!asv.empty()) {
					auto c = asv.pop_frontR();
					*pDst++ = _vl[c.idx];
				}
				std::swap(_vl, nvl);
			}
		}

		GEpa::GEpa(const IModel& m0, const IModel& m1): GSimplex(m0,m1) {
			int nV = _nVtx;
			for(int i=0 ; i<nV ; i++)
				_vl[i] = Vec2x2(_vtx[i], _posB[i]);
			_adjustLoop();
			_epaMethod();
		}
		Vec2x2 GEpa::getNearestPair() const {
			return Vec2x2(_nvec[0], _nvec[1]);
		}
		const Vec2& GEpa::getPVector() const {
			return _pvec;
		}
		const spn::AMat22 GEpa::cs_mRot[2] = {
			{spn::COS45, spn::SIN45,
			-spn::SIN45, spn::COS45},
			{spn::COS45, -spn::SIN45,
			spn::SIN45, spn::COS45}
		};

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