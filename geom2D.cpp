#include "geom2D.hpp"
#include "spinner/bits.hpp"
#include "spinner/misc.hpp"

namespace boom {
	namespace geo2d {
		using spn::Bit;
		template <class CLIP>
		inline Vec2 NearestPoint(const StLineCore& st, const Vec2& p, CLIP clip) {
			Vec2 toP = p - st.pos;
			float d = st.dir.dot(toP);
			return st.pos + st.dir * clip(d);
		}
		template <class CLIP>
		inline Vec2x2 NearestPoint(const StLineCore& st0, const StLineCore& st1, CLIP clip0, CLIP clip1) {
			float st0d = st0.dir.len_sq(),
					st1d = st1.dir.len_sq(),
					st01d = st0.dir.dot(st1.dir);
			float d = st0d * st1d - spn::Square(st01d);
			if(std::fabs(d) < 1e-5f) {
				// 2つの直線は平行
				return Vec2x2(st0.pos, NearestPoint(st1, st0.pos, [](float f){return f;}));
			}
			spn::Mat22 m0(st1d, st01d,
						st01d, st0d);
			spn::Vec2	m1((st1.pos - st0.pos).dot(st0.dir),
							(st0.pos - st1.pos).dot(st1.dir));

			m1 = m0 * m1;
			m1 *= _sseRcp22Bit(d);
			return Vec2x2(st0.pos + st0.dir * clip0(m1.x),
							st1.pos + st1.dir * clip1(m1.y));
		}

		RForce::F& RForce::F::operator += (const F& f) {
			linear += f.linear;
			torque += f.torque;
			return *this;
		}
		RForce& RForce::operator +=(const RForce& rf) {
			sdump += rf.sdump;
			fricD += rf.fricD;
			return *this;
		}

		CircleCore IModel::getBBCircle() const {
			throw std::runtime_error("not supported function");
		}

		// ---------------------- Point ----------------------
		float PointCore::distance(const LineCore& l) const {
			Vec2 dir(l.point[1] - l.point[0]);
			dir.normalize();
			return std::fabs(((*this) - l.point[0]).ccw(dir));
		}
		LNear PointCore::nearest(const LineCore& l) const {
			return l.nearest(*this);
		}
		bool PointCore::hit(const PointCore& p) const {
			return distance(p) <= NEAR_THRESHOLD;
		}

		Vec2 Point::support(const Vec2& dir) const {
			return *this;
		}
		Vec2 Point::center() const {
			return *this;
		}
		CircleCore Point::getBBCircle() const {
			// 円の半径が0だと点同士の時にヒットしないので微量含める
			return CircleCore(*this, 1e-6f);
		}

		// ---------------------- Circle ----------------------
		CircleCore::CircleCore(const Vec2& c, float r): center(c), radius(r) {}
		float CircleCore::area() const {
			throw std::runtime_error("");
		}
		Vec2 CircleCore::support(const Vec2& dir) const {
			return dir * radius + center;
		}
		bool CircleCore::hit(const Vec2& pt) const {
			return center.dist_sq(pt) <= spn::Square(radius);
		}

		Vec2 Circle::support(const Vec2& dir) const { return CircleCore::support(dir); }
		Vec2 Circle::center() const { return CircleCore::center; }
		bool Circle::isInner(const Vec2& pos) const { return CircleCore::hit(pos); }
		CircleCore Circle::getBBCircle() const { return *this; }

		CircleModel::CircleModel(): _rflag(0xff) {}
		CircleModel::CircleModel(const CircleCore& c): _circle(c), _rflag(0xff) {}

		Vec2 CircleModel::support(const Vec2& dir) const { return _circle.support(dir); }
		Vec2 CircleModel::center() const { return _circle.center; }
		bool CircleModel::isInner(const Vec2& pos) const { return _circle.hit(pos); }
		CircleCore CircleModel::getBBCircle() const { return _circle; }

		// ---------------------- StLine ----------------------
		StLineCore::StLineCore(const Vec2& p, const Vec2& d): pos(p), dir(d) {}
		Vec2x2 StLineCore::nearest(const StLineCore& st) const {
			auto fn = [](float f) { return f; };
			return NearestPoint(*this, st, fn, fn);
		}
		Vec2 StLineCore::nearest(const Vec2& p) const {
			return pos + (p-pos).dot(dir);
		}
		float StLineCore::distance(const Vec2& p) const {
			return std::fabs(dir.ccw(p - pos));
		}
		Vec2 StLineCore::placeOnLine(const Vec2& p) const {
			return pos + dir*posDot(p);
		}
		float StLineCore::posDot(const Vec2& p) const {
			Vec2 tv(p-pos);
			return dir.dot(tv);
		}

		// ---------------------- Ray ----------------------
		RayCore::RayCore(const Vec2& p, const Vec2& d): pos(p), dir(d) {}
		const StLineCore& RayCore::asStLine() const { return *reinterpret_cast<const StLineCore*>(this); }
		Vec2x2 RayCore::nearest(const RayCore& r) const {
			auto fn = [](float f) { return std::max(0.f, f); };
			return NearestPoint(asStLine(), r.asStLine(), fn, fn);
		}
		Vec2 RayCore::nearest(const Vec2& p) const {
			return NearestPoint(asStLine(), p, [](float f){ return std::max(0.f,f); });
		}

		// ---------------------- Line ----------------------
		LineCore::LineCore(const Vec2& v0, const Vec2& v1): point{v0,v1} {}
		float LineCore::distance(const LineCore& l) const {
			auto fn = [](float f) { return spn::Saturate(f, 0.f, 1.f); };
			Vec2x2 vp = NearestPoint(toStLine(), l.toStLine(), fn,fn);
			return vp.first.distance(vp.second);
		}
		float LineCore::length() const {
			return point[0].distance(point[1]);
		}
		float LineCore::len_sq() const {
			return point[0].dist_sq(point[1]);
		}
		bool LineCore::hit(const LineCore& l) const {
			return distance(l) < 1e-5f;
		}
		LNear LineCore::nearest(const Vec2& p) const {
			Vec2 toP(p-point[0]),
				toV1(point[1]-point[0]);
			float lenV1 = toV1.length();
			toV1 *= _sseRcp22Bit(lenV1);
			float d = toV1.dot(toP);
			if(d <= 0)
				return LNear(point[0], LINEPOS::BEGIN);
			else if(d >= lenV1)
				return LNear(point[1], LINEPOS::END);
			else
				return LNear(point[0]+toV1*d, LINEPOS::ONLINE);
		}
		float LineCore::ratio(const Vec2& p) const {
			Vec2 toP(p-point[0]),
				toV1(point[1]-point[0]);
			float len = toV1.length();
			return toV1.dot(toP) * _sseRcp22Bit(len);
		}
		StLineCore LineCore::toStLine() const {
			return StLineCore(point[0], (point[1]-point[0]).normalization());
		}
		Vec2 LineCore::support(const Vec2& dir) const {
			float d[2] = {dir.dot(point[0]), dir.dot(point[1])};
			if(d[0] > d[1])
				return point[0];
			return point[1];
		}
		bool LineCore::online(const Vec2& p) const {
			Vec2 toV1(point[1]-point[0]),
				toP(p-point[0]);
			toV1.normalize();
			return spn::IsNear(toV1.dot(toP), toP.length(), 1e-5f);
		}
		LNear LineCore::crossPoint(const LineCore& l) const {
			auto fn = [](float f){ return spn::Saturate(f, 0.f,1.f); };
			Vec2 cp = NearestPoint(l.toStLine(), point[0], fn);
			return LNear(cp, online(cp) ? LINEPOS::ONLINE : LINEPOS::NOTHIT);
		}
		LNear LineCore::crossPoint(const StLineCore& l) const {
			float c0 = l.dir.ccw(point[0]-l.pos),
				c1 = l.dir.ccw(point[1]-l.pos);
			if(c0*c1 <= 0) {
				Vec2 dir(point[1]-point[0]);
				dir.normalize();
				float d = c0 / (c0+c1);
				return LNear(point[0] + dir*d, LINEPOS::ONLINE);
			}
			return LNear(Vec2(), LINEPOS::NOTHIT);
		}

		Vec2 Line::support(const Vec2& dir) const {
			return LineCore::support(dir);
		}
		Vec2 Line::center() const {
			return (point[0] + point[1]) * 0.5f;
		}
		CircleCore LineCore::getBBCircle() const {
			return CircleCore(point[0] + point[1] * 0.5f,
							point[0].distance(point[1]));
		}
		CircleCore Line::getBBCircle() const { return LineCore::getBBCircle(); }

		// ---------------------- Poly ----------------------
		PolyCore::PolyCore(const Vec2& p0, const Vec2& p1, const Vec2& p2): point{p0,p1,p2} {}
		float PolyCore::area() const {
			return CalcArea(point[0], point[1], point[2]);
		}
		Vec2 PolyCore::center() const {
			return (point[0] + point[1] + point[2]) * (1.0f/3);
		}
		float PolyCore::CalcArea(const Vec2& p0, const Vec2& p1, const Vec2& p2) {
			return  (p1-p0).cw(p2-p0) * 0.5f;
		}
		float PolyCore::CalcArea(const Vec2& p0, const Vec2& p1) {
			return p0.cw(p1) * 0.5f;
		}
		Vec2 PolyCore::support(const Vec2& dir) const {
			throw std::runtime_error("");
		}
		void PolyCore::addOffset(const Vec2& ofs) {
			for(int i=0 ; i<3 ; i++)
				point[i] += ofs;
		}
		bool PolyCore::isInTriangle(const Vec2& p) const {
			Vec2 vt(p-point[0]),
				v1(point[1]-point[0]),
				v2(point[2]-point[0]);
			return v1.cw(vt) >= 0 &&
					vt.cw(v2) >= 0;
		}
		int PolyCore::getObtuseCorner() const {
			AVec2 v01(point[1]-point[0]),
				v02(point[2]-point[0]),
				v12(point[2]-point[1]);
			if(v01.dot(v02) < 0)
				return 0;

			v01 *= -1;
			if(v01.dot(v12) < 0)
				return 1;

			v12 *= -1;
			v02 *= -1;
			if(v02.dot(v12) < 0)
				return 2;
			return -1;
		}
		CircleCore PolyCore::getBBCircle() const {
			int id = getObtuseCorner();
			if(id >= 0) {
				// 鈍角を持っていれば直径を使う
				const Vec2 &v0 = point[id],
							&v1 = point[spn::CndSub(id+1, 3)];
				return CircleCore((v0+v1)*0.5f, v0.distance(v1));
			} else {
				// なければ3点の外接円
				StLineCore line0((point[1]+point[0]) * 0.5f,
								(point[1]-point[0]) * cs_mRot90[0]),
							line1((point[2]+point[0]) * 0.5f,
								(point[2]-point[0]) * cs_mRot90[0]);
				auto vp = line0.nearest(line1);
				return CircleCore(vp.first,
									vp.first.distance(point[0]));
			}
		}

		Vec2 Poly::support(const Vec2& dir) const { return PolyCore::support(dir); }
		Vec2 Poly::center() const { return PolyCore::center(); }
		bool Poly::isInner(const Vec2& pos) const { return isInTriangle(pos); }
		CircleCore Poly::getBBCircle() const { return PolyCore::getBBCircle(); }

		PolyModel::PolyModel(): _rflag(RFL_ALL) {}
		PolyModel::PolyModel(const Vec2& p0, const Vec2& p1, const Vec2& p2): _poly{p0,p1,p2}, _rflag(RFL_ALL) {}
		float PolyModel::getArea() const {
			if(Bit::ChClear(_rflag, RFL_AREA))
				_area = _poly.area();
			return _area;
		}
		const AVec2& PolyModel::getCenter() const {
			if(Bit::ChClear(_rflag, RFL_CENTER))
				_center = _poly.center();
			return _center;
		}
		void PolyModel::setPoint(int n, const Vec2& v) {
			_poly.point[n] = v;
			_rflag = RFL_ALL;
		}
		void PolyModel::addOffset(const Vec2& ofs) {
			_poly.addOffset(ofs);
		}
		float PolyModel::getInertia() const {
			if(Bit::ChClear(_rflag, RFL_INERTIA)) {
				_inertia = (1.0f/18) * (_poly.point[0].dot(_poly.point[0])
										+ _poly.point[0].dot(_poly.point[0])
										+ _poly.point[0].dot(_poly.point[0])
										- _poly.point[1].dot(_poly.point[2])
										- _poly.point[2].dot(_poly.point[0])
										- _poly.point[0].dot(_poly.point[1]));
			}
			return _inertia;
		}
		Vec2 PolyModel::support(const Vec2& dir) const {
			return _poly.support(dir);
		}
		Vec2 PolyModel::center() const {
			return Vec2(getCenter());
		}
		bool PolyModel::isInner(const Vec2& pos) const { return _poly.isInTriangle(pos); }
		CircleCore PolyModel::getBBCircle() const { return _poly.getBBCircle(); }

		// ---------------------- Convex ----------------------
		ConvexCore::ConvexCore(const PointL& pl): point(pl) {}
		ConvexCore::ConvexCore(PointL&& pl): point(pl) {}
		ConvexCore::ConvexCore(std::initializer_list<Vec2> v): point(v.size()) {
			auto itrD = point.begin();
			auto itr = v.begin();
			while(itr != v.end())
				*itrD++ = *itr++;
		}
		ConvexCore ConvexCore::FromConcave(const PointL& src) {
			int nV = src.size();
			assert(nV >= 3);

			// X軸についてソート
			PointL tsrc(src);
			std::sort(tsrc.begin(), tsrc.end(), [](const Vec2& v0, const Vec2& v1){ return v0.x < v1.x; });

			PointL pts(nV*2);
			Vec2* pDst = &pts[0];
			*pDst++ = tsrc[0];
			*pDst++ = tsrc[1];
			for(int rc=2 ; rc<nV ; rc++) {
				if(Vec2::Cw(tsrc[rc-2], tsrc[rc-1], tsrc[rc]) < 0)
					--pDst;
				*pDst++ = tsrc[rc];
			}
			*pDst++ = tsrc[nV-1];
			*pDst++ = tsrc[nV-2];
			for(int rc=nV-3 ; rc>=0 ; rc--) {
				if(Vec2::Cw(tsrc[rc+2], tsrc[rc+1], tsrc[rc]) < 0)
					--pDst;
				*pDst++ = tsrc[rc];
			}
			assert(&pts[0]+nV*2 <= pDst);
			pts.resize(pDst - &pts[0]);
			return ConvexCore(std::move(pts));
		}

		float ConvexCore::area() const {
			AreaSum as;
			iterate(as);
			return as.result;
		}
		void ConvexCore::addOffset(const Vec2& ofs) {
			for(auto& p : point)
				p += ofs;
		}
		Vec2 ConvexCore::support(const Vec2& dir) const {
			Vec2 result = point[0];
			float dMax = point[0].dot(dir);
			int nV = point.size();
			for(int i=1 ; i<nV ; i++) {
				float d = point[i].dot(dir);
				if(d > dMax) {
					result = point[i];
					dMax = d;
				}
			}
			return result;
		}
		Vec2 ConvexCore::center() const {
			return std::get<2>(area_inertia_center());
		}
		LineCore ConvexCore::getOuterLine(int n) const {
			return LineCore(point[(n+1)%point.size()], point[n]);
		}
		Convex Convex::getOverlappingConvex(const Convex& cnv, const Vec2& inner) const {
			PointL pt0(getOverlappingPoints(cnv, inner)),			// m0がm1にめり込んでいる頂点リストを出力
					pt1(cnv.getOverlappingPoints(*this, inner));	// m1がm0にめり込んでいる頂点リストを出力
			int nV0 = pt0.size(),
				nV1 = pt1.size();
			if(nV0 == 0) {
				assert(nV1 >= 3);
				// m1のめり込んだ頂点全部がm0の1つの辺に収まっている
				return Convex(std::move(pt1));
			} else if(nV1 == 0) {
				assert(nV0 >= 3);
				// m0のめり込んだ頂点全部がm1の1つの辺に収まっている
				return Convex(std::move(pt0));
			} else {
				assert(nV0>=3 && nV1>=3);
				PointL pt(nV0 + nV1 - 2);
				auto* pDst = &pt[0];
				for(int i=0 ; i<nV0-1 ; i++)
					*pDst++ = pt0[i];
				// 繋ぎ目の処理: m0リストの始点をRayに変換した物とm1の終端で頂点を1つ
				auto res = LineCore(pt0[nV0-2], pt0[nV0-1]).crossPoint(LineCore(pt1[0], pt1[1]));
				assert(res.second == LINEPOS::ONLINE);
				*pDst++ = res.first;
				for(int i=0 ; i<nV1-1 ; i++)
					*pDst++ = pt1[i];
				// m1リストの始点をRayに変換した物とm0の終端で頂点を1つ
				res = LineCore(pt1[nV1-2], pt1[nV1-1]).crossPoint(LineCore(pt0[0], pt0[1]));
				assert(res.second == LINEPOS::ONLINE);
				*pDst = res.first;
				return Convex(std::move(pt));
			}
		}

		PointL ConvexCore::getOverlappingPoints(const IModel& mdl, const Vec2& inner) const {
			auto res = checkPosition(inner);
			if(res.first != POSITION::OUTER) {
				int nV = point.size();
				auto fn = std::bind((int (*)(int,int))spn::CndSub, std::placeholders::_1, nV);
				auto binSearch = [this, &fn, &mdl](int a, int b) -> int {
					for(;;) {
						if(a+1 == b)
							break;
						int c = (a+b)/2;
						if(mdl.isInner(point[fn(c)]))
							b = c;
						else
							a = c;
					}
					return fn(a);
				};

				// 2分探索で衝突が始まる地点を探す
				int a = res.second,
					b = res.second + nV;
				int begI = binSearch(a, b);		// 衝突開始インデックス

				// 衝突が終わる地点を探す
				a = res.second;
				b = begI;
				int endI = binSearch(a,b);

				PointL pts(spn::CndRange(endI - begI, nV));
				auto* pDst = &pts[0];
				// 開始地点から終了地点までを出力
				while(begI != endI) {
					*pDst++ = point[begI];
					begI = fn(begI+1);
				}
				return std::move(pts);
			}
			return PointL();
		}

		ConvexCore::CPos ConvexCore::checkPosition(const Vec2& pos) const {
			// 適当に内部点を算出
			Vec2 inner = (point[0] + point[1] + point[2]) * (1.f/3);
			Vec2 toP(pos - inner);
			if(toP.len_sq() < 1e-6f) {
				// 重心がちょうどposと重なってしまったら少しずらす
				inner.lerp(point[0], 0.5f);
			}

			// 内部のどの三角形に該当するか2分探索
			int nV = point.size(),
				a = 0,
				b = nV;
			while(a+1 < b) {
				int c = (b-a) / 2;
				float dA = (point[a]-inner).cw(toP),
					dB = (point[(b%nV)]-inner).cw(toP),
					dC = (point[c]-inner).cw(toP);
				if(dA >= 0 && dB < 0)
					b = c;
				else if(dB>=0 && dC < 0)
					a = c;
			}
			float d = (point[b] - point[a]).cw(pos - point[a]);
			POSITION posf;
			if(d > 0)
				posf = POSITION::INNER;
			else if(d < 0)
				posf = POSITION::OUTER;
			else
				posf = POSITION::ONLINE;
			return CPos(posf, a);
		}
		bool ConvexCore::isInner(const Vec2& pos) const {
			return checkPosition(pos).first != POSITION::OUTER;
		}
		std::tuple<float,float,Vec2> ConvexCore::area_inertia_center() const {
			int nL = point.size();
			AreaList al(nL);
			iterate(std::ref(al));
			float invarea = _sseRcp22Bit(al.sum);
			float areaInv3 = invarea * (1.0f/3),
					areaInv6 = invarea * (1.0f/6);
			Vec2 center(0,0);
			float inertia = 0;
			iterate([&, areaInv3, areaInv6](int n, const Vec2& p0, const Vec2& p1) {
				center += (p0 + p1) * al.areaL[n] * areaInv3;
				inertia += al.areaL[n] * areaInv6 * (p0.dot(p0) + p0.dot(p1) + p1.dot(p1));
			});
			inertia -= center.len_sq();
			return std::make_tuple(al.sum, inertia, center);
		}
		std::pair<ConvexCore,ConvexCore> ConvexCore::splitTwo(const StLineCore& l) const {
			int nV = point.size();
			PointL pt0(nV+1), pt1(nV+1);		// 最初に最大容量確保しておき、後で縮める
			auto *pDst0 = &pt0[0],
				*pDst1 = &pt1[0];

			constexpr float DOT_THRESHOLD = 1e-5f;
			// ライン上は前回値, プラスは1, マイナスは0を返す
			auto fcheck = [](int prev, float d) -> int{
				if(d < -DOT_THRESHOLD)
					return 0;
				if(d < DOT_THRESHOLD)
					return prev;
				return 1;
			};
			auto fadd = [&pDst0, &pDst1, &l](const Vec2& pPrev, const Vec2& pCur, int flg) {
				switch(flg) {
					case 0x03:		// Left -> Left
						*pDst0++ = pCur;
						break;
					case 0x02: {	// Left -> Right
						auto res = LineCore(pPrev, pCur).crossPoint(l);
						assert(res.second == LINEPOS::ONLINE);
						*pDst0++ = res.first;
						*pDst1++ = pCur;
						break; }
					case 0x01: {	// Right -> Left
						auto res = LineCore(pPrev, pCur).crossPoint(l);
						assert(res.second == LINEPOS::ONLINE);
						*pDst1++ = res.first;
						*pDst0++ = pCur;
						break; }
					case 0x00:		// Right -> Right
						*pDst1++ = pCur;
						break;
				}
			};

			int fchk = fcheck(0, l.dir.ccw(point[0]));
			int prev = fchk;
			int flag = 0;
			for(int i=1 ; i<nV ; i++) {
				// 正数が左側、負数は右側
				const auto& p = point[i];
				prev = fcheck(prev, l.dir.ccw(p));
				flag = ((flag<<1) | prev) & 0x03;

				fadd(point[i-1], p, flag);
			}
			flag = ((flag<<1) | fchk) & 0x03;
			fadd(point[nV-1], point[0], flag);

			pt0.resize(pDst0 - &pt0[0]);
			pt1.resize(pDst1 - &pt1[0]);
			return std::make_pair(ConvexCore(std::move(pt0)),
								  ConvexCore(std::move(pt1)));
		}
		ConvexCore ConvexCore::split(const StLineCore& l) {
			auto res = splitTwo(l);
			std::swap(res.first, *this);
			return std::move(res.second);
		}
		void ConvexCore::splitThis(const StLineCore& l) {
			split(l);
		}
		CircleCore ConvexCore::getBBCircle() const {
			// 多分遅いアルゴリズムだが、今はこれで我慢
			// 全ての3点の組み合わせを調べる
			int nV = point.size();
			assert(nV >= 3);
			CircleCore c;
			c.radius = -1;
			for(int i=0 ; i<nV-2 ; i++) {
				for(int j=i+1 ; j<nV-1 ; j++) {
					for(int k=j+1 ; k<nV ; k++) {
						PolyCore p(point[i], point[j], point[k]);
						auto tc = p.getBBCircle();
						if(c.radius < tc.radius)
							c = tc;
					}
				}
			}
			return c;
		}

		Vec2 Convex::support(const Vec2& dir) const { return ConvexCore::support(dir); }
		Vec2 Convex::center() const { return ConvexCore::center(); }
		bool Convex::isInner(const Vec2& pos) const { return ConvexCore::isInner(pos); }
		PointL Convex::getOverlappingPoints(const IModel& mdl, const Vec2& inner) const { return ConvexCore::getOverlappingPoints(mdl,inner); }
		CircleCore Convex::getBBCircle() const { return ConvexCore::getBBCircle(); }

		ConvexModel::ConvexModel(): _rflag(RFL_ALL) {}
		ConvexModel::ConvexModel(const PointL& pl): _convex(pl), _rflag(RFL_ALL) {}
		ConvexModel::ConvexModel(PointL&& pl): _convex(pl), _rflag(RFL_ALL) {}

		void ConvexModel::_refreshCInertia() const {
			if(Bit::ChClear(_rflag, RFL_CENTER_INERTIA)) {
				auto res = _convex.area_inertia_center();
				_area = std::get<0>(res);
				_inertia = std::get<1>(res);
				_center = std::get<2>(res);
			}
		}
		float ConvexModel::getArea() const {
			_refreshCInertia();
			return _area;
		}
		float ConvexModel::getInertia() const {
			_refreshCInertia();
			return _inertia;
		}
		const AVec2& ConvexModel::getCenter() const {
			_refreshCInertia();
			return _center;
		}
		const PointL& ConvexModel::getPoint() const { return _convex.point; }
		PointL& ConvexModel::refPoint() {
			_rflag = RFL_ALL;
			return _convex.point;
		}
		void ConvexModel::addOffset(const Vec2& ofs) { _convex.addOffset(ofs); }
		Vec2 ConvexModel::support(const Vec2& dir) const { return _convex.support(dir); }
		Vec2 ConvexModel::center() const { return Vec2(getCenter()); }
		bool ConvexModel::isInner(const Vec2& pos) const { return _convex.isInner(pos); }
		PointL ConvexModel::getOverlappingPoints(const IModel& mdl, const Vec2& inner) const { return _convex.getOverlappingPoints(mdl, inner); }
		CircleCore ConvexModel::getBBCircle() const { return _convex.getBBCircle(); }

		Vec3 Dual(const Plane& plane) {
			float d = plane.d;
			assert(plane.d <= 0);
			if(std::fabs(d) < 1e-3f)
				d = (d<0) ? -1e-3f : 1e-3f;

			return plane.getNormal() / -d;
		}
		Plane Dual(const Vec3& p) {
			float invlen = _sseRcp22Bit(p.length());
			return Plane(p*invlen, -invlen);
		}
		StLineCore Dual(const Vec2& v) {
			Vec2 t0(0,-v.y),
				t1(1, v.x-v.y);
			t1 = (t1-t0).normalization();
			return StLineCore(t0, t1);
		}
		Vec2 Dual(const StLineCore& l) {
			const Vec2& t0 = l.pos;
			Vec2 t1 = l.dir + t0;
			return Vec2((t0.y-t1.y) / (t0.x-t1.x),
						-t0.ccw(t1) / (t1.x-t0.x));
		}
		// ----- 領域の積分計算関数 -----
		namespace {
			//! 凸包を三角形に分割して抗力を計算
			RForce CalcRF_Convex(const RPose& rp0, const RPose& rp1, const ConvexCore& cv, const StLineCore& div) {
				class Tmp {
					public:
						struct TmpIn {
							float height;
							float velN;
							float pos;		// 直線に射影した2D頂点は1次元の数値になる
						};

					private:
						TmpIn	_tmp[2];
						int		_swI = 0;
						const RPose			&_rp0,
											&_rp1;
						const StLineCore	&_div;
						Vec2				_nml;

						void _advance(const Vec2& p) {
							_doSwitch();
							auto& cur = _current();
							cur.height = _div.dir.ccw(p);
							cur.velN = _nml.dot(_rp0.getVelocAt(p) - _rp1.getVelocAt(p));
							cur.pos = _div.posDot(p);
						}
						void _doSwitch() { _swI ^= 1; }
						TmpIn& _current() { return _tmp[_swI]; }
						TmpIn& _prev() { return _tmp[_swI^1]; }

					public:
						Tmp(const RPose& r0, const RPose& r1, const StLineCore& div): _rp0(r0), _rp1(r1), _div(div) {
							// 衝突ライン(2D)の法線
							_nml = div.dir * cs_mRot90[0];
							if(_nml.dot(_rp0.getOffset()) < 0)
								_nml *= -1.f;
						}
						void calcForce(RForce& dst, const ConvexCore& c) {
							const auto& pts = c.point;
							int nV = pts.size();
							if(nV < 3)
								return;

							float p_lin = 0,
								p_tor = 0;
							_advance(pts[0]);
							for(int i=1 ; i<nV ; i++) {
								_advance(pts[i]);
								auto& cur = _current();
								auto& pre = _prev();
								// calc spring
								p_lin += (pre.height + cur.height) * 0.5f;
								p_tor += (1.f/3) * (pre.pos*pre.height + (pre.pos*cur.height + cur.pos*pre.height)*0.5f + cur.pos*cur.height);
								// calc dumper
								p_lin += (pre.velN + cur.velN) * 0.5f;
								p_tor += (3.f/2)*pre.pos*pre.velN +
										 (2.f/3)*pre.pos*cur.velN +
										 cur.pos * cur.velN;
								// TODO: calc dynamic-friction
								// TODO: calc static-friction
							}
							dst.sdump.linear += _nml * p_lin;
							dst.sdump.torque += p_tor;
						}
				};

				Tmp tmp(rp0, rp1, div);
				RForce rf = {};
				// 平面で2つに切ってそれぞれ計算
				auto cvt = cv.splitTwo(div);
				tmp.calcForce(rf, cvt.first);
				tmp.calcForce(rf, cvt.second);
				return rf;
			}
			//! 円同士
			RForce CalcOV_Circle2(const Rigid& r0, const Rigid& r1, const Vec2& inner, const StLineCore& div) {
				// 境界線分を計算して共通領域を2つに分けてそれぞれ積分
				// 中身がCircleだと分かっているのでポインタの読み替え
				throw std::runtime_error("not implemented yet");
			}
			RForce CalcOV_Convex2(const Rigid& r0, const Rigid& r1, const Vec2& inner, const StLineCore& div) {
				// 領域算出
				auto &m0 = reinterpret_cast<const Convex&>(r0.getModel()),
					&m1 = reinterpret_cast<const Convex&>(r1.getModel());
				return CalcRF_Convex(r0.getPose(), r1.getPose(), m0.getOverlappingConvex(m1, inner), div);
			}
			//! 円とBox含む多角形
			RForce CalcOV_CircleConvex(const Rigid& r0, const Rigid& r1, const Vec2& inner, const StLineCore& div) {
				PointL _pts;
				// pts[0]とpts[nV-1]の間は円弧を表す
				// 円弧部分は弓部で分けて凸包部分は三角形で計算し、残りは独自式で積分
				throw std::runtime_error("not implemented yet");
			}
		}

		RForce CalcForce(const Rigid& r0, const Rigid& r1, const Vec2& inner, const StLineCore& div) {
			// 実質Convex, Box, Circle専用
			if(r0.getModel().getCID() == Circle::GetCID()) {
				if(r1.getModel().getCID() == Circle::GetCID())
					return CalcOV_Circle2(r0, r1, inner, div);
				return CalcOV_CircleConvex(r0, r1, inner, div);
			}
			return CalcOV_Convex2(r0, r1, inner, div);
		}
	}
}
