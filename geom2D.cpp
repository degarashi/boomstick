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

		CircleCore IModel::getBCircle() const {
			throw std::runtime_error("not supported function"); }
		Convex2 IModel::splitTwo(const StLineCore& line) const {
			throw std::runtime_error("not supported function"); }

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
		Vec2 Point::getCenter() const {
			return *this;
		}
		CircleCore Point::getBCircle() const {
			// 円の半径が0だと点同士の時にヒットしないので微量含める
			return CircleCore(*this, 1e-6f);
		}
		// ---------------------- Box ----------------------
		BoxCore::BoxCore(const Vec2& min_v, const Vec2& max_v): minV(min_v), maxV(max_v) {}

		Vec2 BoxCore::support(const Vec2& dir) const {
			constexpr uint32_t one = 0x3f800000;
			Vec2 v(spn::ReinterpretValue<float>((spn::ReinterpretValue<int>(dir.x) & 0x80000000) | one),
					spn::ReinterpretValue<float>((spn::ReinterpretValue<int>(dir.y) & 0x80000000) | one));
			return v *= (maxV - minV) * 0.5f;
		}

		bool BoxCore::hit(const LineCore& l) const {
			Vec2 dir(minV.x, maxV.y - minV.y);
			if(LineCore(minV, Vec2(minV.x,maxV.y)).crossPoint(l).second == LINEPOS::ONLINE)
				return true;
			if(LineCore(Vec2(maxV.x,minV.y), maxV).crossPoint(l).second == LINEPOS::ONLINE)
				return true;
			return false;
		}
		Vec2 BoxCore::center() const {
			return (minV + maxV) * 0.5f;
		}
		Vec2 BoxCore::nearest(const Vec2& pos) const {
			return pos.getMax(minV).getMin(maxV);
		}
		float BoxCore::inertia() const {
			throw std::runtime_error("not implemented yet");
		}
		float BoxCore::area() const {
			Vec2 sz = maxV - minV;
			return sz.x * sz.y;
		}
		CircleCore BoxCore::bcircle() const {
			// 対角線 = 直径
			return CircleCore((minV + maxV) * 0.5f,
								minV.distance(maxV));
		}

		CircleCore Box::getBCircle() const { return BoxCore::bcircle(); }
		Vec2 Box::support(const Vec2& dir) const { return BoxCore::support(dir); }
		Vec2 Box::getCenter() const { return BoxCore::center(); }

		BoxModel::BoxModel(const BoxCore& b): _box(b) {}
		CircleCore BoxModel::getBCircle() const { return getCache(TagBCircle()); }
		Vec2 BoxModel::support(const Vec2& dir) const { return _box.support(dir); }
		Vec2 BoxModel::getCenter() const { return getCache(TagCenter()); }
		float BoxModel::getArea(bool) const { return getCache(TagArea()); }
		float BoxModel::getInertia(bool) const { return getCache(TagInertia()); }

		// ---------------------- Circle ----------------------
		CircleCore::CircleCore(const Vec2& c, float r): vCenter(c), fRadius(r) {}
		float CircleCore::area() const {
			throw std::domain_error("not implemented yet");
		}
		float CircleCore::inertia() const {
			throw std::domain_error("not implemented yet");
		}
		Vec2 CircleCore::support(const Vec2& dir) const {
			return dir * fRadius + vCenter;
		}
		bool CircleCore::hit(const Vec2& pt) const {
			return vCenter.dist_sq(pt) <= spn::Square(fRadius);
		}
		bool CircleCore::hit(const CircleCore& c) const {
			return vCenter.dist_sq(c.vCenter) <= spn::Square(fRadius + c.fRadius);
		}
		CircleCore CircleCore::operator * (const AMat32& m) const {
			auto& m2 = reinterpret_cast<const spn::AMat22&>(m);
			Vec2 tx(vCenter + Vec2(fRadius,0)),
				ty(vCenter + Vec2(0,fRadius));
			tx = tx * m2 - vCenter;
			ty = ty * m2 - vCenter;
			return CircleCore((vCenter.asVec3(1)*m),
							  spn::_sseSqrt(std::max(tx.len_sq(), ty.len_sq())));
		}

		Vec2 Circle::support(const Vec2& dir) const { return CircleCore::support(dir); }
		Vec2 Circle::getCenter() const { return CircleCore::vCenter; }
		bool Circle::isInner(const Vec2& pos) const { return CircleCore::hit(pos); }
		CircleCore Circle::getBCircle() const { return *this; }
		float Circle::getArea(bool) const { return CircleCore::area(); }
		float Circle::getInertia(bool) const { return CircleCore::inertia(); }

		CircleModel::CircleModel(const CircleCore& c): Cache(c) {}
		Vec2 CircleModel::support(const Vec2& dir) const { return _core.support(dir); }
		Vec2 CircleModel::getCenter() const { return _core.vCenter; }
		bool CircleModel::isInner(const Vec2& pos) const { return _core.hit(pos); }
		CircleCore CircleModel::getBCircle() const { return _core; }
		float CircleModel::getArea(bool) const { return getCache(TagArea()); }
		float CircleModel::getInertia(bool) const { return getCache(TagInertia()); }

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
		StLineCore StLineCore::operator * (const AMat32& m) const {
			return StLineCore{pos.asVec3(1)*m, dir.asVec3(0)*m};
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
		RayCore RayCore::operator * (const AMat32& m) const {
			return RayCore{pos.asVec3(1)*m, dir.asVec3(0)*m};
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
			auto fn = [](float f){ return f; };
			Vec2x2 v2 = NearestPoint(toStLine(), l.toStLine(), fn, fn);
			return LNear(v2.first, online(v2.first) ? LINEPOS::ONLINE : LINEPOS::NOTHIT);
		}
		LNear LineCore::crossPoint(const StLineCore& l) const {
			float c0 = l.dir.ccw(point[0]-l.pos),
				c1 = l.dir.ccw(point[1]-l.pos);
			if(c0*c1 <= 0) {
				Vec2 diff(point[1]-point[0]);
				c0 = std::fabs(c0);
				float d = c0 / (c0+std::fabs(c1));
				return LNear(point[0] + diff*d, LINEPOS::ONLINE);
			}
			return LNear(Vec2(), LINEPOS::NOTHIT);
		}
		LineCore LineCore::operator * (const AMat32& m) const {
			return LineCore{point[0].asVec3(1)*m, point[1].asVec3(1)*m};
		}

		Vec2 Line::support(const Vec2& dir) const {
			return LineCore::support(dir);
		}
		Vec2 Line::getCenter() const {
			return (point[0] + point[1]) * 0.5f;
		}
		CircleCore LineCore::getBCircle() const {
			return CircleCore(point[0] + point[1] * 0.5f,
							point[0].distance(point[1]));
		}
		CircleCore Line::getBCircle() const { return LineCore::getBCircle(); }

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
		float PolyCore::inertia() const {
			return (1.0f/18) * (point[0].dot(point[0])
									+ point[0].dot(point[0])
									+ point[0].dot(point[0])
									- point[1].dot(point[2])
									- point[2].dot(point[0])
									- point[0].dot(point[1]));
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
		CircleCore PolyCore::bcircle() const {
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
		Vec2 Poly::getCenter() const { return PolyCore::center(); }
		bool Poly::isInner(const Vec2& pos) const { return isInTriangle(pos); }
		CircleCore Poly::getBCircle() const { return PolyCore::bcircle(); }

		PolyModel::PolyModel(const Vec2& p0, const Vec2& p1, const Vec2& p2): Cache(p0,p1,p2) {}
		void PolyModel::setPoint(int n, const Vec2& v) {
			_core.point[n] = v;
			Cache::setCacheFlagAll();
		}
		void PolyModel::addOffset(const Vec2& ofs) {
			_core.addOffset(ofs);
		}
		Vec2 PolyModel::support(const Vec2& dir) const {
			return _core.support(dir);
		}
		bool PolyModel::isInner(const Vec2& pos) const {
			return _core.isInTriangle(pos);
		}
		Vec2 PolyModel::getCenter() const { return getCache(TagCenter()); }
		float PolyModel::getArea(bool) const { return getCache(TagArea()); }
		float PolyModel::getInertia(bool) const { return getCache(TagInertia()); }
		CircleCore PolyModel::getBCircle() const { return getCache(TagBCircle()); }

		// ---------------------- Convex ----------------------
		ConvexCore::ConvexCore(const PointL& pl): point(pl) {}
		ConvexCore::ConvexCore(PointL&& pl): point(pl) {}
		ConvexCore::ConvexCore(std::initializer_list<Vec2> v): point(v.size()) {
			auto itrD = point.begin();
			auto itr = v.begin();
			while(itr != v.end())
				*itrD++ = *itr++;
		}
		ConvexCore::ConvexCore(ConvexCore&& c): point(std::forward<PointL>(c.point)) {}
		ConvexCore& ConvexCore::operator = (ConvexCore&& c) {
			point = std::move(c.point);
			return *this;
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
		LineCore ConvexCore::getOuterLine(int n) const {
			return LineCore(point[(n+1)%point.size()], point[n]);
		}
		StLineCore ConvexCore::getOuterStLine(int n) const {
			return StLineCore(point[n], (point[spn::CndSub(n+1, point.size())] - point[n]).normalization());
		}
		std::ostream& ConvexCore::dbgPrint(std::ostream& os) const {
			for(auto& p : point)
				os << '[' << p.x << ',' << p.y << ']' << std::endl;
			return os;
		}

		Convex::Convex(ConvexCore&& c): ConvexCore(std::forward<ConvexCore>(c)) {}
		Convex2 Convex::splitTwo(const StLineCore& l) const {
			return Convex2(std::move(ConvexCore::splitTwo(l)));
		}
		Convex2 ConvexModel::splitTwo(const StLineCore& l) const {
			return Convex2(std::move(_core.splitTwo(l)));
		}
		std::ostream& ConvexModel::dbgPrint(std::ostream& os) const {
			return _core.dbgPrint(os);
		}
		ConvexCore& ConvexCore::operator *= (const AMat32& m) {
			for(auto& p : point)
				p = p.asVec3(1) * m;
			return *this;
		}
		Convex& Convex::operator *= (const AMat32& m) {
			((ConvexCore&)*this) *= m;
			return *this;
		}

		namespace {
			//! c0の1領域に対してc1がめり込んでいる場合の領域抽出
			auto Func(const IModel& c0, const IModel& c1, const Vec2& inner) -> Convex {
				// innerをc0ローカルにして判定
				auto ip = c0.checkPosition(inner);
				const Vec2 &p0 = c0.getPoint(ip.second),
							&p1 = c0.getPoint(spn::CndSub(ip.second+1, c0.getNPoints()));
				// lineをc1ローカルにして分割 -> ワールド座標系に変換
				StLineCore line(p0, (p1-p0).normalization());
				return c1.splitTwo(line).second;
			};
		}
		std::ostream& operator << (std::ostream& os, const IModel& mdl) {
			return mdl.dbgPrint(os);
		}
		std::ostream& Convex::dbgPrint(std::ostream& os) const {
			int nP = getNPoints();
			if(nP > 0) {
				for(int i=0 ; i<nP-1 ; i++)
					os  << getPoint(i) << std::endl;
				os << getPoint(nP-1);
			}
			return os;
		}

		Convex Convex::GetOverlappingConvex(const IModel& m0, const IModel& m1, const Vec2& inner) {
			if(!m0.isInner(inner) || !m1.isInner(inner))
				throw std::runtime_error("invalid inner point");

			// m0がm1にめり込んでいる頂点リストを出力
			auto res0 = m0.getOverlappingPoints(m1,inner);
			if(res0.first)
				return Convex(std::move(res0.second));
			// m1がm0にめり込んでいる頂点リストを出力
			auto res1 = m1.getOverlappingPoints(m0,inner);
			if(res1.first)
				return Convex(std::move(res1.second));

			const PointL &pt0 = res0.second,
						&pt1 = res1.second;
			int nV0 = pt0.size(),
				nV1 = pt1.size();
			if(nV0 == 0) {
				assert(nV1 >= 3);	// これがfalseなら、2つのConvexは重なっていない可能性がある
				// m1のめり込んだ頂点全部がm0の1つの辺に収まっている
				return Func(m0, m1, inner);
			} else if(nV1 == 0) {
				assert(nV0 >= 3);	// これがfalseなら、2つのConvexは重なっていない可能性がある
				// m0のめり込んだ頂点全部がm1の1つの辺に収まっている
				return Func(m1, m0, inner);
			} else {
				assert(nV0>=3 && nV1>=3);
				PointL ptDst(nV0 + nV1 - 2);

				auto* pDst = &ptDst[0];
				// 繋ぎ目の処理: m0リストの始点をRayに変換した物とm1の終端で頂点を1つ
				auto res = LineCore(pt0[0], pt0[1]).crossPoint(LineCore(pt1[nV1-2], pt1[nV1-1]));
				assert(res.second == LINEPOS::ONLINE);
				*pDst++ = res.first;
				for(int i=1 ; i<nV0-1 ; i++)
					*pDst++ = pt0[i];
				// m1リストの始点をRayに変換した物とm0の終端で頂点を1つ
				res = LineCore(pt0[nV0-2], pt0[nV0-1]).crossPoint(LineCore(pt1[0], pt1[1]));
				assert(res.second == LINEPOS::ONLINE);
				*pDst++ = res.first;
				for(int i=1 ; i<nV1-1 ; i++)
					*pDst++ = pt1[i];
				return Convex(std::move(ptDst));
			}
		}
		int Convex::getNPoints() const { return point.size(); }
		Vec2 Convex::getPoint(int n) const { return point[n]; }
		IModel::CPos Convex::checkPosition(const Vec2& pos) const { return ConvexCore::checkPosition(pos); }

		namespace {
			int BinSearch(const PointL& point, int nV, const IModel& mdl, int a, int b, bool flip) {
				for(;;) {
					if(a+1 >= b)
						break;
					int c = (a+b)/2;
					if(mdl.isInner(point[spn::CndRange(c, nV)]) ^ flip)
						b = c;
					else
						a = c;
				}
				return spn::CndRange(a, nV);
			}
		}
		IModel::PosL ConvexCore::getOverlappingPoints(const IModel& mdl, const Vec2& inner) const {
			auto res = checkPosition(inner);
			if(res.first != IModel::POSITION::OUTER) {
				int nV = point.size();
				int a,b, begI;
				// 2分探索で衝突が始まる地点を探す
				if(mdl.isInner(point[res.second])) {
					a = res.second - nV;
					b = res.second;
					// 衝突開始インデックス(これ自体は衝突していない)
					begI = BinSearch(point, nV, mdl, a,b,false);
				} else {
					begI = res.second;
					if(!mdl.isInner(point[spn::CndSub(res.second+1, nV)]))
						return std::make_pair(false, PointL());
				}

				// 衝突が終わる地点を探す
				a = begI + 1;
				b = a + nV;
				int endI = BinSearch(point,nV,mdl,a,b,true);		// 衝突終了インデックス(これ自体衝突している)

				if(begI == endI) {
					// 全部出力
					return std::make_pair(true, point);
				}
				endI = spn::CndSub(endI+1, nV);
				endI = spn::CndAdd(endI, begI, nV)+1;
				PointL pts(endI - begI);
				auto* pDst = &pts[0];
				// 開始地点から終了地点までを出力
				while(begI != endI) {
					*pDst++ = point[spn::CndSub(begI,nV)];
					++begI;
				}
				return std::make_pair(false, std::move(pts));
			}
			return std::make_pair(false, PointL());
		}

		IModel::CPos ConvexCore::checkPosition(const Vec2& pos) const {
			// 適当に内部点を算出
			Vec2 inner = (point[0] + point[1] + point[2]) * (1.f/3);
			Vec2 toP(pos - inner);
			if(toP.len_sq() < 1e-6f) {
				// 重心がちょうどposと重なってしまったら少しずらす
				inner.lerp(point[0], 0.5f);
			}

			// 内部のどの三角形に該当するか2分探索
			int nV = point.size();
			struct Tmp {
				GAP_VECTOR(vec, 2,
						   (int index)
						   (float cw))

				Tmp& operator = (const Tmp& t) {
					_mm_store_ps(vec.m, _mm_load_ps(t.vec.m));
					return *this;
				}
			};
			// tmpA = point[0]へのベクトル
			Tmp tmpA, tmpB, tmpC;
			tmpA.vec = point[0]-inner;
			tmpA.index = 0;
			tmpA.cw = tmpA.vec.cw(toP);
			// tmpB = point[nV]へのベクトル
			tmpB.index = nV;
			tmpB.vec = tmpA.vec;
			tmpB.cw = tmpA.cw;

			while(tmpA.index+1 < tmpB.index) {
				// 捜査対象範囲の半分で区切る => tmpC
				tmpC.index = (tmpA.index+tmpB.index) / 2;
				tmpC.vec = point[tmpC.index] - inner;
				tmpC.cw = tmpC.vec.cw(toP);
				float crAC = tmpA.vec.cw(tmpC.vec);

				if(tmpA.cw >= 0) {
					// posはtmpAの右側
					if(crAC >= 0) {
						// tmpCはtmpAの右側
						if(tmpC.cw <= 0) {
							// posはtmpCの左側
							tmpB = tmpC;
						} else
							tmpA = tmpC;
					} else
						tmpB = tmpC;
				} else {
					// posはtmpAの左側
					if(crAC >= 0) {
						// tmpCはtmpAの右側
						tmpA = tmpC;
					} else {
						if(tmpC.cw >= 0) {
							// posはtmpCの右側
							tmpA = tmpC;
						} else
							tmpB = tmpC;
					}
				}
			}
			tmpB.index = spn::CndSub(tmpB.index, nV);
			float d = (point[tmpB.index] - point[tmpA.index]).cw(pos - point[tmpA.index]);
			auto pA = point[tmpA.index],
					pB = point[tmpB.index];
			IModel::POSITION posf;
			if(d > 1e-6f)
				posf = IModel::POSITION::INNER;
			else if(d < -1e-6f)
				posf = IModel::POSITION::OUTER;
			else
				posf = IModel::POSITION::ONLINE;
			return IModel::CPos(posf, tmpA.index);
		}
		bool ConvexCore::isInner(const Vec2& pos) const {
			return checkPosition(pos).first != IModel::POSITION::OUTER;
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
		ConvexCore2 ConvexCore::splitTwo(const StLineCore& l) const {
			int nV = point.size();
			PointL pt0(nV+1), pt1(nV+1);		// 最初に最大容量確保しておき、後で縮める
			Vec2* pDst[2] = {&pt0[0], &pt1[0]};

			constexpr float DOT_THRESHOLD = 1e-5f;
			// 右側は0, ライン上は前回値, 左側は1を返す
			auto fcheck = [&l](int prev, const Vec2& pt) -> int {
				float c = l.dir.ccw(pt - l.pos);
				if(c >= DOT_THRESHOLD)
					return 0;
				if(c <= -DOT_THRESHOLD)
					return 1;
				return prev;
			};
			// 走査の始点はライン上でない頂点とする
			int cur = 0;
			int flag;
			while(cur != nV) {
				const auto& p = point[cur++];
				flag = fcheck(-1, p);
				if(flag != -1)
					break;
				// とりあえず片方のリストに貯めておく
				*pDst[0]++ = p;
			}
			// 2点以上がライン上なら分割なし
			if(cur > 2) {
				if(flag == 0x01) {
					// すべて右側
					return ConvexCore2({}, point);
				} else {
					// すべて左側
					return ConvexCore2(point, {});
				}
			}
			if(flag == 0x01) {
				std::swap(pt0, pt1);
				std::swap(pDst[0], pDst[1]);
			}
			int prev = flag;

			auto fadd = [&pDst, &l](const Vec2& pPrev, const Vec2& pCur, int flg) {
				switch(flg) {
					case 0x03:		// Left -> Left
						*pDst[1]++ = pPrev;
						break;
					case 0x02: {	// Left -> Right
						auto res = LineCore(pPrev, pCur).crossPoint(l);
						assert(res.second == LINEPOS::ONLINE);
						*pDst[1]++ = pPrev;
						if(pPrev.dist_sq(res.first) >= 1e-6f)
							*pDst[1]++ = res.first;
						*pDst[0]++ = res.first;
						break; }
					case 0x01: {	// Right -> Left
						auto res = LineCore(pPrev, pCur).crossPoint(l);
						assert(res.second == LINEPOS::ONLINE);
						*pDst[0]++ = pPrev;
						if(pPrev.dist_sq(res.first) >= 1e-6f)
							*pDst[0]++ = res.first;
						*pDst[1]++ = res.first;
						break; }
					case 0x00:		// Right -> Right
						*pDst[0]++ = pPrev;
						break;
				}
			};
			for(int i=cur ; i<nV ; i++) {
				// 正数が左側、負数は右側
				const auto& p = point[i];
				prev = fcheck(prev, p);
				flag = ((flag<<1) | prev) & 0x03;
				fadd(point[i-1], p, flag);
			}
			prev = fcheck(prev, point[0]);
			flag = ((flag<<1) | prev) & 0x03;
			fadd(point[nV-1], point[0], flag);

			if(pDst[0] - &pt0[0] < 2)
				pDst[0] = &pt0[0];
			else if(pDst[1] - &pt1[0] < 2)
				pDst[1] = &pt1[0];
			else {
				if(cur == 2) {
					// 最初の1点がライン上だったので頂点を複製
					if((flag & 1) == 1)
						*pDst[1]++ = pt0[0];
					else
						*pDst[0]++ = pt1[0];
				}
			}

			pt0.resize(pDst[0] - &pt0[0]);
			pt1.resize(pDst[1] - &pt1[0]);
			return ConvexCore2(std::move(pt0), std::move(pt1));
		}
		ConvexCore ConvexCore::split(const StLineCore& l) {
			auto res = splitTwo(l);
			std::swap(res.first, *this);
			return std::move(res.second);
		}
		void ConvexCore::splitThis(const StLineCore& l) {
			split(l);
		}
		CircleCore ConvexCore::bcircle() const {
			// 多分遅いアルゴリズムだが、今はこれで我慢
			// 全ての3点の組み合わせを調べる
			int nV = point.size();
			assert(nV >= 3);
			CircleCore c;
			c.fRadius = -1;
			for(int i=0 ; i<nV-2 ; i++) {
				for(int j=i+1 ; j<nV-1 ; j++) {
					for(int k=j+1 ; k<nV ; k++) {
						PolyCore p(point[i], point[j], point[k]);
						auto tc = p.bcircle();
						if(c.fRadius < tc.fRadius)
							c = tc;
					}
				}
			}
			return c;
		}
		std::tuple<bool,Vec2,Vec2> ConvexCore::checkCrossingLine(const StLineCore& l) const {
			//TODO: 効率の良い算出方法を探す
			Vec2 pt[2];
			Vec2* ppt = pt;
			int nV = point.size();
			for(int i=0 ; i<nV ; i++) {
				auto res = LineCore(point[i], point[spn::CndSub(i+1, nV)]).crossPoint(l);
				if(res.second != LINEPOS::NOTHIT) {
					*ppt++ = res.first;
					if(ppt == pt+2)
						break;
				}
			}
			return std::make_tuple(ppt == pt+2, pt[0], pt[1]);
		}

		Vec2 Convex::support(const Vec2& dir) const { return ConvexCore::support(dir); }
		Vec2 Convex::getCenter() const { return std::get<2>(ConvexCore::area_inertia_center()); }
		bool Convex::isInner(const Vec2& pos) const { return ConvexCore::isInner(pos); }
		IModel::PosL Convex::getOverlappingPoints(const IModel& mdl, const Vec2& inner) const { return ConvexCore::getOverlappingPoints(mdl,inner); }
		CircleCore Convex::getBCircle() const { return ConvexCore::bcircle(); }

		ConvexModel::ConvexModel(std::initializer_list<Vec2> v): Cache(v) {}
		ConvexModel::ConvexModel(const PointL& pl): Cache(pl) {}
		ConvexModel::ConvexModel(PointL&& pl): Cache(pl) {}
		float ConvexModel::getArea(bool) const { return getCache(TagArea()); }
		float ConvexModel::getInertia(bool) const { return getCache(TagInertia()); }

		const PointL& ConvexModel::getPoint() const { return _core.point; }
		PointL& ConvexModel::refPoint() {
			setCacheFlagAll();
			return _core.point;
		}
		Vec2 ConvexModel::support(const Vec2& dir) const { return _core.support(dir); }
		Vec2 ConvexModel::getCenter() const { return Cache::getCache(TagCenter()); }
		bool ConvexModel::isInner(const Vec2& pos) const { return _core.isInner(pos); }
		IModel::PosL ConvexModel::getOverlappingPoints(const IModel& mdl, const Vec2& inner) const { return _core.getOverlappingPoints(mdl, inner); }
		CircleCore ConvexModel::getBCircle() const { return Cache::getCache(TagBCircle()); }
		void ConvexModel::addOffset(const Vec2& ofs) {
			_core.addOffset(ofs);
			// 中心座標をずらす
			refCacheNF(TagCenter()) += ofs;
			refCacheNF(TagBCircle()).vCenter += ofs;
		}
		int ConvexModel::getNPoints() const { return _core.point.size(); }
		Vec2 ConvexModel::getPoint(int n) const { return _core.point[n]; }
		IModel::CPos ConvexModel::checkPosition(const Vec2& pos) const { return _core.checkPosition(pos); }

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
		namespace {
			constexpr float MIN_DUAL2DIST = 1e-5f;
			float Clip(float val) {
				if(std::fabs(val) < MIN_DUAL2DIST)
					return val<0 ? -MIN_DUAL2DIST : MIN_DUAL2DIST;
				return val;
			}
		}
		StLineCore Dual2(const Vec2& v) {
			float rlen = _sseRcp22Bit(Clip(v.length()));
			Vec2 pos = v * rlen,
				dir = v * cs_mRot90[0] * rlen;
			return StLineCore(pos, dir);
		}
		Vec2 Dual2(const StLineCore& l) {
			float rlen = _sseRcp22Bit(Clip(l.pos.length()));
			return l.dir * cs_mRot90[1] * rlen;
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
	}
}
