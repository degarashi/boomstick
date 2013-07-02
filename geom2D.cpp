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
		RForce::F& RForce::F::operator *= (float s) {
			linear *= s;
			torque *= s;
			return *this;
		}
		RForce::F RForce::F::operator * (float s) const {
			F ret(*this);
			return ret *= s;
		}

		RForce& RForce::operator +=(const RForce& rf) {
			sdump += rf.sdump;
			fricD += rf.fricD;
			return *this;
		}
		RForce& RForce::operator *= (float s) {
			sdump *= s;
			fricD *= s;
			return *this;
		}
		RForce RForce::operator * (float s) const {
			RForce ret(*this);
			return ret *= s;
		}

		std::ostream& operator << (std::ostream& os, const RForce::F& f) {
			return os << "linear: " << f.linear << std::endl
					  << "torque: " << f.torque;
		}
		std::ostream& operator << (std::ostream& os, const RForce& f) {
			return os << "sdump: " << f.sdump << std::endl
					  << "fricD: " << f.fricD;
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
		Vec2 Point::center() const {
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
		CircleCore BoxCore::getBCircle() const {
			// 対角線 = 直径
			return CircleCore((minV + maxV) * 0.5f,
								minV.distance(maxV));
		}

		CircleCore Box::getBCircle() const { return BoxCore::getBCircle(); }
		Vec2 Box::support(const Vec2& dir) const { return BoxCore::support(dir); }
		Vec2 Box::center() const { return BoxCore::center(); }

		BoxModel::BoxModel(const BoxCore& b): _box(b) {}
		CircleCore BoxModel::getBCircle() const { return _box.getBCircle(); }
		Vec2 BoxModel::support(const Vec2& dir) const { return _box.support(dir); }
		Vec2 BoxModel::center() const { return _box.center(); }

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
		bool CircleCore::hit(const CircleCore& c) const {
			return center.dist_sq(c.center) <= spn::Square(radius + c.radius);
		}
		CircleCore CircleCore::operator * (const AMat32& m) const {
			auto& m2 = reinterpret_cast<const spn::AMat22&>(m);
			Vec2 tx(center + Vec2(radius,0)),
				ty(center + Vec2(0,radius));
			tx = tx * m2 - center;
			ty = ty * m2 - center;
			return CircleCore((center.asVec3(1)*m),
							  spn::_sseSqrt(std::max(tx.len_sq(), ty.len_sq())));
		}

		Vec2 Circle::support(const Vec2& dir) const { return CircleCore::support(dir); }
		Vec2 Circle::center() const { return CircleCore::center; }
		bool Circle::isInner(const Vec2& pos) const { return CircleCore::hit(pos); }
		CircleCore Circle::getBCircle() const { return *this; }

		CircleModel::CircleModel(): _rflag(0xff) {}
		CircleModel::CircleModel(const CircleCore& c): _circle(c), _rflag(0xff) {}

		Vec2 CircleModel::support(const Vec2& dir) const { return _circle.support(dir); }
		Vec2 CircleModel::center() const { return _circle.center; }
		bool CircleModel::isInner(const Vec2& pos) const { return _circle.hit(pos); }
		CircleCore CircleModel::getBCircle() const { return _circle; }

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
				Vec2 diff(point[1]-point[0]);
				c0 = std::fabs(c0);
				float d = c0 / (c0+std::fabs(c1));
				return LNear(point[0] + diff*d, LINEPOS::ONLINE);
			}
			return LNear(Vec2(), LINEPOS::NOTHIT);
		}

		Vec2 Line::support(const Vec2& dir) const {
			return LineCore::support(dir);
		}
		Vec2 Line::center() const {
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
		CircleCore PolyCore::getBCircle() const {
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
		CircleCore Poly::getBCircle() const { return PolyCore::getBCircle(); }

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
		CircleCore PolyModel::getBCircle() const {
			if(Bit::ChClear(_rflag, RFL_BBCIRCLE))
				_bbCircle = _poly.getBCircle();
			return _bbCircle;
		}

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
		Vec2 ConvexCore::center() const {
			return std::get<2>(area_inertia_center());
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
			return Convex2(std::move(_convex.splitTwo(l)));
		}
		std::ostream& ConvexModel::dbgPrint(std::ostream& os) const {
			return _convex.dbgPrint(os);
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
			auto Func2(const IModel& c0, const IModel& c1, const Vec2& inner) -> std::pair<bool,PointL> {
				PointL pt0(c0.getOverlappingPoints(c1, inner));
				int nV0 = pt0.size();
				if(nV0 == c0.getNPoints()) {
					if(nV0 == 3) {
						for(int i=0 ; i<3 ; i++) {
							if(!c1.isInner(pt0[i]))
								return std::make_pair(false, std::move(pt0));
						}
					}
					// m1が全部m0に収まっている
					return std::make_pair(true, std::move(pt0));
				}
				return std::make_pair(false, std::move(pt0));
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
			// m0がm1にめり込んでいる頂点リストを出力
			auto res0 = Func2(m0, m1, inner);
			if(res0.first)
				return Convex(std::move(res0.second));
			// m1がm0にめり込んでいる頂点リストを出力
			auto res1 = Func2(m1, m0, inner);
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

		PointL ConvexCore::getOverlappingPoints(const IModel& mdl, const Vec2& inner) const {
			auto res = checkPosition(inner);
			if(res.first != IModel::POSITION::OUTER) {
				int nV = point.size();
				auto binSearch = [this, nV, &mdl](int a, int b, bool flip) -> int {
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
				};

				int a,b, begI;
				// 2分探索で衝突が始まる地点を探す
				if(mdl.isInner(point[res.second])) {
					a = res.second - nV;
					b = a + nV;
					// 衝突開始インデックス(これ自体は衝突していない)
					begI = binSearch(a, b, false);
				} else {
					begI = res.second;
					if(!mdl.isInner(point[spn::CndSub(res.second+1, nV)]))
						return PointL();
				}

				// 衝突が終わる地点を探す
				a = begI + 1;
				b = a + nV;
				int endI = binSearch(a,b, true);		// 衝突終了インデックス(これ自体衝突している)

				if(begI == endI) {
					// 全部出力
					return point;
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
				return std::move(pts);
			}
			return PointL();
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
			Tmp tmpA, tmpB, tmpC;
			tmpA.vec = point[0]-inner;
			tmpA.index = 0;
			tmpA.cw = tmpA.vec.cw(toP);
			tmpB.index = nV;
			tmpB.vec = tmpA.vec;
			tmpB.cw = tmpA.cw;

			while(tmpA.index+1 < tmpB.index) {
				tmpC.index = (tmpA.index+tmpB.index) / 2;
				tmpC.vec = point[tmpC.index] - inner;
				tmpC.cw = tmpC.vec.cw(toP);
				float crAC = tmpA.vec.cw(tmpC.vec);

				if(tmpA.cw >= 0) {
					if(crAC >= 0) {
						if(tmpC.cw <= 0)
							tmpB = tmpC;
						else
							tmpA = tmpC;
					} else
						tmpA = tmpC;
				} else {
					if(crAC >= 0)
						tmpA = tmpC;
					else {
						if(tmpC.cw >= 0)
							tmpA = tmpC;
						else
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
						*pDst0++ = pPrev;
						break;
					case 0x02: {	// Left -> Right
						auto res = LineCore(pPrev, pCur).crossPoint(l);
						assert(res.second == LINEPOS::ONLINE);
						*pDst0++ = pPrev;
						*pDst0++ = res.first;
						*pDst1++ = res.first;
						break; }
					case 0x01: {	// Right -> Left
						auto res = LineCore(pPrev, pCur).crossPoint(l);
						assert(res.second == LINEPOS::ONLINE);
						*pDst1++ = pPrev;
						*pDst1++ = res.first;
						*pDst0++ = res.first;
						break; }
					case 0x00:		// Right -> Right
						*pDst1++ = pPrev;
						break;
				}
			};

			int fchk = fcheck(0, l.dir.ccw(point[0] - l.pos)),
				prev = fchk,
				flag = prev;
			for(int i=1 ; i<nV ; i++) {
				// 正数が左側、負数は右側
				const auto& p = point[i];
				prev = fcheck(prev, l.dir.ccw(p-l.pos));
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
		CircleCore ConvexCore::getBCircle() const {
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
						auto tc = p.getBCircle();
						if(c.radius < tc.radius)
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
		Vec2 Convex::center() const { return ConvexCore::center(); }
		bool Convex::isInner(const Vec2& pos) const { return ConvexCore::isInner(pos); }
		PointL Convex::getOverlappingPoints(const IModel& mdl, const Vec2& inner) const { return ConvexCore::getOverlappingPoints(mdl,inner); }
		CircleCore Convex::getBCircle() const { return ConvexCore::getBCircle(); }

		ConvexModel::ConvexModel(): _rflag(RFL_ALL) {}
		ConvexModel::ConvexModel(std::initializer_list<Vec2> v): _convex(v), _rflag(RFL_ALL) {}
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
		const ConvexCore& ConvexModel::getConvex() const { return _convex; }
		Vec2 ConvexModel::support(const Vec2& dir) const { return _convex.support(dir); }
		Vec2 ConvexModel::center() const { return Vec2(getCenter()); }
		bool ConvexModel::isInner(const Vec2& pos) const { return _convex.isInner(pos); }
		PointL ConvexModel::getOverlappingPoints(const IModel& mdl, const Vec2& inner) const { return _convex.getOverlappingPoints(mdl, inner); }
		CircleCore ConvexModel::getBCircle() const {
			if(Bit::ChClear(_rflag, RFL_BBCIRCLE))
				_bbCircle = _convex.getBCircle();
			return _bbCircle;
		}
		void ConvexModel::addOffset(const Vec2& ofs) {
			_convex.addOffset(ofs);
			// 中心座標をずらす
			_center += ofs;
			_bbCircle.center += ofs;
		}
		int ConvexModel::getNPoints() const { return _convex.point.size(); }
		Vec2 ConvexModel::getPoint(int n) const { return _convex.point[n]; }
		IModel::CPos ConvexModel::checkPosition(const Vec2& pos) const { return _convex.checkPosition(pos); }

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
			class Tmp {
				public:
					struct TmpIn {
						float height;		// 深度
						float velN, velT;	// 衝突断面の相対速度の垂直、水平成分
						float pos;			// 直線に射影した2D頂点は1次元の数値になる
						float fricD;
					};

				private:
					TmpIn	_tmp[2];
					int		_swI = 0;
					const RPose			&_rp0,
										&_rp1;
					const StLineCore	&_div;
					const RCoeff		&_coeff;
					Vec2				_nml;

					void _advance(const Vec2& p) {
						_doSwitch();
						auto& cur = _current();
						cur.height = std::fabs(_div.dir.ccw(p-_div.pos));

						// 物体Bから見た物体Aの相対速度
						Vec2 vel = _rp1.getVelocAt(p) - _rp0.getVelocAt(p);
						cur.velN = _nml.dot(vel);
						cur.velT = _div.dir.dot(vel);
						// 物体Aの重心からの相対座標 (直線方向に対して)
						cur.pos = _div.dir.dot(p);
						cur.fricD = 0;
					}
					void _doSwitch() { _swI ^= 1; }
					TmpIn& _current() { return _tmp[_swI]; }
					TmpIn& _prev() { return _tmp[_swI^1]; }

				public:
					Tmp(const RPose& r0, const RPose& r1, const StLineCore& div, const RCoeff& coeff): _swI(0), _rp0(r0), _rp1(r1), _div(div), _coeff(coeff) {
						// 衝突ライン(2D)の法線
						_nml = div.dir * cs_mRot90[0];
						if(_nml.dot(_rp0.getOffset() - div.pos) > 0)
							_nml *= -1.f;
					}
					void calcForce(RForce& dst, const ConvexCore& c, float sign0) {
						const auto& pts = c.point;
						int nV = pts.size();
						if(nV < 3)
							return;

						float p_lin = 0,
							p_tor = 0,
							p_fdLin = 0,
							p_fdTor = 0;
						_advance(pts[0]);
						for(int i=1 ; i<=nV ; i++) {
							int idx = spn::CndSub(i,nV);
							_advance(pts[idx]);
							auto& cur = _current();
							auto& pre = _prev();
							float area = (cur.pos - pre.pos) * sign0;		// マイナスの場合も有り得る
							// calc spring
							cur.fricD = (pre.height + cur.height) * 0.5f * area * _coeff.spring;
							p_tor += (1.f/3) * (pre.pos*pre.height + (pre.pos*cur.height + cur.pos*pre.height)*0.5f + cur.pos*cur.height) * area * _coeff.spring;
							// calc dumper
							cur.fricD += (pre.velN + cur.velN) * 0.5f * area * _coeff.dumper;
							p_lin += cur.fricD;
							p_tor += (1.f/3) * (pre.pos*pre.velN + (pre.pos*cur.velN + cur.pos*pre.velN)*0.5f + cur.pos*cur.velN) * area * _coeff.dumper;
							// dynamic-friction
							cur.fricD = (pre.velT + cur.velT) * 0.5f * cur.fricD;
							cur.fricD *= _coeff.fricD;
							p_fdLin += cur.fricD;
							p_fdTor += (1.f/3) * (pre.pos*pre.fricD + (pre.pos*cur.fricD + cur.pos*pre.fricD)*0.5f + cur.pos*cur.fricD) * area * _coeff.fricD;
							// TODO: calc static-friction
						}
						dst.sdump.linear += _nml * p_lin;
						dst.sdump.torque += p_tor;
						dst.fricD.linear += _div.dir * p_fdLin;
						dst.fricD.torque += p_fdTor;
					}
			};
			//! 凸包を三角形に分割して抗力を計算
			RForce CalcRF_Convex(const RPose& rp0, const RPose& rp1, const RCoeff& coeff, const ConvexCore& cv, const StLineCore& div) {
				// 直線方向に向かって左側が表で、右が裏
				// st.dot(line)がプラスなら足し、逆なら引く
				Tmp tmp(rp0, rp1, div, coeff);
				RForce rf = {};
				// 直線(断面)で2つに切ってそれぞれ計算
				auto cvt = cv.splitTwo(div);
				tmp.calcForce(rf, cvt.first, 1.f);
				tmp.calcForce(rf, cvt.second, -1.f);
				return rf;
			}
			//! 円同士
			RForce CalcOV_Circle2(const Rigid& r0, const Rigid& r1, const Vec2& inner, const RCoeff& coeff, const StLineCore& div) {
				// 境界線分を計算して共通領域を2つに分けてそれぞれ積分
				// 中身がCircleだと分かっているのでポインタの読み替え
				throw std::runtime_error("not implemented yet");
			}
			RForce CalcOV_Convex2(const Rigid& r0, const Rigid& r1, const Vec2& inner, const RCoeff& coeff, const StLineCore& div) {
				// 領域算出
				Convex cnv = Convex::GetOverlappingConvex(r0, r1, inner);
				std::cout << "OverlappingConvex:" << std::endl << cnv << std::endl;
				return CalcRF_Convex(r0.getPose(), r1.getPose(), coeff, cnv, div);
			}
			//! 円とBox含む多角形
			RForce CalcOV_CircleConvex(const Rigid& r0, const Rigid& r1, const Vec2& inner, const RCoeff& coeff, const StLineCore& div) {
				PointL _pts;
				// pts[0]とpts[nV-1]の間は円弧を表す
				// 円弧部分は弓部で分けて凸包部分は三角形で計算し、残りは独自式で積分
				throw std::runtime_error("not implemented yet");
			}
		}

		RForce CalcForce(const Rigid& r0, const Rigid& r1, const Vec2& inner, const RCoeff& coeff, const StLineCore& div) {
			// 実質Convex, Box, Circle専用
			if(r0.getCID() == Circle::GetCID()) {
				if(r1.getCID() == Circle::GetCID())
					return CalcOV_Circle2(r0, r1, inner, coeff, div);
				return CalcOV_CircleConvex(r0, r1, inner, coeff, div);
			}
			return CalcOV_Convex2(r0, r1, inner, coeff, div);
		}
	}
}
