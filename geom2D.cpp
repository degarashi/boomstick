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

		// ---------------------- Circle ----------------------
		CircleCore::CircleCore(const Vec2& c, float r): center(c), radius(r) {}
		float CircleCore::area() const {
			throw std::runtime_error("");
		}
		Vec2 CircleCore::support(const Vec2& dir) const {
			return dir * radius + center;
		}

		Vec2 Circle::support(const Vec2& dir) const {
			return CircleCore::support(dir);
		}
		CircleModel::CircleModel(): _rflag(0xff) {}
		CircleModel::CircleModel(const CircleCore& c): _circle(c), _rflag(0xff) {}
		Vec2 Circle::center() const { return CircleCore::center; }
		Vec2 CircleModel::support(const Vec2& dir) const { return _circle.support(dir); }
		Vec2 CircleModel::center() const { return _circle.center; }

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
			return (p - pos).dot(dir);
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
		LineCore::LNear LineCore::nearest(const Vec2& p) const {
			Vec2 toP(p-point[0]),
				toV1(point[1]-point[0]);
			float lenV1 = toV1.length();
			toV1 *= _sseRcp22Bit(lenV1);
			float d = toV1.dot(toP);
			if(d <= 0)
				return LNear(point[0], BEGIN);
			else if(d >= lenV1)
				return LNear(point[1], END);
			else
				return LNear(point[0]+toV1*d, ONLINE);
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

		// ---------------------- Poly ----------------------
		PolyCore::PolyCore(const Vec2& p0, const Vec2& p1, const Vec2& p2): point{p0,p1,p2} {}
		float PolyCore::area() const {
			return CalcArea(point[0], point[1], point[2]);
		}
		Vec2 PolyCore::center() const {
			return (point[0] + point[1] + point[2]) * (1.0f/3);
		}
		float PolyCore::CalcArea(const Vec2& p0, const Vec2& p1, const Vec2& p2) {
			return  (p1-p0).ccw(p2-p0) * 0.5f;
		}
		float PolyCore::CalcArea(const Vec2& p0, const Vec2& p1) {
			return p0.ccw(p1) * 0.5f;
		}
		Vec2 PolyCore::support(const Vec2& dir) const {
			throw std::runtime_error("");
		}
		void PolyCore::addOffset(const Vec2& ofs) {
			for(int i=0 ; i<3 ; i++)
				point[i] += ofs;
		}

		Vec2 Poly::support(const Vec2& dir) const { return PolyCore::support(dir); }
		Vec2 Poly::center() const { return PolyCore::center(); }

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
				if(Vec2::Ccw(tsrc[rc-2], tsrc[rc-1], tsrc[rc]) < 0)
					--pDst;
				*pDst++ = tsrc[rc];
			}
			*pDst++ = tsrc[nV-1];
			*pDst++ = tsrc[nV-2];
			for(int rc=nV-3 ; rc>=0 ; rc--) {
				if(Vec2::Ccw(tsrc[rc+2], tsrc[rc+1], tsrc[rc]) < 0)
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
		Vec2 Convex::support(const Vec2& dir) const {
			return ConvexCore::support(dir);
		}
		Vec2 Convex::center() const {
			return ConvexCore::center();
		}

		ConvexModel::ConvexModel(): _rflag(RFL_ALL) {}
		ConvexModel::ConvexModel(const ConvexCore::PointL& pl): _convex(pl), _rflag(RFL_ALL) {}
		ConvexModel::ConvexModel(ConvexCore::PointL&& pl): _convex(pl), _rflag(RFL_ALL) {}

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
		const ConvexCore::PointL& ConvexModel::getPoint() const { return _convex.point; }
		ConvexCore::PointL& ConvexModel::refPoint() {
			_rflag = RFL_ALL;
			return _convex.point;
		}
		void ConvexModel::addOffset(const Vec2& ofs) {
			_convex.addOffset(ofs);
		}
		Vec2 ConvexModel::support(const Vec2& dir) const {
			return _convex.support(dir);
		}
		Vec2 ConvexModel::center() const {
			return Vec2(getCenter());
		}
	}
}
