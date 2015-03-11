#include "geom2D.hpp"

namespace boom {
	namespace geo2d {
		AABB::AABB(const Vec2& min_v, const Vec2& max_v): minV(min_v), maxV(max_v) {}
		Vec2 AABB::support(const Vec2& dir) const {
			return {(dir.x > 0 ? maxV.x : minV.x),
					(dir.y > 0 ? maxV.y : minV.y)};
		}

		bool AABB::hit(const Vec2& pos, float t) const {
			return spn::IsInRange(pos.x, minV.x-t, maxV.x+t) &&
					spn::IsInRange(pos.y, minV.y-t, maxV.y+t);
		}
		bool AABB::hit(const AABB& ab, float t) const {
			if(ab.maxV.x < minV.x - t)
				return false;
			if(ab.minV.x > maxV.x + t)
				return false;
			if(ab.maxV.y < minV.y - t)
				return false;
			if(ab.minV.y > maxV.y + t)
				return false;
			return true;
		}
		int AABB::_getAreaNumX(float p) const {
			if(p < minV.x)
				return 0x00;
			if(p < maxV.x)
				return 0x01;
			return 0x02;
		}
		int AABB::_getAreaNumY(float p) const {
			if(p < minV.y)
				return 0x00;
			if(p < maxV.y)
				return 0x01;
			return 0x02;
		}
		void AABB::_makeSegmentX(Segment& s, int num) const {
			s.from.x = s.to.x = ((num==0x01) ? minV.x : maxV.x);
			s.from.y = minV.y;
			s.to.y = maxV.y;
		}
		void AABB::_makeSegmentY(Segment& s, int num) const {
			s.from.y = s.to.y = ((num==0x01) ? minV.y : maxV.y);
			s.from.x = minV.x;
			s.to.x = maxV.x;
		}
		bool AABB::_checkHitX(const Segment& s0, int from, int to) const {
			if(to < from)
				std::swap(from, to);
			Segment s1;
			while(from != to) {
				_makeSegmentX(s1, ++from);
				if(s0.hit(s1))
					return true;
			}
			return false;
		}
		bool AABB::_checkHitY(const Segment& s0, int from, int to) const {
			if(to < from)
				std::swap(from, to);
			Segment s1;
			while(from != to) {
				_makeSegmentY(s1, ++from);
				if(s0.hit(s1))
					return true;
			}
			return false;
		}
		std::pair<int,int> AABB::_getAreaNum(const Vec2& p) const {
			return std::make_pair(_getAreaNumX(p.x),
									_getAreaNumY(p.y));
		}
		bool AABB::hit(const Segment& l, float t) const {
			auto f0 = _getAreaNum(l.from),
				f1 = _getAreaNum(l.to);
			auto diffX = std::abs(f0.first - f1.first);
			auto diffY = std::abs(f0.second - f1.second);
			// 片方の点が中央位置にあればヒットしている
			if(f0.first==1 && f0.second==1)
				return true;
			if(f1.first==1 && f1.second==1)
				return true;
			// 片方の軸に変化が無ければ領域のみで判定可
			if(diffX==0)
				return f0.first==0x01 && diffY==2;
			if(diffY==0)
				return f0.second==0x01 && diffX==2;

			// X軸について判定
			if(_checkHitX(l, f0.first, f1.first))
				return true;
			// Y軸について判定
			return _checkHitY(l, f0.second, f1.second);
		}
		Vec2 AABB::nearest(const Vec2& pos) const {
			return pos.getMax(minV).getMin(maxV);
		}
		Vec2 AABB::bs_getCenter() const {
			return (minV + maxV) * 0.5f;
		}
		float AABB::bs_getInertia() const {
			AssertT(Trap, false, (std::domain_error)(const char*), "not implemented yet") throw 0;
		}
		float AABB::bs_getArea() const {
			Vec2 sz = maxV - minV;
			return sz.x * sz.y;
		}
		Circle AABB::bs_getBVolume() const {
			// 対角線 = 直径
			return Circle((minV + maxV) * 0.5f,
								minV.distance(maxV));
		}

		void AABB::setBoundary(const IModel* p) {
			minV.x = p->im_support(Vec2(-1,0)).x;
			maxV.x = p->im_support(Vec2(1,0)).x;
			minV.y = p->im_support(Vec2(0,-1)).y;
			maxV.y = p->im_support(Vec2(0,1)).y;
		}
		void AABB::appendBoundary(const IModel* p) {
			minV.x = std::min(minV.x, p->im_support(Vec2(-1,0)).x);
			maxV.x = std::max(maxV.x, p->im_support(Vec2(1,0)).x);
			minV.y = std::min(minV.y, p->im_support(Vec2(0,-1)).y);
			maxV.y = std::max(maxV.y, p->im_support(Vec2(0,1)).y);
		}
		AABB AABB::Boundary(const IModel* p, size_t n, size_t stride) {
			AssertT(Trap, n>0, (std::invalid_argument)(const char*), "size must not 0")

			AABB a;
			a.setBoundary(p);
			auto pv = reinterpret_cast<uintptr_t>(p);
			pv += stride;
			while(n-- > 1) {
				p = reinterpret_cast<const IModel*>(pv);
				a.appendBoundary(p);
				pv += stride;
			}
			return a;
		}
		std::ostream& operator << (std::ostream& os, const AABB& a) {
			return os << "AABB(2d) [ min: " << a.minV << std::endl
						<< "max: " << a.maxV << ']';
		}
	}
}
