#include "geom2D.hpp"

namespace boom {
	namespace geo2d {
		AABB::AABB(const Vec2& min_v, const Vec2& max_v): minV(min_v), maxV(max_v) {}
		Vec2 AABB::support(const Vec2& dir) const {
			constexpr uint32_t one = 0x3f800000;
			Vec2 v(spn::ReinterpretValue<float>((spn::ReinterpretValue<int>(dir.x) & 0x80000000) | one),
					spn::ReinterpretValue<float>((spn::ReinterpretValue<int>(dir.y) & 0x80000000) | one));
			return v *= (maxV - minV) * 0.5f;
		}

		bool AABB::hit(const Segment& l) const {
			Vec2 dir(minV.x, maxV.y - minV.y);
			if(Segment(minV, Vec2(minV.x,maxV.y)).crossPoint(l).second == LinePos::OnLine)
				return true;
			if(Segment(Vec2(maxV.x,minV.y), maxV).crossPoint(l).second == LinePos::OnLine)
				return true;
			return false;
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
	}
}
