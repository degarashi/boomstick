#include "geom2D.hpp"

namespace boom {
	namespace geo2d {
		Circle::Circle(const Vec2& c, float r): vCenter(c), fRadius(r) {}
		float Circle::bs_getArea() const {
			Assert(Trap, false, "not implemented yet") throw 0;
		}
		float Circle::bs_getInertia() const {
			Assert(Trap, false, "not implemented yet") throw 0;
		}
		const Vec2& Circle::bs_getCenter() const {
			return vCenter;
		}
		const Circle& Circle::bs_getBVolume() const {
			return *this;
		}
		Vec2 Circle::support(const Vec2& dir) const {
			return dir * fRadius + vCenter;
		}
		bool Circle::isInner(const Vec2& pos) const {
			return pos.distance(vCenter) <= fRadius;
		}
		bool Circle::hit(const Vec2& pt) const {
			return vCenter.dist_sq(pt) <= spn::Square(fRadius);
		}
		bool Circle::hit(const Circle& c) const {
			return vCenter.dist_sq(c.vCenter) <= spn::Square(fRadius + c.fRadius);
		}
		Circle Circle::operator * (const AMat32& m) const {
			auto& m2 = reinterpret_cast<const AMat22&>(m);
			Vec2 tx(vCenter + Vec2(fRadius,0)),
				ty(vCenter + Vec2(0,fRadius)),
				ori(vCenter);
			ori = ori.asVec3(1) * m;
			tx = tx * m2 - ori;
			ty = ty * m2 - ori;
			return Circle(ori,
						spn::Sqrt(std::max(tx.len_sq(), ty.len_sq())));
		}
		// 半径のみ倍率をかける
		Circle Circle::operator * (float s) const {
			return Circle(vCenter, fRadius*s);
		}
	}
}
