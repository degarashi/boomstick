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
		bool Circle::hit(const Vec2& pt, float t) const {
			return vCenter.dist_sq(pt) <= spn::Square(fRadius + t);
		}
		bool Circle::hit(const Circle& c, float t) const {
			return vCenter.dist_sq(c.vCenter) <= spn::Square(fRadius + c.fRadius + t);
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
		void Circle::setBoundary(const IModel* p) {
			*this = p->im_getBVolume();
		}
		void Circle::appendBoundary(const IModel* p) {
			Circle c2 = p->im_getBVolume();
			Vec2 toC2 = c2.vCenter - vCenter;
			float len = toC2.normalize(),
				  rd = fRadius - c2.fRadius;
			if(len > rd) {
				// サポート写像で新たな円を算出
				Vec2 tv0 = c2.vCenter + toC2 * c2.fRadius,
					 tv1 = vCenter - toC2 * fRadius;
				vCenter = (tv0 + tv1) * .5f;
				fRadius = tv0.distance(tv1) * .5;
			}
		}
		Circle Circle::Boundary(const IModel* p, size_t n, size_t stride) {
			AssertT(Trap, n>0, (std::invalid_argument)(const char*), "size must not 0")

			Circle c;
			c.setBoundary(p);
			auto pv = reinterpret_cast<uintptr_t>(p);
			pv += stride;
			while(n-- > 1) {
				p = reinterpret_cast<const IModel*>(pv);
				c.appendBoundary(p);
				pv += stride;
			}
			return c;
		}
		// 半径のみ倍率をかける
		Circle Circle::operator * (float s) const {
			return Circle(vCenter, fRadius*s);
		}
		std::ostream& operator << (std::ostream& os, const Circle& c) {
			return os << "Circle(2d) [ center: " << c.vCenter << std::endl
						<< "radius: " << c.fRadius << ']';
		}
	}
}
