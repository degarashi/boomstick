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
		bool Circle::hit(const Segment& s, float t) const {
			Vec2 np = s.nearest(vCenter).first;
			return np.dist_sq(vCenter) <= spn::Square(fRadius + t);
		}
		Circle Circle::operator * (const AMat32& m) const {
			auto& m2 = reinterpret_cast<const AMat22&>(m);
			Vec2 tx(fRadius,0),
				ty(0,fRadius),
				ori(vCenter);
			ori = ori.asVec3(1) * m;
			tx = tx * m2;
			ty = ty * m2;
			return Circle(ori,
						spn::Sqrt(std::max(tx.len_sq(), ty.len_sq())));
		}
		Circle& Circle::operator += (const Vec2& ofs) {
			vCenter += ofs;
			return *this;
		}
		void Circle::setBoundary(const IModel* p) {
			*this = p->im_getBVolume();
		}
		void Circle::appendBoundary(const IModel* p) {
			Circle c2 = p->im_getBVolume();
			Vec2 toC2 = c2.vCenter - vCenter;
			float lensq = toC2.len_sq();
			if(lensq >= ZEROVEC_LENGTH_SQ) {
				toC2 *= spn::RSqrt(lensq);
				Vec2 tv(support(toC2) - c2.support(-toC2));
				float r_min = std::min(fRadius, c2.fRadius);
				if(tv.dot(toC2) < 0 ||
					tv.len_sq() < spn::Square(r_min*2))
				{
					// 新たな円を算出
					Vec2 tv0(vCenter - toC2*fRadius),
						 tv1(c2.vCenter + toC2*c2.fRadius);
					fRadius = tv0.distance(tv1) * .5f;
					vCenter = (tv0+tv1) * .5f;
				} else {
					// 円が内包されている
					if(fRadius < c2.fRadius) {
						fRadius = c2.fRadius;
						vCenter = c2.vCenter;
					}
				}
			} else {
				// 円の中心が同じ位置にある
				fRadius = std::max(c2.fRadius, fRadius);
			}
		}
		// 半径のみ倍率をかける
		Circle Circle::operator * (float s) const {
			return Circle(vCenter, fRadius*s);
		}
		void Circle::distend(float width, float mindist) {
			fRadius = std::max(mindist, fRadius+width);
		}
		std::ostream& operator << (std::ostream& os, const Circle& c) {
			return os << "Circle(2d) [ center: " << c.vCenter << std::endl
						<< "radius: " << c.fRadius << ']';
		}
	}
}
