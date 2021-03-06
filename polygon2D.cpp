#include "geom2D.hpp"

namespace boom {
	namespace geo2d {
		Poly::Poly(const Vec2& p0, const Vec2& p1, const Vec2& p2): point{p0,p1,p2} {}
		float Poly::bs_getArea() const {
			return CalcArea(point[0], point[1], point[2]);
		}
		Vec2 Poly::bs_getCenter() const {
			return (point[0] + point[1] + point[2]) * (1.0f/3);
		}
		float Poly::bs_getInertia() const {
			return (1.0f/18) * (point[0].dot(point[0])
									+ point[0].dot(point[0])
									+ point[0].dot(point[0])
									- point[1].dot(point[2])
									- point[2].dot(point[0])
									- point[0].dot(point[1]));
		}
		Circle Poly::bs_getBCircle() const {
			int id = getObtuseCorner();
			if(id >= 0) {
				// 鈍角を持っていれば直径を使う
				const Vec2 &v0 = point[spn::CndAdd(id-1, 3)],
							&v1 = point[spn::CndSub(id+1, 3)];
				return Circle((v0+v1)*0.5f, v0.distance(v1));
			} else {
				// なければ3点の外接円
				Line line0((point[1]+point[0]) * 0.5f,
							(point[1]-point[0]) * cs_mRot90[0]),
					line1((point[2]+point[0]) * 0.5f,
							(point[2]-point[0]) * cs_mRot90[0]);
				auto vp = line0.nearest(line1);
				return Circle(vp.first,
							vp.first.distance(point[0]));
			}
		}
		AABB Poly::bs_getBBox() const {
			AABB ab(point[0], point[0]);
			for(int i=1 ; i<3 ; i++) {
				ab.minV.selectMin(point[i]);
				ab.maxV.selectMax(point[i]);
			}
			return ab;
		}
		float Poly::CalcArea(const Vec2& p0, const Vec2& p1, const Vec2& p2) {
			return  (p1-p0).cw(p2-p0) * 0.5f;
		}
		float Poly::CalcArea(const Vec2& p0, const Vec2& p1) {
			return p0.cw(p1) * 0.5f;
		}
		Vec2 Poly::support(const Vec2& dir) const {
			int index = -1;
			float d = std::numeric_limits<float>::lowest();
			for(int i=0 ; i<3 ; i++) {
				float td = point[i].dot(dir);
				if(td > d) {
					d = td;
					index = i;
				}
			}
			return point[index];
		}
		void Poly::addOffset(const Vec2& ofs) {
			for(int i=0 ; i<3 ; i++)
				point[i] += ofs;
		}
		int Poly::getObtuseCorner() const {
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
		bool Poly::_isInTriangle(const Vec2& p, float threshold) const {
			return (point[1] - point[0]).cw(p - point[0]) >= -threshold &&
					(point[2] - point[1]).cw(p - point[1]) >= -threshold &&
					(point[0] - point[2]).cw(p - point[2]) >= -threshold;
		}
		bool Poly::isInTriangle(const Vec2& p) const {
			return _isInTriangle(p, 0.f);
		}
		bool Poly::hit(const Vec2& p, float t) const {
			return _isInTriangle(p, t);
		}
		bool Poly::hit(const Poly& p, float t) const {
			for(int i=0 ; i<3 ; i++) {
				const auto& p0 = point[i];
				for(int j=0 ; j<3 ; j++) {
					const auto& p1 = p.point[j];
					Line line(p0, (p1-p0).normalization());
					auto res0 = checkSide(line, t),
						res1 = p.checkSide(line, t);
					if(res0 != res1 &&
						res0 != LineDivision::Bridge &&
						res1 != LineDivision::Bridge)
						return false;
				}
			}
			return true;
		}
		LineDivision Poly::checkSide(const Line& l, float t) const {
			uint32_t res = 0;
			for(int i=0 ; i<3 ; i++)
				res |= l.checkSide(point[i], t);
			return static_cast<LineDivision>(res);
		}
		bool Poly::isCW() const {
			return (point[1]-point[0]).cw(point[2]-point[0]) >= 0;
		}
		void Poly::invert() {
			std::swap(point[1], point[2]);
		}
		Poly Poly::operator * (const spn::AMat32& m) const {
			Poly ret;
			for(int i=0 ; i<3 ; i++)
				ret.point[i] = point[i].asVec3(1) * m;
			return ret;
		}
		Poly& Poly::operator += (const Vec2& ofs) {
			for(auto& v : point)
				v += ofs;
			return *this;
		}
		void Poly::distend(float width, float mindist) {
			Vec2 center = bs_getCenter();
			for(auto& v : point) {
				auto tmp = v - center;
				float dist = tmp.normalize();
				dist = std::max(mindist, dist+width);
				v = center + tmp * dist;
			}
		}
		std::ostream& operator << (std::ostream& os, const Poly& p) {
			return os << "Polygon(2d) [ 0: " << p.point[0] << std::endl
						<< "1: " << p.point[1] << std::endl
						<< "2: " << p.point[2] << ']';
		}
	}
}
