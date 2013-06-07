#include "geom2D.hpp"

namespace spn {
	namespace Geo2D {
		// ---------------------- Poly ----------------------
		Poly::Poly(): _rflag(RFL_ALL) {}
		Poly::Poly(const Vec2& p0, const Vec2& p1, const Vec2& p2): _point{p0,p1,p2}, _rflag(RFL_ALL) {}

		float Poly::getArea() const {
			if(_rflag & RFL_AREA) {
				_rflag &= ~RFL_AREA;
				_area = CalcArea(_point[0], _point[1], _point[2]);
			}
			return _area;
		}
		const AVec2& Poly::getCenter() const {
			if(_rflag & RFL_CENTER) {
				_rflag &= ~RFL_CENTER;
				_center = (_point[0] + _point[1] + _point[2]) * (1.0f/3);
			}
			return _center;
		}
		void Poly::setPoint(int n, const Vec2& v) {
			_point[n] = v;
			_rflag = RFL_ALL;
		}
		void Poly::addOffset(const Vec2& ofs) {
			for(int i=0 ; i<3 ; i++)
				_point[i] += ofs;
		}
		float Poly::getInertia() const {
			if(_rflag & RFL_INERTIA) {
				_rflag &= ~RFL_INERTIA;

				_inertia = (1.0f/18) * (_point[0].dot(_point[0])
										+ _point[0].dot(_point[0])
										+ _point[0].dot(_point[0])
										- _point[1].dot(_point[2])
										- _point[2].dot(_point[0])
										- _point[0].dot(_point[1]));
			}
			return _inertia;
		}

		float Poly::CalcArea(const Vec2& p0, const Vec2& p1, const Vec2& p2) {
			return  (p1-p0).ccw(p2-p0) * 0.5f;
		}
		float Poly::CalcArea(const Vec2& p0, const Vec2& p1) {
			return p0.ccw(p1) * 0.5f;
		}

		// ---------------------- Circle ----------------------
		Circle::Circle(): _rflag(true) {}
		Circle::Circle(const Vec2& c, float r): _center(c), _radius(r), _rflag(true) {}
		float Circle::area() const {
			if(_rflag) {
				_rflag = false;

				_area = 0;
			}
			return _area;
		}
		const AVec2& Circle::getCenter() const { return _center; }
		float Circle::getRadius() const { return _radius; }
		void Circle::setCenter(const Vec2& v) {
			// 面積は変化しない
			_center = v;
		}
		void Circle::setRadius(float r) {
			_rflag = true;
			_radius = r;
		}

		// ---------------------- Convex ----------------------
		Convex::Convex(): _rflag(RFL_ALL) {}
		Convex::Convex(const PointL& pl): _point(pl), _rflag(RFL_ALL) {}
		Convex::Convex(PointL&& pl): _point(pl), _rflag(RFL_ALL) {}

		const Convex::PointL& Convex::getPoint() const { return _point; }
		Convex::PointL& Convex::refPoint() {
			_rflag = RFL_ALL;
			return _point;
		}
		float Convex::getArea() const {
			if(_rflag & RFL_AREA) {
				_rflag &= ~RFL_AREA;

				AreaSum as;
				iterate(as);
				_area = as.result;
			}
			return _area;
		}
		void Convex::_refreshCInertia() const {
			if(_rflag & RFL_CENTER_INERTIA) {
				_rflag &= ~RFL_CENTER_INERTIA;

				int nL = _point.size();
				AreaList al(nL);
				iterate(al);
				al.sum = _sseRcp22Bit(al.sum);
				float areaInv3 = al.sum * (1.0f/3),
						areaInv6 = al.sum * (1.0f/6);
				_center *= 0;
				_inertia = 0;
				iterate([this, &al, areaInv3, areaInv6](int n, const Vec2& p0, const Vec2& p1) {
					_center += al.areaL[n] * areaInv3 * (p0 + p1);
					_inertia += al.areaL[n] * areaInv6 * (p0.dot(p0) + p0.dot(p1) + p1.dot(p1));
				});
				_inertia -= _center.dot(_center);
			}
		}
		const AVec2& Convex::center() const {
			_refreshCInertia();
			return _center;
		}
		float Convex::getInertia() const {
			_refreshCInertia();
			return _inertia;
		}
	}
}
