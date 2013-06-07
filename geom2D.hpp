#pragma once
#include "spinner/vector.hpp"
#include <cassert>
#include <vector>

namespace spn {
	namespace Geo2D {
		class Poly {
			mutable AVec2	_center;
			mutable float	_area,
							_inertia;

			Vec2			_point[3];
			enum REFLAG {
				RFL_CENTER = 0x01,
				RFL_AREA = 0x02,
				RFL_INERTIA = 0x04,
				RFL_ALL = 0xff
			};
			mutable uint32_t	_rflag;

			public:
				Poly();
				Poly(const Vec2& p0, const Vec2& p1, const Vec2& p2);
				Poly(const Poly& p) = default;

				float getInertia() const;
				float getArea() const;
				const AVec2& getCenter() const;

				void setPoint(int n, const Vec2& v);
				void addOffset(const Vec2& ofs);

				static float CalcArea(const Vec2& p0, const Vec2& p1, const Vec2& p2);
				static float CalcArea(const Vec2& p0, const Vec2& p1);
		};
		class Circle {
			AVec2			_center;
			float			_radius;
			mutable float	_area;
			mutable bool	_rflag;

			public:
				Circle();
				Circle(const Circle& c) = default;
				Circle(const Vec2& c, float r);

				float area() const;

				const AVec2& getCenter() const;
				float getRadius() const;
				void setCenter(const Vec2& v);
				void setRadius(float r);
		};
		class Convex {
			public:
				using PointL = std::vector<Vec2>;
				using AreaL = std::vector<float>;

			private:
				mutable AVec2	_center;
				mutable float	_area,
								_inertia;

				PointL			_point;
				enum REFLAG {
					RFL_AREA = 0x01,
					RFL_CENTER_INERTIA = 0x02,
					RFL_ALL = 0xff
				};
				mutable uint32_t _rflag;
				void _refreshCInertia() const;

				struct AreaSum {
					float result;

					AreaSum(): result(0) {}
					void operator()(int n, const Vec2& p0, const Vec2& p1) {
						result += Poly::CalcArea(p0,p1);
					}
				};
				struct AreaList {
					AreaL	areaL;
					float*	ptr;
					float	sum;

					AreaList(int n): areaL(n), ptr(&areaL[0]), sum(0) {}
					void operator()(int n, const Vec2& p0, const Vec2& p1) {
						float a = Poly::CalcArea(p0,p1);
						*ptr++ = a;
						sum += a;
					}
				};

			public:
				Convex();
				Convex(const Convex& cnv) = default;
				Convex(const PointL& pl);
				Convex(PointL&& pl);

				float getArea() const;
				float getInertia() const;
				const PointL& getPoint() const;
				PointL& refPoint();

				template <class CB>
				void iterate(CB cb) const {
					// 頂点数は3つ以上
					int nL = _point.size();
					assert(nL > 2);

					// 先にブリッジの箇所を処理
					cb(nL-1, _point.back(), _point.front());
					for(int i=0 ; i<nL-1 ; i++)
						cb(i, _point[i], _point[i+1]);
				}

				const AVec2& center() const;
				float inertia() const;
				void addOffset(const Vec2& ofs);
		};
	}
}
