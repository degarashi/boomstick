#pragma once
#include "cache.hpp"

namespace boom {
	namespace geo2d {
		using spn::Vec2;

		struct TagArea;
		struct TagInertia;
		struct TagCenter;
		struct TagBCircle;

		template <class Dat>
		struct TagBase {
			using Data = Dat;
			Dat operator()();
		};

		//! CacheTag: area
		struct TagArea : TagBase<float> {
			template <class T>
			decltype(T().area()) get(const T& src) {
				return src.area();
			}
		};
		//! CacheTag: inertia
		struct TagInertia : TagBase<float> {
			template <class T>
			decltype(T().inertia()) get(const T& src) {
				return src.inertia();
			}
		};
		//! CacheTag: center
		struct TagCenter : TagBase<Vec2> {
			template <class T>
			decltype(T().center()) get(const T& src) {
				return src.center();
			}
		};
		struct CircleCore;
		//! CacheTag: bounding volume(circle)
		struct TagBCircle : TagBase<CircleCore> {
			template <class T>
			decltype(T().bcircle()) get(const T& src) {
				return src.bcircle();
			}
		};

		#define DEF_CORE_FUNCS(...) using CT = CType<__VA_ARGS__>; \
			template <class T> \
			using Wrap = CoreRaw<T>; \
			float area() const; \
			float inertia() const; \
			Vec2 center() const; \
			CircleCore bcircle() const;
	}
}