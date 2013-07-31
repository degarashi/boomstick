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
			decltype(std::declval<T>().area()) get(const T& src) {
				return src.area();
			}
			void get(...) {}
		};
		//! CacheTag: inertia
		struct TagInertia : TagBase<float> {
			template <class T>
			decltype(std::declval<T>().inertia()) get(const T& src) {
				return src.inertia();
			}
			void get(...) {}
		};
		//! CacheTag: center
		struct TagCenter : TagBase<Vec2> {
			template <class T>
			decltype(std::declval<T>().center()) get(const T& src) {
				return src.center();
			}
			void get(...) {}
		};
		struct CircleCore;
		//! CacheTag: bounding volume(circle)
		struct TagBCircle : TagBase<CircleCore> {
			template <class T>
			decltype(std::declval<T>().bcircle()) get(const T& src) {
				return src.bcircle();
			}
			void get(...) {}
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