#pragma once
#include "cache.hpp"

namespace boom {
	namespace geo2d {
		using spn::Vec2;

		//! CacheTag: area
		DEF_TAG(TagArea, float, area)
		//! CacheTag: inertia
		DEF_TAG(TagInertia, float, inertia)
		//! CacheTag: center
		DEF_TAG(TagCenter, Vec2, center)

		struct CircleCore;
		//! CacheTag: bounding volume(circle)
		DEF_TAG(TagBCircle, CircleCore, bcircle)

		#undef DEF_TAG

		#define DEF_CORE_FUNCS(...) using CT = CType<__VA_ARGS__>; \
			template <class T> \
			using Wrap = CoreRaw<T>; \
			float area() const; \
			float inertia() const; \
			Vec2 center() const; \
			CircleCore bcircle() const;
	}
}