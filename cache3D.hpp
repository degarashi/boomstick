#pragma once
#include "cache.hpp"

namespace boom {
	namespace geo3d {
		using spn::Vec3;
		using spn::Mat44;

		//! CacheTag: area
		DEF_TAG(TagArea, float, area)
		//! CacheTag: inertia
		DEF_TAG(TagInertia, Mat44, inertia)
		//! CacheTag: center
		DEF_TAG(TagCenter, Vec3, center)

		struct SphereCore;
		//! CacheTag: bounding volume(sphere)
		DEF_TAG(TagBSphere, SphereCore, bsphere)
		#undef DEF_TAG

		#define DEF_CORE_FUNCS3D(...) using CT = CType<__VA_ARGS__>; \
		template <class T> \
		using Wrap = CoreRaw<T>; \
		float area() const; \
		Mat44 inertia() const; \
		Vec3 center() const; \
		SphereCore bsphere() const;
	}
}
