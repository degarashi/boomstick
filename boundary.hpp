#pragma once

#include "geom.hpp"
namespace boom {
	template <class S, class IM>
	spn::Optional<S> MakeBoundary(const IM** p, size_t n, typename IModelNode::Time_t t) {
		S s;
		auto fp = &S::setBoundary;
		while(n-- > 0) {
			if((*p)->imn_refresh(t)) {
				(s.*fp)(*p);
				fp = &S::appendBoundary;
			}
			++p;
		}
		if(fp == &S::setBoundary)
			return spn::none;
		return s;
	}
	template <class S, class IM>
	spn::Optional<S> MakeBoundaryPtr(const void* pp, size_t n, size_t stride, typename IModelNode::Time_t t) {
		S s;
		auto fp = &S::setBoundary;
		auto pv = reinterpret_cast<uintptr_t>(pp);
		while(n-- > 0) {
			auto* p = reinterpret_cast<const IM*>(pv);
			if(p->imn_refresh(t)) {
				(s.*fp)(p);
				fp = &S::appendBoundary;
			}
			pv += stride;
		}
		if(fp == &S::setBoundary)
			return spn::none;
		return s;
	}
}
