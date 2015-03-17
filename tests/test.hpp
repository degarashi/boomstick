#pragma once
#include "spinner/tests/test.hpp"
#include "test.hpp"

namespace boom {
	namespace test {
		template <int N, class RD>
		std::vector<spn::VecT<N,false>> GenRVectors(RD& rd, const int n, const spn::RangeF& rV, const float threshold_sq) {
			auto rv = [&](){ return spn::test::GenRVec<N,false>(rd, rV); };
			using Vt = spn::VecT<N,false>;
			std::vector<Vt> pl(n);
			for(int i=0 ; i<n ; i++) {
				Vt p;
				bool bLoop;
				do {
					bLoop = false;
					p = rv();
					for(int j=0 ; j<i ; j++) {
						if(pl[j].dist_sq(p) <= threshold_sq) {
							bLoop = true;
							break;
						}
					}
				} while(bLoop);
				pl[i] = p;
			}
			return std::move(pl);
		}
	}
}
