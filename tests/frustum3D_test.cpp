#ifdef WIN32
	#include <intrin.h>
#endif
#include "test3D.hpp"

namespace boom {
	namespace test3d {
		using namespace spn::test;
		using geo3d::SphereM;
		using geo3d::FrustumM;
		using geo3d::GSimplex;
		class Frustum3D : public spn::test::RandomTestInitializer {};
		TEST_F(Frustum3D, Hit_AABB) {
			auto rd = getRand();
			auto rdf = rd.template getUniformF<float>();

			auto a = GenRAABB(rdf, {-1e1, 1e1});
			auto f = GenRFrustum(rdf, {-1.f, 1.f}, {1e-1f, 1e1f}, {1e-1f, 1e1f});
			const bool b0 = f.hit(a);
			bool b1 = GSimplex(a, f).getResult();
			if(b0 != b1) {
				// GJKの検査結果にある程度の誤差を許容
				if(b0) {
					a.expand(1e-2f);
					b1 = GSimplex(a,f).getResult();
				} else {
					a.expand(-1e2f);
					b1 = GSimplex(a,f).getResult();
				}
			}
			ASSERT_EQ(b0, b1);
		}
	}
}
