#ifdef WIN32
	#include <intrin.h>
#endif
#include "test3D.hpp"

namespace boom {
	namespace test3d {
		using namespace spn::test;
		using geo3d::SphereM;
		using geo3d::PointM;
		using geo3d::GSimplex;
		class Sphere3D : public spn::test::RandomTestInitializer {};
		TEST_F(Sphere3D, Hit_Point) {
			auto rd = getRand();

			auto s = GenRSphere(rd);
			auto p = GenRPoint(rd);
			bool b0 = s.hit(p);
			bool b1 = GSimplex(s, p).getResult();
			ASSERT_EQ(b0, b1);
		}
	}
}
