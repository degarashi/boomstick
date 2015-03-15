#ifdef WIN32
	#include <intrin.h>
#endif
#include "test2D.hpp"

namespace boom {
	namespace test {
		using namespace spn::test;
		using spn::Vec2;
		using geo2d::PointM;
		using geo2d::GSimplex;
		class Point2D : public spn::test::RandomTestInitializer {};
		TEST_F(Point2D, Hit_Point) {
			auto rd = getRand();

			PointM p0(GenRPoint(rd)),
				   p1(GenRPoint(rd));
			// Point -> Point の、Hit関数とGJK関数で結果が一致するかチェック
			bool b0 = p0.hit(p1);
			GSimplex gs(p0, p1);
			bool b1 = gs.getResult();
			ASSERT_EQ(b0, b1);
		}
	}
}
