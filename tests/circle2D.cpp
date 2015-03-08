#ifdef WIN32
	#include <intrin.h>
#endif
#include "spinner/tests/test.hpp"
#include "geom2D.hpp"

namespace boom {
	namespace test {
		using namespace spn::test;
		using spn::Vec2;
		using geo2d::PointM;
		using geo2d::CircleM;
		using geo2d::GSimplex;

		class Circle2D : public spn::test::RandomTestInitializer {};
		TEST_F(Circle2D, Hit_Point) {
			auto rd = getRand();
			auto rc = [&](){ return rd.template getUniform<float>({-1e4f,1e4f}); };
			auto rv = [&](){ return GenRVec<2,false>(rc); };

			CircleM c(rv(), rd.template getUniform<float>({0, 1e3f}));
			PointM p(rv());

			bool b0 = c.hit(p);
			bool b1 = GSimplex(c, p).getResult();
			ASSERT_EQ(b0, b1);
		}
		TEST_F(Circle2D, Hit_Circle) {
			auto rd = getRand();
			auto rc = [&](){ return rd.template getUniform<float>({-1e4f,1e4f}); };
			auto rv = [&](){ return GenRVec<2,false>(rc); };

			CircleM c0(rv(), rd.template getUniform<float>({0, 1e3f})),
					c1(rv(), rd.template getUniform<float>({0, 1e3f}));

			bool b0 = c0.hit(c1);
			bool b1 = GSimplex(c0,c1).getResult();
			ASSERT_EQ(b0, b1);
		}
		TEST_F(Circle2D, Support) {
			auto rd = getRand();
			auto rc = [&](){ return rd.template getUniform<float>({-1e3f,1e3f}); };
			auto rc11 = [&](){ return rd.template getUniform<float>({-1,1}); };
			auto rv = [&](){ return GenRVec<2,false>(rc); };

			CircleM c(rv(), rd.template getUniform<float>({0, 1e3f}));
			Vec2 dir(GenRDir<2,false>(rc11));
			Vec2 sv = c.support(dir);

			// 中心座標からの距離は全てradiusと等しい筈
			EXPECT_NEAR((sv - c.vCenter).length(), c.fRadius, 1e-3f);
			// 方向ベクトルとの内積は限りなくradiusに近い筈
			EXPECT_NEAR((sv - c.vCenter).dot(dir), c.fRadius, 1e-3);
		}
	}
}
