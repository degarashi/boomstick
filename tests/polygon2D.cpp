#ifdef WIN32
	#include <intrin.h>
#endif
#include "spinner/tests/test.hpp"
#include "geom2D.hpp"

namespace boom {
	namespace test {
		using namespace spn::test;
		using spn::Vec2;
		using geo2d::Poly;
		class Poly2D : public spn::test::RandomTestInitializer {};
		TEST_F(Poly2D, Triangle) {
			auto rd = getRand();
			auto rc = [&](){ return rd.template getUniform<float>({0,1}); };
			auto rv = [&](){ return GenRVec<2,false>(rc); };

			Poly p{rv(), rv(), rv()};
			if(!p.isCW())
				p.invert();

			// 3頂点を合計1になるような係数で合成した座標は必ず三角形の中に入る
			float c[3];
			c[0] = rc();
			c[1] = std::max(.0f, 1.f - c[0]);
			c[2] = std::max(.0f, 1.f - c[0] - c[1]);
			ASSERT_NEAR(c[0]+c[1]+c[2], 1.f, 1e-5f);
			Vec2 pos = spn::MixVector(c, p.point[0], p.point[1], p.point[2]);
			EXPECT_TRUE(p.hit(pos));

			// 頂点順が逆ならば結果も逆になる(辺上の場合は常にflaseだから問題なし)
			p.invert();
			EXPECT_FALSE(p.isInTriangle(pos));
		}
	}
}
