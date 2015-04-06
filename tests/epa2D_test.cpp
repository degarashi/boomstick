#ifdef WIN32
	#include <intrin.h>
#endif
#include "test2D.hpp"

namespace boom {
	namespace test2d {
		using namespace spn::test;
		using geo2d::GSimplex;
		using geo2d::GEpa;
		template <class T>
		class Epa2D : public spn::test::RandomTestInitializer {};
		using Epa2DTypeList = ::testing::Types<geo2d::Circle,
												geo2d::AABB>;
		TYPED_TEST_CASE(Epa2D, Epa2DTypeList);
		namespace {
			template <class RD>
			void GenRShapeEPA(geo2d::CircleM& c, RD& rd) {
				test2d::GenRShape(c, rd, spn::RangeF{-1e2f, 1e2f}, spn::RangeF{0, 1e2f});
			}
			template <class RD>
			void GenRShapeEPA(geo2d::AABBM& a, RD& rd) {
				test2d::GenRShape(a, rd, spn::RangeF{-1e2f, 1e2f});
			}
		}

		TYPED_TEST(Epa2D, Epa) {
			auto rd = this->getRand();

			using ShapeM = geo2d::Model<TypeParam>;
			ShapeM c0, c1;
			GenRShapeEPA(c0, rd);
			GenRShapeEPA(c1, rd);
			constexpr float ErrorAdjust = 5e-3f;
			GEpa gepa(c0, c1, ErrorAdjust);
			if(gepa.getResult()) {
				auto pv = gepa.getPVector();
				// 左辺に指定した物体を衝突回避ベクトル分移動させたら衝突回避出来る筈
				c0 += pv;
				c0.distend(-ErrorAdjust*2);
				ASSERT_FALSE(GSimplex(c0, c1).getResult());
			} else {
				auto np = gepa.getNearestPair();
				// 左辺の物体を最近傍対ベクトル分(p1 - p0)移動させたら衝突する筈
				auto pv = (np.second - np.first);
				c0 += pv;
				c0.distend(ErrorAdjust*2);
				ASSERT_TRUE(GSimplex(c0, c1).getResult());
			}
		}
	}
}
