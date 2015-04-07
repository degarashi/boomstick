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
		using Epa2DTypeList = ::testing::Types<geo2d::Poly,
												geo2d::Circle,
												geo2d::AABB>;
		TYPED_TEST_CASE(Epa2D, Epa2DTypeList);
		namespace {
			const spn::RangeF c_rangeV{-1e2f, 1e2f},
								c_rangeR{0, 1e2f};
			template <class RD>
			void GenRShapeEPA(geo2d::CircleM& c, RD& rd) {
				test2d::GenRShape(c, rd, c_rangeV, c_rangeR);
			}
			template <class S, class RD>
			void GenRShapeEPA(S& s, RD& rd) {
				test2d::GenRShape(s, rd, c_rangeV);
			}
			template <class RD>
			void GenRShapeEPA(geo2d::ConvexM& c, RD& rd) {
				test2d::GenRShape(c, rd, c_rangeV);
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
				c0.distend(-ErrorAdjust);
				GEpa gepa2(c0, c1, ErrorAdjust);
				ASSERT_TRUE(!gepa2.getResult() ||
							gepa2.getPVector().length() < ErrorAdjust);
			} else {
				auto np = gepa.getNearestPair();
				// 左辺の物体を最近傍対ベクトル分(p1 - p0)移動させたら衝突する筈
				auto pv = (np.second - np.first);
				c0 += pv;
				c0.distend(ErrorAdjust);
				GEpa gepa2(c0, c1, ErrorAdjust);
				bool b = gepa2.getResult();
				if(!b) {
					auto np2 = gepa2.getNearestPair();
					float d2 = np2.first.distance(np2.second);
					// 誤差の関係で割と範囲を広めにとる
					ASSERT_LT(d2, 1e-1);
				}
			}
		}
	}
}
