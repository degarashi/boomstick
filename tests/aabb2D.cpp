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
		using geo2d::SegmentM;
		using geo2d::AABBM;
		using geo2d::GSimplex;

		namespace {
			template <class RV>
			AABBM GenRAABB(RV& rv) {
				Vec2 v0 = rv(),
					 v1 = rv(),
					 tmp = v0;
				v0.selectMin(v1);
				v1.selectMax(tmp);
				return AABBM(v0, v1);
			}
		}
		class AABB2D : public spn::test::RandomTestInitializer {};
		TEST_F(AABB2D, Hit_Point) {
			auto rd = getRand();
			auto rc = [&](){ return rd.template getUniform<float>({-1e4f,1e4f}); };
			auto rv = [&](){ return GenRVec<2,false>(rc); };
			PointM p(rv());
			AABBM ab(GenRAABB(rv));
			ASSERT_LE(ab.minV.x, ab.maxV.x);
			ASSERT_LE(ab.minV.y, ab.maxV.y);
			// AABB -> Point のHit関数をGJK関数で結果を比較
			bool b0 = ab.hit(p);
			GSimplex gs(ab, p);
			bool b1 = gs.getResult();
			ASSERT_EQ(b0, b1);
		}
		TEST_F(AABB2D, Hit_Segment) {
			auto rd = getRand();
			auto rc = [&](){ return rd.template getUniform<float>({-1e3f,1e3f}); };
			auto rv = [&](){ return GenRVec<2,false>(rc); };
			AABBM ab(GenRAABB(rv));
			SegmentM s(rv(), rv());
			bool b0 = ab.hit(s);
			GSimplex gs(ab, s);
			bool b1 = gs.getResult();
			ASSERT_EQ(b0, b1);
		}
		TEST_F(AABB2D, Hit_AABB) {
			auto rd = getRand();
			auto rc = [&](){ return rd.template getUniform<float>({-1e3f,1e3f}); };
			auto rv = [&](){ return GenRVec<2,false>(rc); };
			AABBM ab0(GenRAABB(rv)),
				  ab1(GenRAABB(rv));
			bool b0 = ab0.hit(ab1);
			GSimplex gs(ab0, ab1);
			bool b1 = gs.getResult();
			ASSERT_EQ(b0, b1);
		}
	}
}
