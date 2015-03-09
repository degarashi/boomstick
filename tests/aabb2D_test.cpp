#ifdef WIN32
	#include <intrin.h>
#endif
#include "test.hpp"

namespace boom {
	namespace test {
		using namespace spn::test;
		using spn::Vec2;
		using geo2d::PointM;
		using geo2d::SegmentM;
		using geo2d::AABBM;
		using geo2d::GSimplex;

		class AABB2D : public spn::test::RandomTestInitializer {};
		TEST_F(AABB2D, Hit_Point) {
			auto rd = getRand();
			PointM p(GenRPoint(rd));
			AABBM ab(GenRAABB(rd));
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
			spn::RangeF rV{-1e3f, 1e3f};
			AABBM ab(GenRAABB(rd));
			SegmentM s(GenRSegment(rd));
			bool b0 = ab.hit(s);
			GSimplex gs(ab, s);
			bool b1 = gs.getResult();
			ASSERT_EQ(b0, b1);
		}
		TEST_F(AABB2D, Hit_AABB) {
			auto rd = getRand();
			AABBM ab0(GenRAABB(rd)),
				  ab1(GenRAABB(rd));
			bool b0 = ab0.hit(ab1);
			GSimplex gs(ab0, ab1);
			bool b1 = gs.getResult();
			ASSERT_EQ(b0, b1);
		}
	}
}
