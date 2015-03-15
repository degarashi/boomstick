#ifdef WIN32
	#include <intrin.h>
#endif
#include "test2D.hpp"

namespace boom {
	namespace test {
		using namespace spn::test;
		using spn::Vec2;
		using geo2d::PointM;
		using geo2d::SegmentM;
		using geo2d::AABB;
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
		TEST_F(AABB2D, Boundary) {
			auto rd = getRand();

			int n = rd.template getUniform<int>({1,64});
			std::vector<AABBM> av;
			for(int i=0 ; i<n ; i++)
				av.push_back(GenRAABB(rd));
			// ランダムな個数のAABBで境界(代表)AABBを算出
			AABBM abound(AABB::Boundary(av.data(), n, sizeof(AABBM)));

			for(int i=0 ; i<0x100 ; i++) {
				AABBM atest = GenRAABB(rd);
				bool b0 = atest.hit(abound),
					b1 = false;
				for(int j=0 ; j<n ; j++) {
					if((b1 = atest.hit(av[j])))
						break;
				}
				// 0b00: 境界AABBがFalseなら個別判定もFalse
				// 0b01: 不正
				// 0b10: 境界AABBがTrueで個別判定がFalseはあり得る
				// 0b11: 個別判定がTrueなら境界判定もTrue
				uint32_t flag = (b0 ? 0b10 : 0b00) | (b1 ? 0b01 : 0b00);
				ASSERT_NE(flag, 0b01);
			}
		}
	}
}
