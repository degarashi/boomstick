#ifdef WIN32
	#include <intrin.h>
#endif
#include "test2D.hpp"

namespace boom {
	namespace test2d {
		using namespace spn::test;
		using spn::Vec2;
		using geo2d::PointM;
		using geo2d::SegmentM;
		using geo2d::AABB;
		using geo2d::AABBM;
		using geo2d::GSimplex;
		using geo2d::CircleM;

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
		TEST_F(AABB2D, Hit_Circle) {
			auto rd = getRand();
			AABBM ab(GenRAABB(rd));
			CircleM c(GenRCircle(rd));

			bool b0 = false;
			// 内包判定(Circle in AABB)
			auto &amin = ab.minV,
				&amax = ab.maxV;
			if(c.support({1,0}).x <= amax.x
				&& c.support({-1,0}).x >= amin.x
				&& c.support({0,1}).y <= amax.y
				&& c.support({0,-1}).y >= amin.y)
				b0 = true;
			// 内包判定(AABB in Circle)
			if(ab.support({1,0}).x <= c.support({1,0}).x
				&& ab.support({-1,0}).x >= c.support({-1,0}).x
				&& ab.support({0,1}).y <= c.support({0,1}).y
				&& ab.support({0,-1}).y >= c.support({0,-1}).y)
				b0 = true;
			// 4辺をセグメント判定
			if(c.hit(SegmentM(amin, Vec2{amin.x, amax.y}))
				|| c.hit(SegmentM(Vec2{amin.x, amax.y}, amax))
				|| c.hit(SegmentM(amax, Vec2{amax.x, amin.y}))
				|| c.hit(SegmentM(Vec2{amax.x, amin.y}, amin)))
				b0 = true;
			bool b1 = GSimplex(ab, c).getResult();
			ASSERT_EQ(b0, b1);
		}
		TEST_F(AABB2D, Boundary) {
			auto rd = getRand();

			int n = rd.template getUniform<int>({1,64});
			std::vector<AABBM> av;
			for(int i=0 ; i<n ; i++)
				av.push_back(GenRAABB(rd));
			// ランダムな個数のAABBで境界(代表)AABBを算出
			AABBM abound(*MakeBoundaryPtr<AABB, geo2d::IModel>(av.data(), n, sizeof(AABBM), 0));

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
