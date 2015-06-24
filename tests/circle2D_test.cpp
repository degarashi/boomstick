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
		using geo2d::CircleM;
		using geo2d::Circle;
		using geo2d::GSimplex;

		class Circle2D : public spn::test::RandomTestInitializer {};
		TEST_F(Circle2D, Hit_Point) {
			auto rd = getRand();
			auto rdf = rd.template getUniformF<float>();

			CircleM c(GenRCircle(rdf));
			PointM p(GenRPoint(rdf));

			bool b0 = c.hit(p);
			bool b1 = GSimplex(c, p).getResult();
			ASSERT_EQ(b0, b1);
		}
		TEST_F(Circle2D, Hit_Circle) {
			auto rd = getRand();
			auto rdf = rd.template getUniformF<float>();

			CircleM c0(GenRCircle(rdf)),
					c1(GenRCircle(rdf));

			bool b0 = c0.hit(c1);
			bool b1 = GSimplex(c0,c1).getResult();
			ASSERT_EQ(b0, b1);
		}
		TEST_F(Circle2D, Hit_Segment) {
			auto rd = getRand();
			auto rdf = rd.template getUniformF<float>();

			CircleM c(GenRCircle(rdf));
			SegmentM s(GenRSegment(rdf));

			bool b0 = c.hit(s);
			bool b1 = GSimplex(c,s).getResult();
			ASSERT_EQ(b0, b1);
		}
		TEST_F(Circle2D, Support) {
			auto rd = getRand();
			auto rdf = rd.template getUniformF<float>();

			CircleM c(GenRCircle(rdf, {-1e3f,1e3f}));
			Vec2 dir(Vec2::RandomDir(rdf));
			Vec2 sv = c.support(dir);

			// 中心座標からの距離は全てradiusと等しい筈
			EXPECT_NEAR((sv - c.vCenter).length(), c.fRadius, 1e-3f);
			// 方向ベクトルとの内積は限りなくradiusに近い筈
			EXPECT_NEAR((sv - c.vCenter).dot(dir), c.fRadius, 1e-3);
		}
		TEST_F(Circle2D, Boundary) {
			auto rd = getRand();
			auto rdf = rd.template getUniformF<float>();

			int n = rd.template getUniform<int>({1,64});
			std::vector<CircleM> cv;
			for(int i=0 ; i<n ; i++)
				cv.push_back(GenRCircle(rdf));
			// ランダムな個数の円で境界(代表)円を算出
			CircleM cbound(*MakeBoundaryPtr<Circle, geo2d::IModel>(cv.data(), n, sizeof(CircleM), 0));

			for(int i=0 ; i<0x100 ; i++) {
				CircleM ctest = GenRCircle(rdf);
				bool b0 = ctest.hit(cbound),
					 b1 = false;
				for(int j=0 ; j<n ; j++) {
					if((b1 = ctest.hit(cv[j])))
						break;
				}
				// 0b00: 境界円がFalseなら個別判定もFalse
				// 0b01: 不正
				// 0b10: 境界円がTrueで個別判定がFalseはあり得る
				// 0b11: 個別判定がTrueなら境界判定もTrue
				uint32_t flag = (b0 ? 0b10 : 0b00) | (b1 ? 0b01 : 0b00);
				ASSERT_NE(flag, 0b01);
			}
		}
	}
}
