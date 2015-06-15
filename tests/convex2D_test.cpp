#ifdef WIN32
	#include <intrin.h>
#endif
#include "test2D.hpp"

namespace boom {
	namespace test2d {
		using namespace spn::test;
		using spn::Vec2;
		using geo2d::PointL;
		using geo2d::PointM;
		using geo2d::Convex;
		using geo2d::ConvexM;
		using geo2d::GSimplex;
		using geo2d::Poly;
		using geo2d::Circle;

		class Convex2D : public spn::test::RandomTestInitializer {};
		TEST_F(Convex2D, BVolume) {
			auto rd = getRand();
			int np = rd.template getUniform<int>({3, 64});
			ConvexM cv = GenRConvex(rd, defval::convex_pos, np);

			Circle c = cv.bs_getBCircle();
			// 全ての点が円の中に入ってるかチェック
			for(auto& p : cv.point) {
				ASSERT_LE(c.vCenter.distance(p), c.fRadius + NEAR_THRESHOLD);
			}
		}
		TEST_F(Convex2D, CheckPosition) {
			auto rd = getRand();
			auto rcf = [&](){ return rd.template getUniform<float>({0, 1.f}); };

			int np = rd.template getUniform<int>({3, 64});
			ConvexM cv = GenRConvex(rd, defval::convex_pos, np);
			np = cv.getNPoints();
			ASSERT_TRUE(cv.checkCW());
			Vec2 center = (cv.getPoint(0) + cv.getPoint(1) + cv.getPoint(2)) / 3.f;
			for(int i=0 ; i<np ; i++) {
				// ポリゴンを構成する3頂点から係数を合計1で算出した座標は必ずcheckPositionでその位置を返す
				Poly poly{center,
						cv.getPoint(i),
						cv.getPoint((i+1)%np)};
				ASSERT_TRUE(poly.isCW());
				float c[3] = {rcf(), rcf(), rcf()};
				float sum = c[0] + c[1] + c[2];
				for(float& f : c)
					f /= sum;
				Vec2 tmp = spn::MixVector(c, poly.point[0], poly.point[1], poly.point[2]);
				auto res = cv.checkPosition(tmp);
				EXPECT_NE(ConvexPos::Outer, res.first);
				EXPECT_EQ(i, res.second);
			}
		}
		TEST_F(Convex2D, Hit_Point) {
			auto rd = getRand();
			auto rv = [&](){ return spn::Vec2::Random(rd, {-1e3f,1e3f}); };

			// 凸包の頂点数(3〜64)
			int np = rd.template getUniform<int>({3, 64});
			ConvexM c = GenRConvex(rd, defval::convex_pos, np);
			PointM p(GenRPoint(rd, {-1e3f, 1e3f}));
			bool b0 = c.hit(p);
			bool b1 = GSimplex(c,p).getResult();
			ASSERT_EQ(b0, b1);
		}
		TEST_F(Convex2D, Hit_Convex) {
			auto rd = getRand();
			auto rnp = [&](){ return rd.template getUniform<int>({3, 64}); };
			// 凸包の頂点数(3〜64)
			ConvexM c0 = GenRConvex(rd, {-1e2f, 1e2f}, rnp()),
					c1 = GenRConvex(rd, {-1e2f, 1e2f}, rnp());
			int nc0 = c0.getNPoints(),
				nc1 = c1.getNPoints();
			// ポリゴンを総当りで地道に(確実に)判定
			bool b0 = false;
			for(int i=0 ; i<nc0-2 ; i++) {
				Poly poly0{c0.getPoint(0), c0.getPoint(i+1), c0.getPoint(i+2)};
				ASSERT_TRUE(poly0.isCW());
				for(int j=0 ; j<nc1-2 ; j++) {
					Poly poly1{c1.getPoint(0), c1.getPoint(j+1), c1.getPoint(j+2)};
					ASSERT_TRUE(poly1.isCW());
					if(poly0.hit(poly1, 1e-5f)) {
						b0 = true;
						i=nc0;
						break;
					}
				}
			}
			// GJKにて判定
			bool b1 = GSimplex(c0, c1).getResult();
			// 結果を比較
			ASSERT_EQ(b0, b1);
		}
	}
}
