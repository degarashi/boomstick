#ifdef WIN32
	#include <intrin.h>
#endif
#include "spinner/tests/test.hpp"
#include "spinner/size.hpp"
#include "geom2D.hpp"

namespace boom {
	namespace test {
		using namespace spn::test;
		using spn::Vec2;
		using geo2d::PointM;
		using geo2d::SegmentM;
		using geo2d::GSimplex;

		template <int N, bool A, class RD>
		auto MakeRVec(RD& rd, const spn::RangeF& s) {
			auto rc = [&rd, s](){ return rd.template getUniform<float>(s); };
			return [rc](){
				return GenRVec<N, A>(rc);
			};
		}
		template <int N, bool A, class RD>
		auto MakeRDir(RD& rd) {
			auto rc = [&rd](){ return rd.template getUniform<float>({-1.f, 1.f}); };
			return [rc](){
				return GenRDir<N,A>(rc);
			};
		}
		class Segment2D : public spn::test::RandomTestInitializer {};
		TEST_F(Segment2D, Hit_Point) {
			auto rd = getRand();
			auto rc = [&rd](){ return rd.template getUniform<float>({0, 1.f}); };
			auto rv = MakeRVec<2,false>(rd, {-1e3f, 1e3f});

			SegmentM s{rv(), rv()};
			PointM v(rv());
			// 点との判定がGJKと一致するか
			auto chk = [&s,&v](){
				bool b0 = s.hit(v);
				GSimplex gs(s,v);
				bool b1 = gs.getResult();
				EXPECT_EQ(b0, b1);
			};
			chk();

			// 両端の点が同じケース
			s.to = s.from;
			chk();

			// 線分上の点は必ずヒットする
			float c[2];
			c[0] = rc();
			c[1] = 1.f-c[0];
			if(s.from.dist_sq(s.to) < geo2d::NEAR_THRESHOLD_SQ)
				static_cast<Vec2&>(v) = s.from;
			else
				static_cast<Vec2&>(v) = s.from*c[0] + s.to*c[1];
			EXPECT_TRUE(s.hit(v));
		}
		TEST_F(Segment2D, Nearest_Point) {
			auto rd = getRand();
			auto rv = MakeRVec<2,false>(rd, {-1e3f, 1e3f});
			SegmentM s{rv(), rv()};
			PointM p(rv());
			constexpr auto NS = geo2d::NEAR_THRESHOLD_SQ;
			// 最近傍点がきちんと線分上に乗っているかのテスト
			auto res = s.nearest(p);
			switch(res.second) {
				case LinePos::OnLine:
					EXPECT_TRUE(s.hit(res.first));
					break;
				case LinePos::Begin:
					EXPECT_LE(s.from.dist_sq(res.first), NS);
					break;
				case LinePos::End:
					EXPECT_LE(s.to.dist_sq(res.first), NS);
					break;
				default:
					EXPECT_TRUE(false);
			}
		}
		TEST_F(Segment2D, Hit_Segment) {
			// 線分との判定結果をGJKと比較
			auto rd = getRand();
			auto rv = MakeRVec<2,false>(rd, {-1e3f, 1e3f});
			SegmentM s0{rv(), rv()},
					 s1{rv(), rv()};
			bool b0 = s0.hit(s1);
			GSimplex gs(s0, s1);
			bool b1 = gs.getResult();
			if(b0 != b1) {
				bool b0 = s0.hit(s1);
				GSimplex gs(s0, s1);
				bool b1 = gs.getResult();
			}
			EXPECT_EQ(b0, b1);
		}
		TEST_F(Segment2D, Support) {
			auto rd = getRand();
			auto rv = MakeRVec<2,false>(rd, {-1e3f, 1e3f});
			auto rdir = MakeRDir<2,false>(rd);

			SegmentM s{rv(), rv()};
			auto dir = rdir();
			auto sv = s.support(dir);
			float dist0 = sv.dist_sq(s.from),
				  dist1 = sv.dist_sq(s.to),
				  dot0 = dir.dot(s.from),
				  dot1 = dir.dot(s.to);

			// サポート座標は線分の両端のどちらかと一致している筈
			if(dist0 < geo2d::NEAR_THRESHOLD_SQ) {
				// 更に両端の内、内積が大きい方と一致しているか
				ASSERT_GE(dot0, dot1);
			} else if(dist1 < geo2d::NEAR_THRESHOLD_SQ) {
				ASSERT_GE(dot1, dot0);
			} else {
				ASSERT_TRUE(false);
			}
		}
	}
}
