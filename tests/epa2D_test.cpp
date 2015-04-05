#ifdef WIN32
	#include <intrin.h>
#endif
#include "test2D.hpp"

namespace boom {
	namespace test2d {
		using namespace spn::test;
		using geo2d::GSimplex;
		using geo2d::GEpa;
		class Epa2D : public spn::test::RandomTestInitializer {};
		TEST_F(Epa2D, Epa) {
			auto rd = getRand();

			const spn::RangeF rv{-1e3f,1e3f},
								rr{0, 1e2f};
			auto c0 = test2d::GenRCircle(rd, rv, rr),
				 c1 = test2d::GenRCircle(rd, rv, rr);
			float dist = c0.vCenter.distance(c1.vCenter);
			float radius = c0.fRadius+c1.fRadius;
			LogOutput("%1%\n%2%\ndistance=%3%\nradius=%4%", c0, c1, dist, radius);
			GEpa gepa(c0, c1, 1e-4f);
			constexpr float ErrorAdjust = 5e-3f;
			if(gepa.getResult()) {
				auto pv = gepa.getPVector();
				LogOutput("Hit: %1% length=%2%, diff=%3%", pv, pv.length(), dist-radius+pv.length());
				// 左辺に指定した物体を衝突回避ベクトル分移動させたら衝突回避出来る筈
				c0.vCenter += pv;
				float kusoge = c0.vCenter.distance(c1.vCenter);
				// 誤差の関係で確実に判定するために半径を調整
				c0.fRadius *= 1.f - ErrorAdjust;
				bool b0 = c0.hit(c1);
				ASSERT_FALSE(GSimplex(c0, c1).getResult());
			} else {
				auto np = gepa.getNearestPair();
				float dist2 = np.first.distance(np.second);
				LogOutput("NotHit: %1% -> %2% length=%3%, diff=%4%", np.first, np.second,
						dist2, dist-radius-dist2);
				// 左辺の物体を最近傍対ベクトル分(p1 - p0)移動させたら衝突する筈
				c0.vCenter += (np.second - np.first);
				// 誤差の関係で確実に判定するために半径を調整
				c0.fRadius *= 1.f + ErrorAdjust;
				ASSERT_TRUE(GSimplex(c0, c1).getResult());
			}
		}
	}
}
