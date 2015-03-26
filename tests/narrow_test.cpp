#ifdef WIN32
	#include <intrin.h>
#endif
#include "test2D.hpp"

namespace boom {
	namespace test2d {
		void Narrow::SetUp() {
			base_t::SetUp();
			Narrow_t::Initialize();
		}
	}
	namespace test {
		using geo2d::PointM;
		using Narrow = test2d::Narrow;
		// 階層構造を含めたNarrowPheseテスト (2D)
		TEST_F(Narrow, RandomTree2D) {
			auto rd = getRand();

			// ランダムに当たり判定階層構造を作る
			auto c0 = test2d::MakeRandomTree(rd, 64, 8);
			auto c1 = test2d::MakeRandomTree(rd, 64, 8);
			auto v0 = test2d::CollectLeaf(c0),
				 v1 = test2d::CollectLeaf(c1);
			// 個別に判定した結果
			bool b0 = false;
			for(auto* p0 : v0) {
				for(auto* p1 : v1) {
					if((b0 |= Narrow_t::Hit(p0, p1, 0)))
						break;
				}
			}
			// ツリー構造で判定した結果
			bool b1 = Narrow_t::Hit(c0.get(), c1.get(), 0);
			// 両者は一致する筈
			ASSERT_EQ(b0, b1);
		}
	}
}
