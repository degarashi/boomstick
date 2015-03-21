#ifdef WIN32
	#include <intrin.h>
#endif
#include "test2D.hpp"

namespace boom {
	namespace test {
		using geo2d::PointM;
		class Narrow : public spn::test::RandomTestInitializer {
			protected:
				using base_t = spn::test::RandomTestInitializer;
				using Types = ::boom::geo2d::Types;
				using Narrow_t = Types::Narrow;
				void SetUp() override {
					base_t::SetUp();
					Narrow_t::Initialize();
				}
		};
		// 階層構造を含めたNarrowPheseテスト (2D)
		TEST_F(Narrow, Dim2) {
			auto rd = getRand();

			;
			// TODO: ランダムに当たり判定階層構造を作る
			// とりあえずは単純な階層構造でテスト
			auto fnN = [](){ return std::make_shared<geo2d::TfNode_Static<geo2d::Circle>>(); };
			auto fnL = [&rd](){ return std::make_shared<geo2d::TfLeaf<geo2d::Circle>>(test2d::GenRCircle(rd)); };
			auto fnLen = [&]() { return rd.template getUniform<int>({1, 64}); };

			auto c0 = fnN(),
				 c1 = fnN();
			int n = fnLen();
			for(int i=0 ; i<n ; i++)
				c0->addChild(fnL());
			n = fnLen();
			for(int i=0 ; i<n ; i++)
				c1->addChild(fnL());

			// 個別に判定した結果
			bool b0 = false;
			auto v0 = c0->plainPtr(),
				 v1 = c1->plainPtr();
			v0.erase(v0.begin());
			v1.erase(v1.begin());
			for(auto* p0 : v0) {
				for(auto* p1 : v1) {
					if((b0 |= Narrow_t::Hit(p0, p1)))
						break;
				}
			}
			// ツリー構造で判定した結果
			bool b1 = Narrow_t::Hit(c0.get(), c1.get());
			// 両者は一致する筈
			ASSERT_EQ(b0, b1);
		}
	}
}
