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

		template <class T>
		class TfNode : public Narrow {
			public:
				using Narrow::Narrow;
		};
		using TfNodeTypeList = ::testing::Types<geo2d::Circle,
												geo2d::Segment,
												geo2d::Line,
												geo2d::Ray,
												geo2d::Poly,
												geo2d::AABB>;
		TYPED_TEST_CASE(TfNode, TfNodeTypeList);

		using namespace spn::test;
		// 単一ノードによる姿勢変換テスト
		TYPED_TEST(TfNode, TFNode2D) {
			auto rd = this->getRand();

			using Shape = TypeParam;
			using ShapeM = geo2d::Model<Shape>;
			// 基本の形状に姿勢変換を掛けた物(=A)と
			// 変換後の座標で直接生成した物(=B)の2種類を用意
			ShapeM s;
			test2d::GenRShape(s, rd);
			spn::Pose2D ps(GenR2Vec(rd, {-1e3f, 1e3f}),
							spn::DegF(rd.template getUniform<float>({-180.f, 180.f})),
							GenR2VecAbs(rd, {1e-1f, 1e2f}, {1e-1f, 1e2f}));
			auto spA = std::make_shared<geo2d::TfLeaf<Shape>>(s * ps.getToWorld()),
				 spB = std::make_shared<geo2d::TfLeaf<Shape>>(s);
			spB->setPose(ps);

			int nCheck = 500;
			while(--nCheck >= 0) {
				// 衝突試験用の形状を1つ用意
				ShapeM check;
				test2d::GenRShape(check, rd);
				// 結果は同じになる筈
				using this_t = TfNode<TypeParam>;
				bool bA = this_t::Narrow_t::Hit(spA.get(), &check, 0),
					 bB = this_t::Narrow_t::Hit(spB.get(), &check, 0);
				ASSERT_EQ(bA, bB);
			}
		}
	}
}
