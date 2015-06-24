#ifdef WIN32
	#include <intrin.h>
#endif
#include "test2D.hpp"

namespace boom {
	namespace test2d {
		using namespace spn::test;
		using geo2d::CapsuleM;
		using geo2d::PointM;
		using geo2d::GSimplex;

		class Capsule2D : public spn::test::RandomTestInitializer {};
		TEST_F(Capsule2D, Hit_Point) {
			auto rd = getRand();
			auto rdf = rd.template getUniformF<float>();

			CapsuleM c(GenRCapsule(rdf));
			PointM p(GenRPoint(rdf));

			bool b0 = c.hit(p);
			bool b1 = GSimplex(c, p).getResult();
			ASSERT_EQ(b0, b1);
			bool b2 = c.im_hitPoint(p);
			ASSERT_EQ(b0, b2);
		}
		TEST_F(Capsule2D, Hit_Capsule) {
			auto rd = getRand();
			auto rdf = rd.template getUniformF<float>();

			auto range = spn::RangeF{0, 1e2f};
			CapsuleM c0(GenRCapsule(rdf, range, range)),
					 c1(GenRCapsule(rdf, range, range));
			bool b0 = c0.hit(c1);
			bool b1 = GSimplex(c0, c1).getResult();
			geo2d::Segment s0(c0.from, c0.to),
							s1(c1.from, c1.to);
			float d = s0.distance(s1);
			float r = c0.radius + c1.radius;
			if(b0 != b1) {
				auto dir0 = (c0.to-c0.from).normalization();
				auto dir1 = (c1.to-c1.from).normalization();
				float dd = dir0.dot(dir1);
				bool b0 = c0.hit(c1);
				bool b1 = GSimplex(c0, c1).getResult();
			}
			ASSERT_EQ(b0, b1);
		}
	}
}
