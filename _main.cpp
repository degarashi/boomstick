//! テストコードの記述
#include "geom3D.hpp"
#include "geom2D.hpp"
#include "bc_roundrobin.hpp"
#include "convex.hpp"
using namespace spn;
int main() {
// 	boom::geo3d::Narrow na;
// 	boom::geo3d::Narrow::Initialize();
// 	boom::geo3d::SphereM sp;
// 	boom::geo3d::CapsuleM cp;
// 	boom::geo3d::Narrow::Hit(&sp, &cp);
// 	boom::geo2d::Narrow na2;
// 	boom::geo2d::Narrow::Initialize();
// 	boom::geo2d::CircleM m;

	boom::geo3d::SphereM m(Vec3(0,0,0), 10.f);
	boom::geo3d::SegmentM s(Vec3(-1,0,0), Vec3(-1,1,0));
	auto vv = boom::geo3d::Polygon(Vec3(0,0,0), Vec3(0,1,0), Vec3(1,0,0)).nearest(Vec3(0,0,0)).first;
	boom::geo3d::Types::Narrow narrow;
	decltype(narrow)::Initialize();
	boom::geo3d::GSimplex gs(m, s);
	auto inp = gs.getInterPoint();
	boom::geo3d::EConvex ec(gs, m, s);
	auto av = ec.getAvoidVector();
	av.first *= av.second;

	{
		boom::geo2d::CircleM m({0,0}, 1.f);
		boom::geo2d::AABBM a({0,1},{1,2});
		boom::geo2d::Types::Narrow narrow;
		decltype(narrow)::Initialize();
		boom::geo2d::GSimplex gs(m, a);
		gs.getResult();
		gs.getInner();
	}
/*
	boom::geo3d::ModelMgr mmgr;
	boom::geo3d::HLMdl hlMdl = mmgr.emplace(new boom::geo3d::SphereM(spn::Vec3(1,2,3), 0));
	boom::ColMgr<boom::BroadC_RoundRobin,
				boom::geo3d::Types,
				uint32_t> col;
	auto ch = col.addCol(0x01, hlMdl);
	auto ch2 = col.addCol(0x01, hlMdl);
	auto ch3 = col.addCol(0x01, hlMdl);
	col.update();
	std::cout << "-------------------------------" << std::endl;
	ch.ref().getCollision([](const decltype(col)::Hist& h) {
		std::cout << h.hCol.getIndex() << ',' << h.nFrame << std::endl;
	});
	ch.ref().getEndCollision([](const decltype(col)::Hist& h) {
		std::cout << "End:" << std::hex << h.hCol.getIndex() << ',' << h.nFrame << std::endl;
	});
	ch2.release();
	ch3.release();

	std::cout << "-------------------------------" << std::endl;
	col.update();
	ch.ref().getCollision([](const decltype(col)::Hist& h) {
		std::cout << h.hCol.getIndex() << ',' << h.nFrame << std::endl;
	});
	ch.ref().getEndCollision([](const decltype(col)::Hist& h) {
		std::cout << "End:" << h.hCol.getIndex() << ',' << h.nFrame << std::endl;
	});*/

	return 0;
}
