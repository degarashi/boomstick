#include "spinner/vector.hpp"
#include "geom2D.hpp"
#include "spinner/assoc.hpp"
#include <iostream>
using namespace spn;
using namespace boom::geo2d;

struct Pair {
	unsigned int first, second;
	bool operator < (const Pair& p) const {
		return first < p.first;
	}
};
std::ostream& operator << (std::ostream& os, const Pair& p) {
	os << '[' << p.first << ',' << p.second << ']' << std::endl;
	return os;
}

void TestAssoc() {
	AssocVec<Pair> asd;
	std::cout << asd.insert(Pair{100, 0xffff}) << std::endl;
	std::cout << asd.insert(Pair{200, 0xdead}) << std::endl;
	std::cout << asd.insert(Pair{300, 0xfffffff1}) << std::endl;
	std::cout << asd.insert(Pair{50, 0xfffff}) << std::endl;

	asd.erase(0);
	asd.pop_back();

	for(const auto& p : asd)
		std::cout << p;
	std::cout << std::endl;
}
void TestConvex() {
	Convex cnv0{Vec2(0,0), Vec2(0,1), Vec2(1,1), Vec2(1,0)},
			cnv1{Vec2(0,0.5f), Vec2(0.5f,1), Vec2(1,1), Vec2(1,0)};
	cnv1.addOffset(Vec2(1.2f,0));
	GEpa ga(cnv0, cnv1);
	auto np = ga.getNearestPair();
	auto ip = ga.getInner();
}
void TestRigid() {
	RigidMgr rm(IItg::sptr(new itg::Eular));

	IModel::sptr spModel(new ConvexModel({Vec2(0,0), Vec2(0,1), Vec2(1,1), Vec2(1,0)}));
	Rigid::sptr spR(new Rigid(spModel));
	rm.add(spR);
	rm.add(spR);
	rm.add(IResist::sptr(new resist::Gravity(Vec2(0,-9.8f))));
	for(int i=0 ; i<10 ; i++)
		rm.simulate(1.f);
}

int main(int argc, char **argv) {
	auto sz = sizeof(RPose);
	auto sz2 = sz * 128;
	TestRigid();
	return 0;
}
