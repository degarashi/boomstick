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
	// 頂点を定義して、重心を求めそこを中心にして座標を変換
	ConvexModel* c0 = ConvexModel::New({Vec2(0,0), Vec2(0,1), Vec2(1,1), Vec2(1,0)});
	ConvexModel* c1 = ConvexModel::New({Vec2(0.5f,0.5f), Vec2(3,1), Vec2(3,0)});
	AVec2 ofs0 = c0->getCenter(),
		ofs1 = c1->getCenter();
	c0->addOffset(-ofs0);
	c1->addOffset(-ofs1);

	Rigid r0((IModel::sptr(c0))),
			r1((IModel::sptr(c1)));
	RPose &p0 = r0.refPose(),
			&p1 = r1.refPose();
	p0.setOfs(ofs0);
	p1.setOfs(ofs1);
	p0.setVelocity(Vec2(0,1));
	p1.setVelocity(Vec2(0,-1));

	GSimplex gs(r0, r1);
	auto bHit = gs.getResult();
	auto inner = gs.getInner();
	if(bHit) {
		RCoeff coeff{1,1,1,1,1};
		auto rf = CalcForce(r0, r1, inner, coeff, StLineCore(Vec2(1,1), Vec2(0,1)));
		std::cout << rf << std::endl;
	}
	return;
//	GEpa ga(cnv0, cnv1);
//	auto np = ga.getNearestPair();
//	auto ip = ga.getInner();
}
void TestRigid() {
	RigidMgr rm(IItg::sptr(new itg::Eular));

	IModel::sptr spModel(new ConvexModel({Vec2(0,0), Vec2(0,1), Vec2(1,1), Vec2(1,0)}));
	Rigid::sptr spR(new Rigid(spModel));
	rm.add(spR);
	rm.add(spR);
	spR->addR(IResist::sptr(new resist::Gravity(Vec2(0,-9.8f))));
	RCoeff rc = {1,1,1,1,1};
	spR->addR(IResist::sptr(new resist::Impact(rc)));
	for(int i=0 ; i<10 ; i++)
		rm.simulate(1.f);
}

int main(int argc, char **argv) {
	TestConvex();
//	TestRigid();
	return 0;
}
