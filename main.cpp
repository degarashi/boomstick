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
void TestRigid() {
	boom::RCoeff rc = {1,1,1,1,1};
	RigidMgr rm(IItg::sptr(new itg::Eular), rc);
	// 頂点を定義して、重心を求めそこを中心にして座標を変換
	ConvexModel *c0 = ConvexModel::New({Vec2(0,0), Vec2(0,1), Vec2(1,1), Vec2(1,0)}),
				*c1 = ConvexModel::New(*c0);
	c1->addOffset(Vec2(0.75f));
	AVec2 ofs0 = c0->getCenter(),
		ofs1 = c1->getCenter();
	c0->addOffset(-ofs0);
	c1->addOffset(-ofs1);

	Rigid::sptr rg[2] = {Rigid::sptr(Rigid::New(IModel::sptr(c0))),
						Rigid::sptr(Rigid::New(IModel::sptr(c1)))};
	RPose* pp[2] = {&rg[0]->refPose(),
					&rg[1]->refPose()};
	pp[0]->setOfs(ofs0);
	pp[1]->setOfs(ofs1);
	pp[0]->setVelocity(Vec2(1,0));
	pp[1]->setVelocity(Vec2(-1,0));

	auto spGrav = IResist::sptr(new resist::Gravity(Vec2(0,-9.8f)));
	auto spCol = IResist::sptr(new resist::Impact);
	for(int i=0 ; i<2 ; i++) {
		rm.addA(rg[i]);
		rg[i]->addR(spGrav);
		rg[i]->addR(spCol);
	}
	for(int i=0 ; i<10 ; i++)
		rm.simulate(1.f);
}

int main(int argc, char **argv) {
	TestRigid();
	return 0;
}
