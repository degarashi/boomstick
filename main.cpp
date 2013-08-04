#include "spinner/vector.hpp"
#include "rigid2D.hpp"
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
	std::unique_ptr<RigidRes> rmgr(new RigidRes);

	boom::RCoeff rc = {1,1,1,1,1};
	RigidMgr rm(IItg::sptr(new itg::Eular));
	rm.setCoeff(rc);
	// 頂点を定義して、重心を求めそこを中心にして座標を変換
	constexpr float fc0 = 0.4f,
					fc1 = 2.0f;
	auto	c0 = ConvexModel::NewUF({Vec2(0,0), Vec2(0,fc0), Vec2(fc0,fc0), Vec2(fc0,0)}),
			c1 = ConvexModel::NewUF({Vec2(-fc1,0), Vec2(fc1,0), Vec2(fc1,-fc1), Vec2(-fc1,-fc1)});
	AVec2 ofs0 = c0->getCenter(),
		ofs1 = c1->getCenter();
	c0->addOffset(-ofs0);
	c1->addOffset(-ofs1);
	HLMdl lhC0 = mgr_rigid.acquireModel(std::move(c0)),
			lhC1 = mgr_rigid.acquireModel(std::move(c1));
	HLRig lhR[2] = {mgr_rigid.acquireRigid(Rigid::NewUF(lhC0)),
					mgr_rigid.acquireRigid(Rigid::NewUF(lhC1))};
	RPose* pp[2] = {&lhR[0].ref()->refPose(),
					&lhR[1].ref()->refPose()};
	pp[0]->setOfs(ofs0);
	pp[1]->setOfs(ofs1);
	pp[0]->setVelocity(Vec2(1,0));
	pp[1]->setVelocity(Vec2(0,0));

	auto spGrav = IResist::sptr(new resist::Gravity(Vec2(0,-4.f)));
	auto spCol = IResist::sptr(new resist::Impact);
	for(int i=0 ; i<2 ; i++) {
		rm.addA(lhR[i].get());
		auto& r = *lhR[i].ref();
		r.addR(spGrav, 0xdead);
		r.addR(spGrav);
		r.addR(spCol);
		r.remR(0xdead);
		r.remR(Rigid::DEFAULT_ID);
	}

	for(int i=0 ; i<10 ; i++)
		rm.simulate(0.01f);
}
int main(int argc, char **argv) {
	TestRigid();
	return 0;
}
