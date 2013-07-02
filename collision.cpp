#include "geom2D.hpp"

namespace boom {
	bool geo2d::NarrowC_Model::operator()(const geo2d::IModel& mdl0, const geo2d::IModel& mdl1) {
		auto c0 = mdl0.getBCircle(),
			c1 = mdl1.getBCircle();
		// collision check by bounding circle
		if(c0.hit(c1)) {
			// collision check by GJK algorithm
			GSimplex gs(mdl0, mdl1);
			if(gs.getResult()) {
				_inner = gs.getInner();
				return true;
			}
		}
		return false;
	}
	// -------------------------- RigidCR --------------------------
	geo2d::RigidCR::RigidCR(const RCoeff& c): _coeff(c) {}
	void geo2d::RigidCR::operator ()(int id0, int id1, const Rigid& r0, const Rigid& r1, const Vec2& inner) {
		//TODO: 衝突平面のちゃんとした計算
		StLineCore st(Vec2(0,0), Vec2(0,1));
		auto f = CalcForce(r0, r1, inner, _coeff, st);
		setFlagWithInfo(id0, id1, f);
	}
}
