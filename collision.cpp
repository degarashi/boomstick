#include "geom2D.hpp"
#include "rigid2D.hpp"

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
}
