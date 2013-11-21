//! 3D形状同士の当たり判定に関する実装 (hit関数)
#include "geom3D.hpp"

namespace boom {
	namespace geo3d {
// 		Narrow::CFunc Narrow::GetFunc(const IModel* mdl0, const IModel* mdl1) {
// 			return cs_cfunc[(mdl0->getCID() << 4) | mdl1->getCID()];
// 		}
// 		bool Narrow::Hit(const IModel* mdl0, const IModel* mdl1) {
// 			auto cf = GetCFunc(mdl0, mdl1);
// 			if(cf(mdl0->getCore(), mdl1->getCore()))
// 				return hitL(mdl0, mdl1);
// 			return false;
// 		}
// 		bool Narrow::HitL(const IModel* mdl0, const IModel* mdl1, bool bFinal) {
// 			auto in = mdl0->getInner();
// 			if(in.first != in.second) {
// 				auto cf = GetCFunc(in.first, mdl1);
// 				do {
// 					if(cf(in.first, mdl1))
// 						return hitL(mdl1, in.first, false);
// 				} while(++in.first != in.second);
// 			}
// 			if(bFinal)
// 				return false;
// 			return hitL(mdl1, mdl0, true);
// 		}
	}
}
