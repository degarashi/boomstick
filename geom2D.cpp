#include "geom2D.hpp"

namespace boom {
	namespace geo2d {
		bool IsCW(const PointL& pts) {
			int nP = pts.size();
			if(nP < 3)
				return false;
			for(int i=0 ; i<nP-2 ; i++) {
				if((pts[i+1] - pts[i]).cw(pts[i+2] - pts[i+1]) < 0)
					return false;
			}
			return (pts[0] - pts[nP-1]).cw(pts[1] - pts[0]) >= 0;
		}
		bool IsCrossing(const Line& ls0, const Line& ls1, float len0, float len1, float t) {
			bool res0 = false,
				 res1 = false;
			auto fn0 = [&res0, len0, t](float f){
				res0 = spn::IsInRange(f, -t, len0+t);
				return f;
			};
			auto fn1 = [&res1, len1, t](float f){
				res1 = spn::IsInRange(f, -t, len1+t);
				return f;
			};
			NearestPoint(ls0, ls1, fn0, fn1);
			return res0 & res1;
		}

		// ------------- MdlItr -------------
		MdlItr::MdlItr(const TfBase_SP& sp): _sp(sp) {}
		MdlItr& MdlItr::operator ++ () {
			_sp = _sp->getSibling();
			return *this;
		}
		bool MdlItr::operator == (const MdlItr& m) const {
			return _sp == m._sp;
		}
		bool MdlItr::operator != (const MdlItr& m) const {
			return !(operator == (m));
		}
		MdlItr::operator bool () const {
			return static_cast<bool>(_sp);
		}
		const TfBase* MdlItr::get() const {
			return _sp.get();
		}
	}
}
