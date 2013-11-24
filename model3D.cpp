#include "geom3D.hpp"

namespace boom {
	namespace geo3d {
		// --------------------- IModel ---------------------
		void IModel::notifyChange() {
			if(pParent)
				pParent->notifyChange();
		}
		void IModel::applyChange() {}
		MdlIP IModel::getInner() const { return MdlIP(); }
		void* IModel::getUserData() const { return nullptr; }
	}
}
