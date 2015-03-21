#include "geom2D.hpp"

namespace boom {
	namespace geo2d {
		// --------------- TfBase ---------------
		void TfBase::onChildAdded(const SP& /*node*/) {}
		void TfBase::onChildRemove(const SP& /*node*/) {}
		void TfBase::onParentChange(const SP& /*from*/, const SP& /*to*/) {}
		void* TfBase::_getUserData(void*, std::true_type) {
			if(auto sp = getParent())
				return sp->getUserData();
			return nullptr;
		}
		void* TfBase::_getUserData(void* udata, std::false_type) {
			return udata;
		}
		MdlItr TfBase::getInner() const {
			return MdlItr(getChild());
		}
		bool TfBase::hasInner() const {
			return static_cast<bool>(getChild());
		}

		// --------------- TfLeaf_base ---------------
		uint32_t TfLeaf_base::_refresh(spn::AMat33& m, Global*) const {
			getNodeAccum();
			auto& ps = getPose();
			ps.getToWorld().convert(m);
			return 0;
		}
		std::ostream& operator << (std::ostream& os, const TfLeaf_base& node) {
			return os << "TfLeaf2D [ pose: " << node.getPose() << std::endl
						<< "node accum: " << node.getNodeAccum() << ']';
		}
	}
}
