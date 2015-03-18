#include "geom2D.hpp"

namespace boom {
	namespace geo2d {
		// --------------- TfBase ---------------
		void TfBase::onChildAdded(const SP& /*node*/) {
			++this->refNodeAccum();
		}
		void TfBase::onChildRemove(const SP& /*node*/) {
			++this->refNodeAccum();
		}
		uint32_t TfBase::_refresh(spn::AMat33& m, Global*) const {
			getNodeAccum();
			auto& ps = getPose();
			ps.getToWorld().convert(m);
			return 0;
		}
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
		std::ostream& operator << (std::ostream& os, const TfBase& node) {
			return os << "TfNode2D [ pose: " << node.getPose() << std::endl
						<< "node accum: " << node.getNodeAccum() << ']';
		}
	}
}
