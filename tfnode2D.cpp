#include "geom2D.hpp"

namespace boom {
	namespace geo2d {
		// --------------- TfBase ---------------
		void TfBase::onChildAdded(const SP& /*node*/) {}
		void TfBase::onChildRemove(const SP& /*node*/) {}
		void TfBase::onParentChange(const SP& /*from*/, const SP& /*to*/) {}
		void TfBase::_setAsChanged() {}
		void TfBase::setAsChanged() {
			_setAsChanged();
			if(auto sp = getParent())
				sp->setAsChanged();
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
		spn::Pose2D& TfBase::tf_refPose() {
			Assert(Trap, false, "invalid function call")
			throw 0;
		}
		const spn::Pose2D& TfBase::tf_getPose() const {
			return const_cast<TfBase*>(this)->tf_getPose();
		}
		Model_SP TfBase::im_clone() const {
			return clone();
		}
	}
}
