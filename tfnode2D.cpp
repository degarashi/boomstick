#include "geom2D.hpp"

namespace boom {
	namespace geo2d {
		void TfNode::onParentChange(const SP&, const SP&) {
			++this->refNodeAccum();
		}
		uint32_t TfNode::_refresh(spn::AMat33& m, Global*) const {
			getNodeAccum();
			auto& ps = getPose();
			ps.getToWorld().convert(m);
			return 0;
		}
		std::ostream& operator << (std::ostream& os, const TfNode& node) {
			return os << "TfNode2D [ pose: " << node.getPose() << std::endl
						<< "node accum: " << node.getNodeAccum() << ']';
		}
	}
}
