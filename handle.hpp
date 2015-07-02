#pragma once
#include "spinner/resmgr.hpp"

namespace boom {
	namespace geo2d {
		class TfBase;
		using TfBase_SP = std::shared_ptr<TfBase>;
		class TfLeafBase;
		using TfLeaf_SP = std::shared_ptr<TfLeafBase>;
		class TfNode;
		class TfMgr;
		using TfNode_SP = std::shared_ptr<TfNode>;
		DEF_AHANDLE(TfMgr, Tf, TfNode_SP, TfNode_SP)

		struct IModel;
		using Model_SP = std::shared_ptr<IModel>;
	}
	namespace geo3d {
		struct IModel;
		using Model_SP = std::shared_ptr<IModel>;
	}
}

