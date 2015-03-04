#pragma once
#include "spinner/resmgr.hpp"

namespace boom {
	namespace geo2d {
		class TfNode;
		class TfMgr;
		using TfNode_SP = std::shared_ptr<TfNode>;
		DEF_AHANDLE(TfMgr, Tf, TfNode_SP, TfNode_SP)

		struct IModel;
		using UPModel = std::unique_ptr<IModel>;
		class ModelMgr;
		DEF_AHANDLE(ModelMgr, Mdl, UPModel, UPModel)
	}
	namespace geo3d {
		struct IModel;
		using UPModel = std::unique_ptr<IModel>;
		class ModelMgr;
		DEF_AHANDLE(ModelMgr, Mdl, UPModel, UPModel)
	}
}

