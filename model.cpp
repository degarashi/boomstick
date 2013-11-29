#include "geom2D.hpp"

namespace boom {
	const AMat22 cs_mRot90[2] = {
		{spn::COS90, -spn::SIN90,
		spn::SIN90, spn::COS90},
		{spn::COS90, spn::SIN90,
		-spn::SIN90, spn::COS90}
	};
	float Area_x2(const Vec2& v0, const Vec2& v1) {
		return std::fabs(v0.cw(v1));
	}
	float Area_x2(const Vec3& v0, const Vec3& v1) {
		return v0.cross(v1).length();
	}
	// --------------------- IModelNode ---------------------
	void IModelNode::notifyChange() {
		if(pParent)
			pParent->notifyChange();
	}
	void IModelNode::applyChange() {}
	MdlIP IModelNode::getInner() const { return MdlIP(); }
	bool IModelNode::hasInner() const { return false; }
	void* IModelNode::getUserData() const { return nullptr; }

	namespace geo2d {
/*
		// -------------------------- TModel --------------------------
		template <class MDL, class BASE>
		TModel<MDL,BASE>::TModel(const MDL &mdl): _model(mdl) {}

		#define DEF_TMODEL(typ)	template <class MDL,class BASE> typ TModel<MDL,BASE>
		DEF_TMODEL(const IModel&)::_getModel(const IModel& mdl) const { return mdl; }
		DEF_TMODEL(const IModel&)::_getModel(const HLMdl& mdl) const { return *mdl.cref().get(); }
		DEF_TMODEL(Vec2)::support(const Vec2& dir) const {
			// dirをローカルに座標変換してサポート写像した後、またワールド座標系に戻す
			Vec2 pos = _getModel(_model).support(BASE::toLocalDir(dir));
			return BASE::toWorld(pos);
		}
		DEF_TMODEL(Vec2)::getCenter() const {
			Vec2 res = _getModel(_model).getCenter();
			return BASE::toWorld(res);
		}
		DEF_TMODEL(bool)::isInner(const Vec2& pos) const {
			Vec2 tpos = BASE::toLocal(pos);
			return _getModel(_model).isInner(tpos);
		}
		DEF_TMODEL(const MDL&)::getModel() const {
			return _model;
		}
		DEF_TMODEL(uint32_t)::getCID() const {
			return _getModel(_model).getCID();
		}
		DEF_TMODEL(Circle)::getBCircle() const {
			Circle c = _getModel(_model).bs_getBCircle();
			return c * BASE::getToWorld();
		}
		DEF_TMODEL(float)::getArea(bool bInv) const {
			return _getModel(_model).getArea(bInv);
		}
		DEF_TMODEL(float)::getInertia(bool bInv) const {
			// 返すのは常に重心周りの慣性モーメントなので座標変換は不要
			return _getModel(_model).getInertia(bInv);
		}
		DEF_TMODEL(IModel::PosL)::getOverlappingPoints(const IModel& mdl, const Vec2& inner) const {
			// innerとmdlをこちらのローカル座標系に変換
			TModelR<TR_Mat> t_mdl(mdl, (const BASE&)*this, TR_Mat::TagInverse);
			Vec2 t_inner(BASE::toLocal(inner));
			// 結果をグローバル座標系へ変換
			auto ret = _getModel(_model).getOverlappingPoints(t_mdl, t_inner);
			const auto& mat = BASE::getToWorld();
			for(auto& p : ret.second)
				p = p.asVec3(1) * mat;
			return std::move(ret);
		}
		DEF_TMODEL(int)::getNPoints() const { return _getModel(_model).getNPoints(); }
		DEF_TMODEL(Vec2)::getPoint(int n) const {
			Vec2 p = _getModel(_model).getPoint(n);
			return BASE::toWorld(p);
		}
		DEF_TMODEL(IModel::CPos)::checkPosition(const Vec2& pos) const {
			// posをローカル座標系に変換
			Vec2 tpos = BASE::toLocal(pos);
			return _getModel(_model).checkPosition(tpos);
		}
		DEF_TMODEL(Convex2)::splitTwo(const Line& ls) const {
			Line tline(BASE::toLocal(ls.pos),
						BASE::toLocalDir(ls.dir));
			auto res = _getModel(_model).splitTwo(tline);
			const auto& mat = BASE::getToWorld();
			res.first *= mat;
			res.second *= mat;
			return std::move(res);
		}
		DEF_TMODEL(std::ostream&)::dbgPrint(std::ostream& os) const {
			int nP = getNPoints();
			if(nP > 0) {
				for(int i=0 ; i<nP-1 ; i++)
					os << getPoint(i) << std::endl;
				os << getPoint(nP-1);
			}
			return os;
		}
		#undef DEF_TMODEL
		template class TModel<HLMdl, TR_Mat>;
		template class TModel<const IModel&, TR_Mat>;*/
	}
}
