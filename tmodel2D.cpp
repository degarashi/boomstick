#include "geom2D.hpp"

namespace boom {
	namespace geo2d {
// 		TModel::TModel(HMdl hMdl, HTf hTf):
// 				_hlMdl(hMdl),
// 				_hlTf(hTf)
// 		{}
// 		Circle TModel::im_getBVolume() const {
// 			auto bv = _hlMdl->get()->im_getBVolume();
// 			bv.vCenter = toWorld(bv.vCenter);
// 			const auto& m = tm_getToWorld().convertA22();
// 			bv.fRadius *= m.calcDeterminant();
// 			return bv;
// 		}
// 		float TModel::im_getInertia() const {
// 			const auto& m = tm_getToWorld().convertA22();
// 			return _hlMdl->get()->im_getInertia() * m.calcDeterminant();
// 		}
// 		float TModel::im_getArea() const {
// 			return _hlMdl->get()->im_getArea();
// 		}
// 		Vec2 TModel::im_getCenter() const {
// 			return toWorld(_hlMdl->get()->im_getCenter());
// 		}
// 		Vec2 TModel::im_support(const Vec2& dir) const {
// 			return toWorldDir(_hlMdl->get()->im_support(toLocalDir(dir)));
// 		}
// 		bool TModel::im_hitPoint(const Vec2& p, float t) const {
// 			return _hlMdl->get()->im_hitPoint(toLocal(p), t);
// 		}
// 		Vec2 TModel::toLocal(const Vec2& v) const {
// 			return v.asVec3(1) * *getToLocal();
// 		}
// 		Vec2 TModel::toLocalDir(const Vec2& v) const {
// 			return v.asVec3(0) * *getToLocal();
// 		}
// 		Vec2 TModel::toWorld(const Vec2& v) const {
// 			return v.asVec3(1) * *getToWorld();
// 		}
// 		Vec2 TModel::toWorldDir(const Vec2& v) const {
// 			return v.asVec3(0) * *getToWorld();
// 		}
// 		const AMat32& TModel::tm_getToLocal() const {
// 			return _mToLocal = _hlTf->get()->getPose().getToLocal();
// 		}
// 		const AMat32& TModel::tm_getToWorld() const {
// 			return _mToWorld = _hlTf->get()->getPose().getToWorld();
// 		}
// 		spn::Optional<const AMat32&> TModel::getToLocal() const {
// 			return tm_getToLocal();
// 		}
// 		spn::Optional<const AMat32&> TModel::getToWorld() const {
// 			return tm_getToWorld();
// 		}
// 		uint32_t TModel::getCID() const {
// 			return _hlMdl->get()->getCID();
// 		}
// 		const void* TModel::getCore() const {
// 			return _hlMdl->get()->getCore();
// 		}
	}
}
