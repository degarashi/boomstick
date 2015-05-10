#include "tfleaf2D.hpp"

namespace boom {
	namespace geo2d {
		TfLeafBase::TfLeafBase(const Model_SP& m): _model(m) {
			_init();
		}
		void* TfLeafBase::getCore() {
			++refModelAccum();
			return TfBase::getCore();
		}
		const void* TfLeafBase::getCore() const {
			return TfBase::getCore();
		}
		void TfLeafBase::setAsChanged() {
			++refModelAccum();
			TfBase::setAsChanged();
		}
		void TfLeafBase::setModelSource(const Model_SP& m) {
			_model = m;
			_init();
		}
		void TfLeafBase::_init() {
			_rflag.template set<ShapeAccum>(0);
			_rflag.resetAll();
		}
		bool TfLeafBase::isLeaf() const {
			return true;
		}
		void TfLeafBase::im_transform(void* dst, const AMat32& m) const {
			auto m2 = getPose().getToWorld().convertA33() * m;
			_model->im_transform(dst, m2);
		}
		void TfLeafBase::im_getBVolume(Circle& c) const {
			c = getBCircle();
		}
		void TfLeafBase::im_getBVolume(AABB& a) const {
			a = getBBox();
		}
		float TfLeafBase::im_getInertia() const {
			return getInertia();
		}
		float TfLeafBase::im_getArea() const {
			return getArea();
		}
		Vec2 TfLeafBase::im_getCenter() const {
			return getCenter();
		}
		Vec2 TfLeafBase::im_support(const Vec2& dir) const {
			auto& sc = getPose().getScale();
			AssertP(Trap, std::abs(sc.x - sc.y) < 1e-4f);
			return toWorld(_model->im_support(toLocalDir(dir)));
		}
		bool TfLeafBase::im_hitPoint(const Vec2& p, float t) const {
			return _model->im_hitPoint(toLocal(p), t);
		}
		Vec2 TfLeafBase::toLocal(const spn::Vec2& v) const {
			return v.asVec3(1) * getLocal();
		}
		Vec2 TfLeafBase::toLocalDir(const Vec2& v) const {
			return (v.asVec3(0) * getLocal()).normalization();
		}
		Vec2 TfLeafBase::toWorld(const Vec2& v) const {
			return v.asVec3(1) * getGlobal();
		}
		Vec2 TfLeafBase::toWorldDir(const Vec2& v) const {
			return (v.asVec3(0) * getGlobal()).normalization();
		}
		spn::Optional<const AMat32&> TfLeafBase::getToLocal() const {
			return getLocal();
		}
		spn::Optional<const AMat32&> TfLeafBase::getToWorld() const {
			return getGlobal();
		}
		spn::Pose2D& TfLeafBase::tf_refPose() {
			return refPose();
		}
		const spn::Pose2D& TfLeafBase::tf_getPose() const {
			return getPose();
		}
		const Model_SP& TfLeafBase::getModelSource() const {
			return _model;
		}
		std::ostream& TfLeafBase::print(std::ostream& os) const {
			return _model->print(os);
		}
		uint32_t TfLeafBase::getCID() const {
			return _model->getCID();
		}
		std::ostream& operator << (std::ostream& os, const TfLeafBase& node) {
			return os << "TfLeaf2D [ pose: " << node.getPose() << std::endl
						<< "node accum: " << node.getNodeAccum() << ']';
		}
		uint32_t TfLeafBase::_refresh(uint32_t& a, ShapeAccum*) const {
			++a;
			getModelAccum();
			getPose();
			return 0;
		}
		uint32_t TfLeafBase::_refresh(AMat32& m, Global*) const {
			getNodeAccum();
			auto& ps = getPose();
			{
				auto& sc = ps.getScale();
				AssertP(Trap, std::abs(sc.x - sc.y) < 1e-4f);
			}
			ps.getToWorld().convert(m);
			return 0;
		}
		uint32_t TfLeafBase::_refresh(AMat32& m, Local*) const {
			Mat33 tm;
			getGlobal().convert33().inversion(tm);
			tm.convert(m);
			return 0;
		}
		uint32_t TfLeafBase::_refresh(float& d, Determinant*) const {
			d = getGlobal().convertA22().calcDeterminant();
			return 0;
		}
		uint32_t TfLeafBase::_refresh(float& f, Inertia*) const {
			getModelAccum();
			f = _model->im_getInertia() * getDeterminant();
			return 0;
		}
		uint32_t TfLeafBase::_refresh(float& f, Area*) const {
			getModelAccum();
			f = _model->im_getArea() * getDeterminant();
			return 0;
		}
		uint32_t TfLeafBase::_refresh(spn::Vec2& v, Center*) const {
			getModelAccum();
			v = toWorld(_model->im_getCenter());
			return 0;
		}
		uint32_t TfLeafBase::_refresh(Circle& c, BCircle*) const {
			getModelAccum();
			_model->im_getBVolume(c);
			c = c * getPose().getToWorld();
			return 0;
		}
		uint32_t TfLeafBase::_refresh(AABB& a, BBox*) const {
			getModelAccum();
			auto m = getPose().getToLocal();
			auto& scale = getPose().getScale();
			a.maxV.x = toWorld(_model->im_support(toLocalDir({1,0}))).x;
			a.maxV.y = toWorld(_model->im_support(toLocalDir({0,1}))).y;
			a.minV.x = toWorld(_model->im_support(toLocalDir({-1,0}))).x;
			a.minV.y = toWorld(_model->im_support(toLocalDir({0,-1}))).y;
			return 0;
		}
	}
}
