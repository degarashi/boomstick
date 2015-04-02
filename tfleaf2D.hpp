#pragma once
#include "geom2D.hpp"

namespace boom {
	namespace geo2d {
		// スケーリングはX,Yとも同じ比率のみ許可
		template <class Shape, class Ud=spn::none_t>
		class TfLeafBase : public spn::CheckAlign<16, TfLeafBase<Shape,Ud>>,
							public TfBase,
							public Model<Shape>
		{
			protected:
				using model_t = Model<Shape>;
				Ud					_udata;
				#define SEQ_TFLEAF \
					((Pose)(spn::Pose2D)) \
					((NodeAccum)(uint32_t)) \
					((Global)(spn::AMat32)(Pose)(NodeAccum)) \
					((Local)(spn::AMat32)(Global)) \
					((Determinant)(float)(Global)) \
					((ModelAccum)(uint32_t)) \
					((Inertia)(float)(ModelAccum)(Determinant)) \
					((Area)(float)(ModelAccum)(Determinant)) \
					((Center)(spn::Vec2)(ModelAccum)(Global)) \
					((GCenter)(spn::Vec2)(ModelAccum)(Global)) \
					((BVolume)(Circle)(ModelAccum)(Determinant)(Global)) \
					((ShapeAccum)(uint32_t)(ModelAccum)(Global))
				RFLAG_S(TfLeafBase, SEQ_TFLEAF)
			private:
				void _init() {
					_rflag.template set<ShapeAccum>(0);
					_rflag.resetAll();
				}

			public:
				RFLAG_SETMETHOD_S(SEQ_TFLEAF)
				RFLAG_GETMETHOD_S(SEQ_TFLEAF)
				RFLAG_REFMETHOD_S(SEQ_TFLEAF)
				#undef SEQ_TFLEAF
				TfLeafBase() { _init(); }
				TfLeafBase(const model_t& m): model_t(m) { _init(); }
				TfLeafBase(model_t&& m): model_t(std::move(m)) { _init(); }

				/*! ユーザーデータがvoidの時は親ノードのデータを返す */
				void* getUserData() override {
					return _getUserData(&_udata, std::is_same<spn::none_t, Ud>());
				}
				bool isLeaf() const override {
					return true;
				}
				Circle im_getBVolume() const override { return getBVolume(); }
				float im_getInertia() const override { return getInertia(); }
				float im_getArea() const override { return getArea(); }
				Vec2 im_getCenter() const override { return getCenter(); }
				Vec2 im_support(const Vec2& dir) const override {
					auto& sc = getPose().getScale();
					AssertP(Trap, std::abs(sc.x - sc.y) < 1e-4f);
					return toWorld(model_t::support(toLocalDir(dir)));
				}
				bool im_hitPoint(const Vec2& p, float t=NEAR_THRESHOLD) const override {
					return model_t::hit(toLocal(p), t);
				}
				Vec2 toLocal(const spn::Vec2& v) const override {
					return v.asVec3(1) * getLocal();
				}
				Vec2 toLocalDir(const Vec2& v) const override {
					return (v.asVec3(0) * getLocal()).normalization();
				}
				Vec2 toWorld(const Vec2& v) const override {
					return v.asVec3(1) * getGlobal();
				}
				Vec2 toWorldDir(const Vec2& v) const override {
					return (v.asVec3(0) * getGlobal()).normalization();
				}
				spn::Optional<const AMat32&> getToLocal() const override {
					return getLocal();
				}
				spn::Optional<const AMat32&> getToWorld() const override {
					return getGlobal();
				}
				friend std::ostream& operator << (std::ostream&, const TfLeafBase&);
		};
		template <class Shape, class Ud>
		std::ostream& operator << (std::ostream& os, const TfLeafBase<Shape,Ud>& node) {
			return os << "TfLeaf2D [ pose: " << node.getPose() << std::endl
						<< "node accum: " << node.getNodeAccum() << ']';
		}

		template <class Shape, class Ud>
		uint32_t TfLeafBase<Shape,Ud>::_refresh(uint32_t& a, ShapeAccum*) const {
			++a;
			getModelAccum();
			getGlobal();
			return 0;
		}
		template <class Shape, class Ud>
		uint32_t TfLeafBase<Shape,Ud>::_refresh(AMat32& m, Global*) const {
			getNodeAccum();
			auto& ps = getPose();
			{
				auto& sc = ps.getScale();
				AssertP(Trap, std::abs(sc.x - sc.y) < 1e-4f);
			}
			ps.getToWorld().convert(m);
			return 0;
		}
		template <class Shape, class Ud>
		uint32_t TfLeafBase<Shape,Ud>::_refresh(AMat32& m, Local*) const {
			Mat33 tm;
			getGlobal().convert33().inversion(tm);
			tm.convert(m);
			return 0;
		}
		template <class Shape, class Ud>
		uint32_t TfLeafBase<Shape,Ud>::_refresh(float& d, Determinant*) const {
			d = getGlobal().convertA22().calcDeterminant();
			return 0;
		}
		template <class Shape, class Ud>
		uint32_t TfLeafBase<Shape,Ud>::_refresh(float& f, Inertia*) const {
			getModelAccum();
			f = model_t::bs_getInertia() * getDeterminant();
			return 0;
		}
		template <class Shape, class Ud>
		uint32_t TfLeafBase<Shape,Ud>::_refresh(float& f, Area*) const {
			getModelAccum();
			f = model_t::bs_getArea() * getDeterminant();
			return 0;
		}
		template <class Shape, class Ud>
		uint32_t TfLeafBase<Shape,Ud>::_refresh(spn::Vec2& v, Center*) const {
			getModelAccum();
			v = toWorld(model_t::bs_getCenter());
			return 0;
		}
		template <class Shape, class Ud>
		uint32_t TfLeafBase<Shape,Ud>::_refresh(spn::Vec2& v, GCenter*) const {
			getModelAccum();
			v = toWorld(model_t::bs_getGCenter());
			return 0;
		}
		template <class Shape, class Ud>
		uint32_t TfLeafBase<Shape,Ud>::_refresh(Circle& c, BVolume*) const {
			getModelAccum();
			c = model_t::bs_getBVolume();
			auto& ps = getPose();
			c.vCenter = toWorld(c.vCenter);
			auto m = ps.getToWorld().convertA22();
			c.fRadius *= getDeterminant();
			return 0;
		}

		// Convex以外
		template <class Shape, class Ud=spn::none_t>
		class TfLeaf : public spn::CheckAlign<16, TfLeaf<Shape,Ud>>,
						public TfLeafBase<Shape, Ud>
		{
			private:
				using base_t = TfLeafBase<Shape, Ud>;
				using ModelAccum = typename base_t::ModelAccum;
				using Global = typename base_t::Global;

				mutable uint32_t _shapeAccum = 0;
				mutable Shape	_tfShape;
				const Shape& getTfShape() const {
					auto accum = base_t::getShapeAccum();
					if(accum != _shapeAccum) {
						_shapeAccum = accum;
						_tfShape = *this * base_t::getGlobal();
					}
					return _tfShape;
				}
			public:
				using base_t::base_t;

				Vec2 im_support(const Vec2& dir) const override {
					// 予め変換しておいた形状でサポート写像
					return getTfShape().support(dir);
				}
				const void* getCore() const override {
					getTfShape();
					return &_tfShape;
				}
				void* getCore() override {
					getTfShape();
					return &_tfShape;
				}
				typename base_t::SP clone() const override {
					return std::make_shared<TfLeaf>(*this);
				}
		};

		//! 座標変換ありのモデル基底
		template <class Ud>
		class TfLeaf<Convex, Ud> : public spn::CheckAlign<16, TfLeaf<Convex,Ud>>,
								public TfLeafBase<Convex, Ud>
		{
			private:
				using base_t = TfLeafBase<Convex, Ud>;
			public:
				using base_t::base_t;
				bool canCacheShape() const override {
					return false;
				}
				const void* getCore[[noreturn]]() const override {
					// canCacheShapeがfalseのオブジェクトに対して呼んではいけない
					throw std::domain_error("calling getCore() is not valid for TfLeaf<Convex>");
				}
				void* getCore[[noreturn]]() override {
					throw std::domain_error("calling getCore() is not valid for TfLeaf<Convex>");
				}
				typename base_t::SP clone() const override {
					return std::make_shared<TfLeaf>(*this);
				}
		};
	}
}
