#pragma once
#include "geom2D.hpp"

namespace boom {
	namespace geo2d {
		// スケーリングはX,Yとも同じ比率のみ許可
		class TfLeafBase : public spn::CheckAlign<16, TfLeafBase>,
							public TfBase
		{
			protected:
				Model_SP			_model;
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
					((BCircle)(Circle)(ModelAccum)(Pose)) \
					((BBox)(AABB)(ModelAccum)(Pose)) \
					((ShapeAccum)(uint32_t)(ModelAccum)(Pose))
				RFLAG_S(TfLeafBase, SEQ_TFLEAF)
			private:
				void _init();
			public:
				RFLAG_SETMETHOD_S(SEQ_TFLEAF)
				RFLAG_GETMETHOD_S(SEQ_TFLEAF)
				RFLAG_REFMETHOD_S(SEQ_TFLEAF)
				#undef SEQ_TFLEAF
				TfLeafBase(const Model_SP& m);

				bool isLeaf() const override;
				void im_transform(void* dst, const AMat32& m) const override;
				void im_getBVolume(Circle& c) const override;
				void im_getBVolume(AABB& a) const override;
				float im_getInertia() const override;
				float im_getArea() const override;
				Vec2 im_getCenter() const override;
				Vec2 im_support(const Vec2& dir) const override;
				bool im_hitPoint(const Vec2& p, float t=NEAR_THRESHOLD) const override;
				Vec2 toLocal(const spn::Vec2& v) const override;
				Vec2 toLocalDir(const Vec2& v) const override;
				Vec2 toWorld(const Vec2& v) const override;
				Vec2 toWorldDir(const Vec2& v) const override;
				spn::Optional<const AMat32&> getToLocal() const override;
				spn::Optional<const AMat32&> getToWorld() const override;
				spn::Pose2D& tf_refPose() override;
				const spn::Pose2D& tf_getPose() const override;
				void setModelSource(const Model_SP& m);
				const Model_SP& getModelSource() const;
				std::ostream& print(std::ostream& os) const override;
				uint32_t getCID() const override;
				friend std::ostream& operator << (std::ostream&, const TfLeafBase&);
		};
		//! 座標変換ありのモデル基底
		template <class Ud=spn::none_t>
		class TfLeaf : public TfLeafBase {
			private:
				using base_t = TfLeafBase;
				Ud	_udata;
			public:
				using base_t::base_t;
				bool canCacheShape() const override {
					return false;
				}
				/*! ユーザーデータがvoidの時は親ノードのデータを返す */
				void* getUserData() override {
					return _getUserData(&_udata, std::is_same<spn::none_t, Ud>());
				}
				const void* getCore[[noreturn]]() const override {
					// canCacheShapeがfalseのオブジェクトに対して呼んではいけない
					AssertT(Trap, false, (std::domain_error)(const char*), "calling getCore() is not valid for TfLeaf")
					throw 0;
				}
				void* getCore[[noreturn]]() override {
					static_cast<const TfLeaf*>(this)->getCore();
				}
				TfBase_SP clone() const override {
					auto sp = std::make_shared<TfLeaf>(*this);
					// ModelSpのクローン
					sp->setModelSource(this->getModelSource()->im_clone());
					return std::move(sp);
				}
		};
		//! キャッシュ有りのTfLeaf
		template <class Shape, class Ud=spn::none_t>
		class TfLeafC : public TfLeaf<Ud> {
			private:
				using base_t = TfLeaf<Ud>;
				using ModelAccum = typename base_t::ModelAccum;
				using Global = typename base_t::Global;

				mutable uint32_t _shapeAccum = 0;
				mutable Shape	_tfShape;
				const Shape& getTfShape() const {
					auto accum = base_t::getShapeAccum();
					if(accum != _shapeAccum) {
						_shapeAccum = accum;
						this->getModelSource()->im_transform(&_tfShape, base_t::getGlobal());
					}
					return _tfShape;
				}
				Ud	_udata;
			public:
				using base_t::base_t;
				bool canCacheShape() const override {
					return true;
				}
				Vec2 im_support(const Vec2& dir) const override {
					// 予め変換しておいた形状でサポート写像
					return getTfShape().support(dir);
				}
				void im_getBVolume(Circle& c) const override {
					c = getTfShape().bs_getBCircle();
				}
				void im_getBVolume(AABB& a) const override {
					a = getTfShape().bs_getBBox();
				}
				const void* getCore() const override {
					getTfShape();
					return &_tfShape;
				}
				void* getCore() override {
					getTfShape();
					// データが改変されるかも知れないので更新フラグを立てる
					++base_t::refModelAccum();
					return &_tfShape;
				}
				TfBase_SP clone() const override {
					auto sp = std::make_shared<TfLeafC>(*this);
					// ModelSpのクローン
					sp->setModelSource(this->getModelSource()->im_clone());
					return std::move(sp);
				}
		};
	}
}
