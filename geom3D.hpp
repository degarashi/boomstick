//! 3D形状の定義
#pragma once
#include "geom.hpp"
#include "spinner/resmgr.hpp"
#include <memory>

namespace boom {
	namespace geo3d {
		using LNear = std::tuple<Vec3,LINEPOS>;
		struct Point;
		struct Sphere;
		struct Capsule;
		struct Line;
		struct AABB;
		struct OABB;
		struct Cone;
		struct Poly;
		struct Convex;
		struct Frustum;
// 		using CTGeo = spn::CType<Point, Sphere, Capsule, Line, AABB, OABB, Cone, Poly, Convex, Frustum>;
		using CTGeo = spn::CType<Point, Line, Sphere, Capsule, Frustum>;
		template <class T>
		using ITagP = ITagP_base<T, CTGeo>;

		struct Sphere : ITagP<Sphere> {
			Vec3 center;
			float radius;

// 			std::tuple<TagGCenter, TagArea> calcboth() const;
			const Vec3& bs_getGCenter() const;
			const Vec3& bs_getCenter() const;
			const Sphere& bs_getBSphere() const;
			float bs_getArea() const;
			Mat33 bs_getInertia() const;
			Vec3 support(const Vec3& dir) const;

			spn::none_t hit(...) const;
			bool hit(const Sphere& s) const { return false; }
			bool hit(const Capsule& c) const { return false; }
			Sphere() = default;
			Sphere(const Vec3& c, float r);
		};

		struct IModel {
			virtual std::tuple<IModel*, IModel*> getInner() const { return std::make_pair(nullptr, nullptr); }
			virtual const void* getCore() const { INVOKE_ERROR }
			virtual uint32_t getCID() const = 0;

			virtual Sphere im_getBSphere() const { INVOKE_ERROR }
			virtual Mat33 im_getInertia() const { INVOKE_ERROR }
			virtual float im_getArea() const { INVOKE_ERROR }
			virtual Vec3 im_getCenter() const { INVOKE_ERROR }
			virtual Vec3 im_getGCenter() const { INVOKE_ERROR }
			virtual Vec3 support() const { INVOKE_ERROR }

			friend std::ostream& operator << (std::ostream& os, const IModel& mdl);

			virtual spn::Optional<Mat43> getTransform() const { return spn::none; }
			virtual Vec3 toLocal(const Vec3& v) const { return v; }
			virtual Vec3 toLocalDir(const Vec3& v) const { return v; }
			virtual Vec3 toWorld(const Vec3& v) const { return v; }
			virtual Vec3 toWorldDir(const Vec3& v) const { return v; }
			static AMat43 s_identityMat;
			virtual const AMat43& getToLocal() const { return s_identityMat; }
			virtual const AMat43& getToWorld() const { return s_identityMat; }
		};
		using UPModel = std::unique_ptr<IModel>;
		#define mgr_model3d (::boom::geo3d::ModelMgr::_ref())
		class ModelMgr : public spn::ResMgrA<UPModel, ModelMgr> {};
		DEF_HANDLE(ModelMgr, Mdl, UPModel)

		template <class T>
		struct Model : IModelP_base<T, IModel> {
			// 各種ルーチンの中継
			Sphere im_getBSphere() const override { return T::bs_getBSphere(); }
			Mat33 im_getInertia() const override { return T::bs_getInertia(); }
			float im_getArea() const override { return T::bs_getArea(); }
			Vec3 im_getCenter() const override { return T::bs_getCenter(); }
			Vec3 im_getGCenter() const override { return T::bs_getGCenter(); }
			Vec3 support(const Vec3& dir) const override { return T::support(dir); }
		};

		// ------------------ Cache system ------------------
		DEF_CACHETAG(TagGCenter, spn::Vec3, bs_getGCenter)
		DEF_CACHETAG(TagCenter, spn::Vec3, bs_getCenter)
		DEF_CACHETAG(TagInertia, spn::Mat33, bs_getInertia)
		DEF_CACHETAG(TagArea, float, bs_getArea)
		DEF_CACHETAG(TagBSphere, Sphere, bs_getBSphere)

		using CTTag = CCType<TagGCenter, TagCenter, TagInertia, TagArea, TagBSphere>;
		template <class CORE>
		class Cache : public ::boom::Cache<CTTag, CacheBase<CORE>>, public ITagP<CORE>, public IModel {
			public:
				using base = ::boom::Cache<CTTag, CacheBase<CORE>>;
				DEF_GETMETHOD(base, getBSphere, TagBSphere)
				DEF_GETMETHOD(base, getInertia, TagInertia)
				DEF_GETMETHOD(base, getArea, TagArea)
				DEF_GETMETHOD(base, getCenter, TagCenter)
				DEF_GETMETHOD(base, getGCenter, TagGCenter)

				const void* getCore() const override { return &base::getCoreRef(); }
				uint32_t getCID() const override { return ITagP<CORE>::GetCID(); }
		};
		// --------------------------------------------------
// 		struct SphereC0 : CacheBase<Sphere> {
// 			DEF_SHARE(Sphere, calcboth, (TagGCenter)(TagArea))
// 		};

		struct Point : Vec3, ITagP<Point> {
			constexpr static float NEAR_THRESHOLD = 1e-5f;

			using Vec3::Vec3;
			using Vec3::distance;
			float distance(const Line& l) const;
			LNear nearest(const Line& l) const;

			spn::none_t hit(...) const;
			bool hit(const Point& p) const { return false; }
		};
		using PointM = Model<Point>;
		using PointC = Cache<Point>;

		struct Line : ITagP<Line> {
			Vec3	from, to;

			Vec3 bs_getGCenter() const;
			Vec3 bs_getCenter() const;
			Sphere bs_getBSphere() const;
			float bs_getArea() const;
			Mat33 bs_getInertia() const;
			Vec3 suppoer(const Vec3& dir) const;

			spn::none_t hit(...) const;
			bool hit(const Sphere& s) const { return false; }
			bool hit(const Capsule& c) const { return false; }
		};
		using LineM = Model<Line>;
		using LineC = Cache<Line>;

		using SphereM = Model<Sphere>;
		class SphereC : public Cache<Sphere> {
			public:
				void setCenter(const Vec3& c) {
					getCoreRef().center = c;
					_invalidate<TagArea>();
				}
		};

		struct Capsule : ITagP<Capsule> {
			Vec3	from, to;
			float	radius;

			float bs_getArea() const;
			Mat33 bs_getInertia() const;
			Vec3 bs_getCenter() const;
			Vec3 bs_getGCenter() const;
			Sphere bs_getBSPhere() const;
			Vec3 support(const Vec3& dir) const;

			spn::none_t hit(...) const;
			bool hit(const Point& p) const { return false; }
			bool hit(const Sphere& s) const { return false; }
			bool hit(const Capsule& c) const { return false; }
			Capsule() = default;
			Capsule(const Capsule&) = default;
			Capsule(const Vec3& beg, const Vec3& end, float r);
		};
		using CapsuleM = Model<Capsule>;
		using CapsuleC = Cache<Capsule>;

		struct Frustum : ITagP<Frustum> {
			Vec3 from, to;
			float fovV, aspect;

			Vec3 bs_getGCenter() const;
			Vec3 bs_getCenter() const;
			Sphere bs_getBSphere() const;
			float bs_getArea() const;
			Mat33 bs_getInertia() const;
			Vec3 support(const Vec3& dir) const;

			spn::none_t hit(...) const;
		};
		using FrustumM = Model<Frustum>;
		using FrustumC = Cache<Frustum>;

		using Narrow = ::boom::Narrow<CTGeo, IModel>;
	}
}
