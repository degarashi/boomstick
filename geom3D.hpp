#pragma once
#include "spinner/matrix.hpp"
#include "spinner/type.hpp"
#include "spinner/misc.hpp"
#include "spinner/pose.hpp"
#include "spinner/withbool.hpp"
#include "spinner/resmgr.hpp"
#include "collision.hpp"

namespace boom {
	template<class CORE >
	class Cache;
	using spn::Vec3;
	using spn::AVec3;
	using spn::Mat33;
	using spn::Mat43;
	using spn::Mat44;
	using apn::AMat43;
	using spn::Plane;

	namespace geo3d {
		struct Point;
		struct Capsule;
		struct Line;
		struct AABB;
		struct OABB;
		struct Cone;
		struct Poly;
		struct Convex;
		struct Frustum;
		using CTGeo = CType<Point, Capsule, Line, AABB, OABB, Cone, Poly, Convex, Frustum>;
		template <class T>
		using ITagP = ITagPB<T, CTGeo>;
		template <class T>
		using IModelP = IModelPB<T, CTGeo>;

		using PointL = std::vector<Vec3>;
		using LineL = std::vector<Line>;
		using LNear = std::pair<Vec3, LINEPOS>;

		struct IModel {
			static AMat43 s_identityMat;
			// 2DのIModelと殆ど同じなのでコメントなどは2Dの方を参照

			virtual Vec3 support(const Vec3& dir) const = 0;
			virtual Vec3 getCenter() const = 0;
			virtual uint32_t getCID() const = 0;
			virtual bool isInner(const Vec3& /*pos*/) const { return false; }
			virtual Sphere getBSphere() const;
			virtual float getArea(bool /*bInv*/) const { INVOKE_ERROR }
			virtual Mat33 getInertia(bool /*bInv*/) const { INVOKE_ERROR }
			virtual bool hit(const IModel& m) const { INVOKE_ERROR }

			friend std::ostream& operator << (std::ostream& os, const IModel& mdl);

			virtual Vec3 toLocal(const Vec3& v) const { return v; }
			virtual Vec3 toLocalDir(const Vec3& v) const { return v; }
			virtual Vec3 toWorld(const Vec3& v) const { return v; }
			virtual Vec3 toWorldDir(const Vec3& v) const { return v; }
			virtual const AMat43& getToLocal() const { return s_identityMat; }
			virtual const AMat43& getToWorld() const { return s_identityMat; }
		};
		using UPModel = std::unique_ptr<IModel>;
		#define mgr_model3d	boom::geo3d::ModelMgr::_ref()
		class ModelMgr : public spn::ResMgrA<UPModel, ModelMgr> {};
		DEF_HANDLE(ModelMgr, Mdl, UPModel)

		#define DEF_IMODEL_FUNCS3D \
			Vec3 support(const Vec3& dir) const override; \
			Vec3 getCenter() const override;
		#define DEF_IMODEL_MASS3D \
			float getArea(bool bInv=false) const override; \
			Mat33 getInertia(bool bInv=false) const override; \
			Sphere getBSphere() const override;

		struct Point : Vec3, ITagP<Point> {
			constexpr static float NEAR_THRESHOLD = 1e-5f;

			using Vec3::Vec3;
			using Vec3::distance;
			float distance(const Line& l) const;
			LNear nearest(const Line& l) const;
			bool hit(const Point& p) const;
		};
		struct PointM : IModelP<Point> {
			using Point::Point;
			DEF_IMODEL_FUNCS3D
		};
		struct Sphere : ITagP<Sphere> {
			using CT = CType<TagArea, TagInertia>;
			template <class T>
			using Wrap = CoreRaw<T>;

			float area() const;
			Mat33 inertia() const;
			const Vec3& center() const;
			const Sphere& bsphere() const;
			// --------------------------------------------
			Vec3	vCenter;
			float	fRadius;

			Sphere() = default;
			Sphere(const Vec3& c, float r);

			Vec3 support(const Vec3& dir) const;
			bool hit(const Vec3& pt) const;
			bool hit(const Sphere& s) const;
		};
		struct SphereM : IModelP<Sphere> {
			using Sphere::Sphere;
			DEF_IMODEL_FUNCS3D
		};
		class SphereC : public IModelP<Sphere>, public Cache<Sphere> {
			Sphere	_sp;
			public:
				SphereC() = default;
				SphereC(const Sphere& c);
				DEF_IMODEL_MASS3D
		};
		template <class BASE>
		class TModel : public BASE {
			public:
				using BASE::BASE;
		};
		using SphereCT = TModel<SphereC>;

		// TODO: 適時必要になり次第実装する
		struct Line {};
		struct Ray {};
		struct StLine {};
		struct Capsule {};
		struct Poly {};
		struct Convex {};
		struct Cone {};
		struct AABB {};
		struct OABB {};

		struct FrusFloat {
			// Near and far plane distances
			float		nearDistance;
			float		farDistance;

			// Precalculated normal components
			float		leftRightX;
			float		leftRightZ;
			float		topBottomY;
			float		topBottomZ;

			void init(float nearZ, float farZ, float aspect, float fov);
		};
		//! (主に)可視カリング用視錐台構造体
		struct Frustum : ITagP<Frustum> {
			struct Points {
				enum POS {LT, RT, LB, RB, NUM_POS};
				Vec3 center,
					pt[NUM_POS];
				bool hit(const spn::Plane& p) const;
			};

			static Frustum FromView(const Vec3& ori, const Vec3& dir, const Vec3& up, float fov, float ndist, float dist, float aspect, float expandF);
			static Frustum FromLookAt(const Vec3& ori, const Vec3& at, const Vec3& up, float fov, float aspect, float expandF);
			static Frustum FromPtDir(const Vec3& origin, const Vec3& dir, const Vec3& up, float w, float h, float length);

			std::tuple<Vec3,Vec3,Vec3> getAxis() const;
			Points getPoints() const;

			bool hit(const Vec3& p) const;
			bool hit(const Plane& p) const;
			bool hit(const Cone& cone) const;
			bool hit(const Frustum& fr) const;
			bool hit(const Sphere& sp) const;
			Vec3 getCenter() const;
		};
		struct FrustumM : IModelP<Frustum> {
			using Frustum::Frustum;
		};
		class FrustumC : public IModelP<Frustum>, public Cache<Frustum> {
			Frustum		_frus;
			public:
				FrustumC() = default;
				FrustumC(const FrustumC& f);
				Rect get2DRect(const Mat43& mV, const Mat44& mP) const;
		};
	}
}
