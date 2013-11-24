//! 3D形状の定義
#pragma once
#include "geom.hpp"
#include "spinner/resmgr.hpp"
#include "spinner/size.hpp"
#include <memory>

namespace boom {
	namespace geo3d {
		using LNear = std::pair<Vec3,LINEPOS>;
		struct Sphere;
		struct Line;
		struct Ray;
		struct Segment;
		struct Capsule;
		struct Frustum;
		struct AABB;
		struct Convex;
		struct Cone;
		using CTGeo = spn::CType<Sphere, Line, Ray, Segment, Capsule, Frustum, AABB, Cone>;//, Convex>;
		template <class T>
		using ITagP = ITagP_base<T, CTGeo>;

		struct IModel;
		using MdlItr = PtrItr<IModel>;
		using MdlIP = std::pair<MdlItr, MdlItr>;
		struct Sphere : ITagP<Sphere> {
			Vec3 center;
			float radius;

			Sphere() = default;
			Sphere(const Vec3& c, float r): center(c), radius(r) {}
			// <キャッシュ対象のパラメータが同時に算出される場合のメモ>
			// std::tuple<TagCenter, TagArea> calcboth() const;
			// ---- cacheable functions ----
			const Vec3& bs_getGCenter() const;
			const Vec3& bs_getCenter() const;
			const Sphere& bs_getBSphere() const;
			float bs_getArea() const;
			Mat33 bs_getInertia() const;
			// -----------------------------
			Vec3 support(const Vec3& dir) const;
			Sphere operator * (const AMat43& m) const;

			spn::none_t hit(...) const;
			bool hit(const Line& l) const;
			bool hit(const Ray& r) const;
			bool hit(const Sphere& s) const;
			bool hit(const Capsule& c) const;

			static Sphere Cover(MdlItr mI, MdlItr mE);
		};
		struct IModel {
			IModel* pParent = nullptr;
			//! 子に変更があった事を親ノードに(あれば)伝える
			virtual void notifyChange();
			virtual void applyChange();
			//! 子ノードの取得
			virtual MdlIP getInner() const;
			//! モデルの実体
			virtual const void* getCore() const = 0;
			//! 最寄りのユーザーデータを取得
			/*! このノードが持っていればそれを返し、無ければ親を遡って探す */
			virtual void* getUserData() const;
			//! 形状の通し番号
			virtual uint32_t getCID() const = 0;

			friend std::ostream& operator << (std::ostream& os, const IModel& mdl);

			// ---- cacheable functions ----
			virtual Sphere im_getBSphere() const = 0;
			virtual Mat33 im_getInertia() const = 0;
			virtual float im_getArea() const = 0;
			virtual Vec3 im_getCenter() const = 0;
			virtual Vec3 im_getGCenter() const = 0;
			// -----------------------------
			virtual Vec3 im_support(const Vec3& dir) const = 0;
		};
		using UPModel = std::unique_ptr<IModel>;
		#define mgr_model3d (::boom::geo3d::ModelMgr::_ref())
		class ModelMgr : public spn::ResMgrA<UPModel, ModelMgr> {};
		DEF_HANDLE(ModelMgr, Mdl, UPModel)

		template <class T>
		struct Model : IModelP_base<T, IModel> {
			using IModelP_base<T, IModel>::IModelP_base;
			//TODO: 対応したルーチンを持っているか自動で判定して、無ければ例外を投げる仕組み
			// 各種ルーチンの中継
			Sphere im_getBSphere() const override { return T::bs_getBSphere(); }
			Mat33 im_getInertia() const override { return T::bs_getInertia(); }
			float im_getArea() const override { return T::bs_getArea(); }
			Vec3 im_getCenter() const override { return T::bs_getCenter(); }
			Vec3 im_getGCenter() const override { return T::bs_getGCenter(); }
			Vec3 im_support(const Vec3& dir) const override { return T::support(dir); }
		};
		//! 子ノードを含むIModelインタフェース
		template <class M, class UD=spn::none_t>
		struct ModelCh : IModelP_base<M, IModel> {
			using ChL = std::vector<M>;
			Sphere	_cover = Sphere(Vec3(0),0);
			ChL		_chL;
			bool	_bChange = false;
			UD		_udata;

			template <class MC>
			void addChild(MC&& mc) {
				_chL.push_back(std::forward<MC>(mc));
				_bChange = true;
			}
			M& operator [](int n) { return _chL[n]; }
			const M& operator [](int n) const { return _chL[n]; }

			MdlIP getInner() const override {
				if(_chL.empty())
					return MdlIP();
				int nC = _chL.size();
				return MdlIP(MdlItr(&_chL[0], sizeof(M)), MdlItr(&_chL[nC-1], sizeof(M)));
			}
			void applyChange() override {
				if(_bChange) {
					_bChange = false;
					for(auto& c : _chL)
						c.applyChange();
				}
			}
			virtual void* getUserData() {
				return _getUserData(std::is_same<spn::none_t, UD>());
			}
			void* _getUserData(std::true_type) {
				return (IModel::pParent) ? IModel::pParent->getUserData() : nullptr;
			}
			void* _getUserData(std::false_type) {
				return &_udata;
			}
		};
		// ------------------ Cache system ------------------
		DEF_CACHETAG(TagCenter, spn::Vec3, bs_getCenter)
		// DEF_CACHETAG(TagGCenter, spn::Vec3, bs_getGCenter)
		// DEF_CACHETAG(TagInertia, spn::Mat33, bs_getInertia)
		// DEF_CACHETAG(TagArea, float, bs_getArea)
		// DEF_CACHETAG(TagBSphere, Sphere, bs_getBSphere)
		// using CTTag = CCType<TagGCenter, TagCenter, TagInertia, TagArea, TagBSphere>;
		// template <class CORE>
		// class Cache : public ::boom::Cache<CTTag, CacheBase<CORE>>, public ITagP<CORE>, public IModel {
		// 	public:
		// 		using base = ::boom::Cache<CTTag, CacheBase<CORE>>;
		// 		DEF_GETMETHOD(base, getBSphere, TagBSphere)
		// 		DEF_GETMETHOD(base, getInertia, TagInertia)
		// 		DEF_GETMETHOD(base, getArea, TagArea)
		// 		DEF_GETMETHOD(base, getCenter, TagCenter)
		// 		DEF_GETMETHOD(base, getGCenter, TagGCenter)
		//
		// 		Sphere im_getBSphere() const override { return getBSphere(); }
		// 		Mat33 im_getInertia() const override { return getInertia(); }
		// 		float im_getArea() const override { return getArea(); }
		// 		Vec3 im_getCenter() const override { return getCenter(); }
		// 		Vec3 im_getGCenter() const { return getGCenter(); }
		// 		Vec3 support(const Vec3& dir) const { return getCoreRef().support(dir); }
		// 		const void* getCore() const override { return &base::getCoreRef(); }
		// 		uint32_t getCID() const override { return ITagP<CORE>::GetCID(); }
		// };
		using CTTag = CCType<TagCenter>;
		template <class CORE>
		class Cache : public ::boom::Cache<CTTag, CORE>, public ITagP<CORE>, public IModel {
			public:
				using base = ::boom::Cache<CTTag, CORE>;
				DEF_GETMETHOD(base, getCenter, TagCenter)
				Vec3 im_getCenter() const override { return getCenter(); }
				const void* getCore() const override { return &base::getCoreRef(); }
				uint32_t getCID() const override { return ITagP<CORE>::GetCID(); }
		};
		// <キャッシュ対象のパラメータが同時に算出される場合のメモ>
		// struct SphereC : CacheBase<Sphere> {
		//	DEF_SHARE(Sphere, calcboth, (TagCenter)(TagArea))
		// };

		#define DEF_INVALID_BSFUNCS \
			const Vec3& bs_getGCenter() const { INVOKE_ERROR } \
			const Vec3& bs_getCenter() const { INVOKE_ERROR } \
			const Sphere& bs_getBSphere() const { INVOKE_ERROR } \
			const float& bs_getArea() const { INVOKE_ERROR } \
			const Mat33& bs_getInertia() const { INVOKE_ERROR } \
			const Vec3& support(const Vec3& dir) const { INVOKE_ERROR }
		//! 点
		struct Point : Vec3, ITagP<Point> {
			constexpr static float NEAR_THRESHOLD = 1e-5f;

			using Vec3::Vec3;
			using Vec3::distance;
			// ---- cacheable functions ----
			const Vec3& bs_getGCenter() const;
			const Vec3& bs_getCenter() const;
			Sphere bs_getBSphere() const;
			const float& bs_getArea() const;
			const Mat33& bs_getInertia() const;
			// -----------------------------
			const Vec3& support(const Vec3& dir) const;

			spn::none_t hit(...) const;
			bool hit(const Point& p) const;
		};
		//! 直線
		struct Line : ITagP<Line> {
			Vec3	pos, dir;
			Line() = default;
			Line(const Vec3& p, const Vec3& d);
			static Line FromPoints(const Vec3& p0, const Vec3& p1);

			DEF_INVALID_BSFUNCS
			Line operator * (const AMat43& m) const;

			Vec3x2 nearestPoint(const Line& s) const;
			Vec3 nearest(const Vec3& p) const;
			float dist_sq(const Vec3& p) const;

			spn::none_t hit(...) const;
		};
		//! 半直線
		struct Ray : ITagP<Ray> {
			Vec3	pos, dir;
			Ray() = default;
			Ray(const Vec3& p, const Vec3& d);
			static Ray FromPoints(const Vec3& from, const Vec3& to);

			DEF_INVALID_BSFUNCS
			Ray operator * (const AMat43& m) const;

			Vec3 getPt(float len) const;
			Vec3x2 nearest(const Ray& r) const;
			Vec3 nearest(const Vec3& p) const;
			const Line& asLine() const;

			spn::none_t hit(...) const;
		};
		//! 線分
		struct Segment : ITagP<Segment> {
			Vec3	from, to;
			Segment() = default;
			Segment(const Vec3& p0, const Vec3& p1);

			// ---- cacheable functions ----
			Vec3 bs_getGCenter() const;
			Vec3 bs_getCenter() const;
			Sphere bs_getBSphere() const;
			const float& bs_getArea() const { INVOKE_ERROR }
			// -----------------------------
			const Vec3& support(const Vec3& dir) const;
			Segment operator * (const AMat43& m) const;

			float dist_sq(const Vec3& p) const;
			float getRatio(const Vec3& pos) const;
			LNear nearest(const Vec3& l) const;
			Vec3x2 nearest(const Segment& l) const;
			float getLength() const;
			Vec3 getDir() const;
			Line asLine() const;
			Ray asRay() const;

			/*! \return 交差したかの符号判定, from側距離の絶対値, to側距離の絶対値 */
			std::tuple<bool, float,float> _crossPoint(const Plane& plane) const;
			//! 平面との交差判定及び交差点の算出
			/*! 線分が面を横切っていない場合は接触なしと判定される */
			spn::Optional<Vec3> crossPoint(const Plane& plane) const;
			/*! 片方が面に接触している場合はその座標、両方なら中間点を出力 */
			spn::Optional<Vec3> crossPointFit(const Plane& plane, float threshold) const;

			spn::none_t hit(...) const;
			bool hit(const Sphere& s) const;
			bool hit(const Capsule& c) const;
			bool hit(const Plane& p) const;
		};
		//! カプセル
		struct Capsule : ITagP<Capsule> {
			Segment seg;
			float	radius;

			// ---- cacheable functions ----
			Vec3 bs_getCenter() const;
			Vec3 bs_getGCenter() const;
			Sphere bs_getBSphere() const;
			float bs_getArea() const;
			Mat33 bs_getInertia() const;
			// -----------------------------
			Vec3 support(const Vec3& dir) const;
			Capsule operator * (const AMat43& m) const;

			spn::none_t hit(...) const;
			bool hit(const Point& p) const { return false; }
			bool hit(const Sphere& s) const { return false; }
			bool hit(const Capsule& c) const { return false; }
		};
		using CapsuleM = Model<Capsule>;
		using CapsuleC = Cache<CacheBase<Capsule>>;

		using SphereM = Model<Sphere>;

		//! 視点がZ方向を向いている時の視錐台パラメータ
		struct VFrusFloat {
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
		//! 四角錐
		/*! 距離1,幅1,高さ1の四角錐を基本としFovの調整はスケーリングで対応 */
		struct Frustum : ITagP<Frustum>, public Pose3D {
			struct Points {
				enum POS { LT,RT,LB,RB, NUM_POS };
				union {
					struct {
						Vec3 center,
							pt[NUM_POS];
					};
					Vec3 ar[NUM_POS+1];
				};
				Points() = default;
				Points(const Points& pts);
				Points& operator = (const Points& pts);
				bool hit(const Plane& p) const;
			};

			Frustum() = default;
			using Pose3D::Pose3D;
			Frustum(const Vec3& ori, const Vec3& dir, const Vec3& up, float fov, float dist, float aspect);
			Frustum(const Vec3& ori, const Vec3& at, const Vec3& up, float fov, float aspect);

			// ---- cacheable functions ----
			Vec3 bs_getGCenter() const;
			Vec3 bs_getCenter() const;
			Sphere bs_getBSphere() const;
			float bs_getArea() const;
			Mat33 bs_getInertia() const;
			// -----------------------------
			Vec3 support(const Vec3& dir) const;
			Frustum operator * (const AMat43& m) const;

			spn::none_t hit(...) const;
			bool hit(const Sphere& s) const;
			bool hit(const Frustum& f) const;
			bool hit(const Cone& cone) const;
			bool hit(const Plane& plane) const;

			//! 5頂点と引数の平面との距離をコールバックで順次受け取る
			void iterateLocalPlane(const Plane& plane, std::function<void (float)> cb) const;
			//! 上、右の平面法線を計算(但しローカル)
			Vec3x2 getUpRightLocalDir() const;
			//! 辺がAを通過しているか
			bool crossHit_LocalB(const Points& pts) const;
			//! 点がAに入っているか
			bool inHit_LocalB(const Points& pts) const;
			//! 面との相対位置
			int checkSide(const Plane& plane, float threshold) const;
			//! 四錐台の4 + 1頂点を計算
			/*! \param[in] bWorld ワールド座標系出力フラグ
				\param[in] */
			Points getPoints(bool bWorld) const;
			//! 指定された行列で変換した時の透視平面上での大きさ
			spn::Rect get2DRect(const AMat43& mV, const AMat44& mP) const;
		};
		//! 円錐
		struct Cone : ITagP<Cone> {
			Vec3	center, dir;
			float	radius, length;

			// ---- cacheable functions ----
			Vec3 bs_getGCenter() const;
			Vec3 bs_getCenter() const;
			Sphere bs_getBSphere() const;
			float bs_getArea() const;
			Mat33 bs_getInertia() const;
			// -----------------------------
			float getAngle() const;
			//! 平面に最も近い座標，または深度を計算
			/*! \return 最深点, 深度距離 */
			std::pair<Vec3,float> nearestPoint(const Plane& plane) const;

			spn::none_t hit(...) const;
// 			bool hit(const Vec3& p) const;
// 			bool hit(const Segment& s) const;
// 			bool hit(const Plane& p) const;

// 			Vec3 support(const Vec3& dir) const;
// 			Cone& operator * (const AMat43& m) const;
		};
		//! AxisAligned Box
		struct AABB : ITagP<AABB> {
			Vec3	vmin, vmax;
			AABB() = default;
			AABB(const Vec3& v_min, const Vec3& v_max);
			static AABB FromPoints(const Vec3* v, int n);

			// ---- cacheable functions ----
			Vec3 bs_getGCenter() const;
			Vec3 bs_getCenter() const;
			Sphere bs_getBSphere() const;
			float bs_getArea() const;
			Mat33 bs_getInertia() const;
			// -----------------------------
			Vec3 support(const Vec3& dir) const;
			AABB operator * (const AMat43& m) const;

			spn::none_t hit(...) const;
			bool hit(const AABB& ab) const;
			bool hit(const Plane& plane) const;
		};

		//! Narrow Phase判定
		using Narrow = ::boom::Narrow<CTGeo, IModel>;

		//! 点から一番近い直線(線分)上の点
		template <class CLIP>
		inline Vec3 NearestPoint(const Line& ls, const Vec3& p, CLIP clip) {
			Vec3 toP = p - ls.pos;
			float d = ls.dir.dot(toP);
			return ls.pos + ls.dir * clip(d);
		}
		//! 直線(線分)の最近傍点対
		template <class CLIP0, class CLIP1>
		inline Vec3x2 NearestPoint(const Line& ls0, const Line& ls1, CLIP0 clip0, CLIP1 clip1) {
			return Vec3x2();
		}
	}
}
