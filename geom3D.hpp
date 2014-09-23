//! 3D形状の定義
#pragma once
#include "geom.hpp"
#include "spinner/resmgr.hpp"
#include "spinner/size.hpp"
#include "spinner/layerbit.hpp"
#include "spinner/structure/freeobj.hpp"
#include <memory>
#include <unordered_set>
#include <map>

namespace boom {
	namespace geo2d {
		struct Convex;
	}
	namespace geo3d {
		struct ConvexUD_Col;
		template <class, class>
		class Convex;
		using ColCv = Convex<Vec3, ConvexUD_Col>;

		constexpr float DOT_TOLERANCE = 1e-4f,
						CREATE_PLANE_DIST = 5e2f,
						PORTALBOX_VOLUME = 5e2f,
						VTX_NEAR_TOLERANCE = 1e-5f,
						NEARPLANE_TOLERANCE = 1e-3f,
						LINEAR_TOLERANCE = 1.0f - 1e-6f;
		using LNear = std::pair<Vec3,LinePos>;
		struct Point;
		struct Sphere;
		struct Line;
		struct Ray;
		struct Segment;
		struct Capsule;
		struct Frustum;
		struct AABB;
		struct Cone;
		struct Cylinder;
		class ConvexP;
		using CTGeo = spn::CType<ConvexP, AABB, Frustum, Cone, Cylinder, Capsule, Sphere, Segment, Ray, Line, Point>;
		template <class T>
		using ITagP = ITagP_base<T, CTGeo>;

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
			const Sphere& bs_getBVolume() const;
			float bs_getArea() const;
			Mat33 bs_getInertia() const;
			// -----------------------------
			Vec3 support(const Vec3& dir) const;
			Sphere operator * (const AMat43& m) const;

			spn::none_t hit(...) const;
			bool hit(const Ray& r) const;
			bool hit(const Sphere& s) const;
			bool hit(const Capsule& c) const;
			bool hit(const Line& ls) const;
			bool hit(const Segment& s) const;

			static Sphere Cover(PtrItr mI, PtrItr mE);
		};
		struct IModel : ::boom::IModelNode {
			//! 形状の通し番号
			virtual uint32_t getCID() const = 0;
			// ---- cacheable functions ----
			virtual Sphere im_getBVolume() const = 0;
			virtual Mat33 im_getInertia() const = 0;
			virtual float im_getArea() const = 0;
			virtual Vec3 im_getCenter() const = 0;
			virtual Vec3 im_getGCenter() const = 0;
			// -----------------------------
			virtual Vec3 im_support(const Vec3& dir) const = 0;

			virtual Vec3 toLocal(const Vec3& v) const { return v; }
			virtual Vec3 toLocalDir(const Vec3& v) const { return v; }
			virtual Vec3 toWorld(const Vec3& v) const { return v; }
			virtual Vec3 toWorldDir(const Vec3& v) const { return v; }
			virtual spn::Optional<const AMat43&> getToLocal() const { return spn::none; }
			virtual spn::Optional<const AMat43&> getToWorld() const { return spn::none; }
			const static AMat43 cs_idMat;
			const AMat43& getToLocalI() const;
			const AMat43& getToWorldI() const;
		};
		using UPModel = std::unique_ptr<IModel>;
		#define mgr_model3d (::boom::geo3d::ModelMgr::_ref())
		class ModelMgr : public spn::ResMgrA<UPModel, ModelMgr> {};
		DEF_AHANDLE(ModelMgr, Mdl, UPModel, UPModel)

		template <class T>
		struct Model : IModelP_base<T, IModel> {
			using IModelP_base<T, IModel>::IModelP_base;
			//TODO: 対応したルーチンを持っているか自動で判定して、無ければ例外を投げる仕組み
			// 各種ルーチンの中継
			Sphere im_getBVolume() const override { return T::bs_getBVolume(); }
			Mat33 im_getInertia() const override { return T::bs_getInertia(); }
			float im_getArea() const override { return T::bs_getArea(); }
			Vec3 im_getCenter() const override { return T::bs_getCenter(); }
			Vec3 im_getGCenter() const override { return T::bs_getGCenter(); }
			Vec3 im_support(const Vec3& dir) const override { return T::support(dir); }
		};
		//! IModelに行列変換をプラス
		template <class T>
		class TModel : public ITagP<T>, public IModel {
			using VMdl = VModel<ModelMgr>;
			VMdl		_model;
			AMat43 		_mToWorld,
						_mToLocal;
			void _calcInv() {
				_mToWorld.convertA44().inversion(reinterpret_cast<AMat44&>(_mToLocal));
			}
			public:
				TModel(const TModel&) = default;
				TModel(const IModel& mdl, const AMat43& mw): _model(mdl), _mToWorld(mw) { _calcInv(); }
				TModel(HMdl hMdl, const AMat43& mw): _model(hMdl), _mToWorld(mw) { _calcInv(); }

				Sphere im_getBVolume() const override {
					return _model.get().im_getBVolume() * _mToWorld; }
				Mat33 im_getInertia() const override {
					return _model.get().im_getInertia().convertA44() * _mToWorld; }
				float im_getArea() const override {
					return _model.get().im_getArea(); }
				Vec3 im_getCenter() const override {
					return toWorld(_model.get().im_getCenter()); }
				Vec3 im_getGCenter() const override {
					return toWorld(_model.get().im_getGCenter()); }
				Vec3 im_support(const Vec3& dir) const override {
					return toWorldDir(_model.get().im_support(toLocalDir(dir))); }

				Vec3 toLocal(const Vec3& v) const override {
					return v.asVec4(1) * _mToLocal; }
				Vec3 toLocalDir(const Vec3& v) const override {
					return v.asVec4(0) * _mToLocal; }
				Vec3 toWorld(const Vec3& v) const override {
					return v.asVec4(1) * _mToWorld; }
				Vec3 toWorldDir(const Vec3& v) const override {
					return v.asVec4(1) * _mToWorld; }
				spn::Optional<const AMat43&> getToLocal() const override { return _mToLocal; }
				spn::Optional<const AMat43&> getToWorld() const override { return _mToWorld; }
				const void* getCore() const override { return _model.get().getCore(); }
				uint32_t getCID() const override { return ITagP<T>::GetCID(); }
		};
		// ------------------ Cache system ------------------
		DEF_CACHETAG(TagCenter, spn::Vec3, bs_getCenter)
		// DEF_CACHETAG(TagGCenter, spn::Vec3, bs_getGCenter)
		// DEF_CACHETAG(TagInertia, spn::Mat33, bs_getInertia)
		// DEF_CACHETAG(TagArea, float, bs_getArea)
		// DEF_CACHETAG(TagBSphere, Sphere, bs_getBVolume)
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
		// 		Sphere im_getBVolume() const override { return getBSphere(); }
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
			const Sphere& bs_getBVolume() const { INVOKE_ERROR } \
			const float& bs_getArea() const { INVOKE_ERROR } \
			const Mat33& bs_getInertia() const { INVOKE_ERROR } \
			const Vec3& support(const Vec3& /*dir*/) const { INVOKE_ERROR }
		//! 点
		struct Point : Vec3, ITagP<Point> {
			constexpr static float NEAR_THRESHOLD = 1e-5f;

			using Vec3::Vec3;
			using Vec3::distance;
			// ---- cacheable functions ----
			const Vec3& bs_getGCenter() const;
			const Vec3& bs_getCenter() const;
			Sphere bs_getBVolume() const;
			const float& bs_getArea() const;
			const Mat33& bs_getInertia() const;
			// -----------------------------
			const Vec3& support(const Vec3& dir) const;

			spn::none_t hit(...) const;
			bool hit(const Vec3& p) const;
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
			bool hit(const Vec3& p) const;
			bool hit(const Line& ls) const;
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
			Sphere bs_getBVolume() const;
			Mat33 bs_getInertia() const { INVOKE_ERROR }
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
			bool hit(const Plane& p) const;
		};
		using SegmentM = Model<Segment>;

		//! カプセル
		struct Capsule : ITagP<Capsule>, Segment {
			float	radius;

			// ---- cacheable functions ----
			Vec3 bs_getCenter() const;
			Vec3 bs_getGCenter() const;
			Sphere bs_getBVolume() const;
			float bs_getArea() const;
			Mat33 bs_getInertia() const;
			// -----------------------------
			Vec3 support(const Vec3& dir) const;
			Capsule operator * (const AMat43& m) const;

			spn::none_t hit(...) const;
			bool hit(const Vec3& p) const;
			bool hit(const Segment& s) const;
			bool hit(const Capsule& c) const;

			using ITagP<Capsule>::GetCID;
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

			void init(float nearZ, float farZ, float aspect, spn::RadF fov);
		};
		//! 四角錐
		/*! 距離1,幅1,高さ1の四角錐を基本としFovの調整はスケーリングで対応 */
		struct Frustum : ITagP<Frustum>, public Pose3D {
			struct Points {
				enum Pos {
					Center,
					LeftTop,
					RightTop,
					LeftBottom,
					RightBottom,
					NumPos
				};
				Vec3 point[NumPos];

				Points() = default;
				Points(const Vec3& cen, const Vec3& lt, const Vec3& rt, const Vec3& lb, const Vec3& rb);
				Points(const Points& pts);
				bool hit(const Plane& p) const;
				Points& operator = (const Points& pts);
				Points operator * (const AMat43& m) const;
			};
			struct Planes {
				enum Pos {
					Left,
					Right,
					Bottom,
					Top,
					Front,
					NumPlane
				};
				Plane plane[NumPlane];

				Planes() = default;
				Planes(const Plane& pL, const Plane& pR, const Plane& pB, const Plane& pT, const Plane& pF);
				Planes(const Planes& ps);
				Planes& operator = (const Planes& ps);
				Planes operator * (const AMat43& m) const;
			};
			const static Points cs_pointsLocal;
			const static Planes cs_planesLocal;

			Frustum() = default;
			using Pose3D::Pose3D;
			Frustum(const Vec3& ori, const Vec3& dir, const Vec3& up, spn::RadF fov, float dist, float aspect);
			Frustum(const Vec3& ori, const Vec3& at, const Vec3& up, spn::RadF fov, float aspect);

			// ---- cacheable functions ----
			Vec3 bs_getGCenter() const;
			Vec3 bs_getCenter() const;
			Sphere bs_getBVolume() const;
			float bs_getArea() const;
			Mat33 bs_getInertia() const;
			// -----------------------------
			Vec3 support(const Vec3& dir) const;
			Frustum operator * (const AMat43& m) const;

			spn::none_t hit(...) const;
			bool hit(const Sphere& s) const;
			bool hit(const Cone& c) const;
			bool hit(const Frustum& f) const;
			bool hit(const Plane& p) const;

			//! 5頂点と引数の平面との距離をコールバックで順次受け取る
			/*! \param[in] cb [(0=中央,1以降はPoints::POSの順番と同じ), 面との距離] */
			void iterateLocalPlane(const Plane& plane, std::function<void (int, float)> cb) const;
			//! 上、右の平面法線を計算(但しローカル)
			Vec3x2 getUpRightLocalDir() const;
			//! 辺がAを通過しているか
			bool crossHit_LocalB(const Points& pts) const;
			//! 点がAに入っているか
			bool inHit_LocalB(const Points& pts) const;
			//! 面との相対位置
			int checkSide(const Plane& plane, float threshold) const;
			//! 四錐台の4 + 1頂点を計算
			/*! \param[in] bWorld ワールド座標系出力フラグ */
			Points getPoints(bool bWorld) const;
			Planes getPlanes(bool bWorld) const;
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
			Sphere bs_getBVolume() const;
			float bs_getArea() const;
			Mat33 bs_getInertia() const;
			// -----------------------------
			float getAngle() const;
			//! 平面に最も近い座標，または深度を計算
			/*! \return 最深点, 深度距離 */
			std::pair<Vec3,float> nearestPoint(const Plane& plane) const;

			spn::none_t hit(...) const;
			bool hit(const Vec3& p) const;
			bool hit(const Plane& p) const;
			bool hit(const Segment& s) const;

			Vec3 support(const Vec3& dir) const;
			Cone operator * (const AMat43& m) const;
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
			Sphere bs_getBVolume() const;
			float bs_getArea() const;
			Mat33 bs_getInertia() const;
			// -----------------------------
			Vec3 support(const Vec3& dir) const;
			AABB operator * (const AMat43& m) const;

			spn::none_t hit(...) const;
			bool hit(const AABB& ab) const;
			bool hit(const Plane& p) const;
		};
		//! 各座標を格納するまでもないけどポリゴンを管理したい時のクラス
		/*! CnvPとセットで使用する */
		struct Idx3 {
			Idx3*	_neighbor[3];	//!< 隣接するポリゴン. エッジ順
			Plane	_plane;
			uint8_t	_id[3],			//!< 構成する頂点のインデックス
					_pid;			//!< ポリゴンインデックス
			public:
				Idx3();
				//! インデックスだけ初期化
				Idx3(int i0, int i1, int i2, int pid);
				//! 初期化と同時にポリゴン平面を計算
				Idx3(int i0, int i1, int i2, const Vec3List& vtx, int pid);
				void flip();
				void iterateEdge(std::function<void (int,int,Idx3*)> cb);
				const Plane& getPlane() const;
				int getID(int n) const;
				std::pair<int,int> getEdge(int n) const;
				void rewriteNeighbor(int idx, Idx3* p, bool bOther=false);
				void initNeighbor(Idx3* p0, Idx3* p1, Idx3* p2);
				Idx3* getNeighbor(int n);
				int findEdge(int id) const;
		};
		using Idx3List = std::vector<Idx3>;
		using IdxP3Set = std::unordered_set<Idx3*>;
		//! quickhullをする際の補助クラス
		class CnvP {
			class PolyF : public spn::FreeObj<Idx3,false> {
				public:
					PolyF(int n): FreeObj(n) {}
					template <class... Args>
					Ptr get(Args&&... args) {
						return FreeObj::get(std::forward<Args>(args)..., getNextID());
					}
			};
			//! ポリゴン張り直しをする場合のリストクラス
			struct PolyEdge {
				int 	_srcIdx[2];
				Idx3* 	_dstPoly;			//!< 接続先ポリゴン (_srcIdx[1]から始まるエッジ)
				PolyEdge(Idx3* poly, int idx0, int idx1) {
					_dstPoly = poly;
					_srcIdx[0] = idx0;
					_srcIdx[1] = idx1;
				}
			};
			using PEVec = std::vector<PolyEdge>;
			using PtrL = std::array<PolyF::Ptr, 0x100>;
			using VCand = std::array<uint8_t, 0x100>;
			using UseVID = std::array<uint8_t, 0x100>;

			PolyF		_poly;			//!< ポリゴンの実体
			PtrL		_pcand,			//!< ポリゴン候補リスト
						_pface;			//!< ポリゴン確定リスト
			VCand		_vcand;			//!< 頂点候補リスト
			UseVID		_useVID;		//!< 現在凸殻に使われている頂点の参照回数)
			int			_pcCur,
						_pfCur,
						_vcCur;
			const Vec3List& _vtx;
			//! 頂点カウントのイテレーションを高速に行うためのビット配列 (レイヤー)
			spn::LayerBitArray<0x100, uint8_t>		_bfVID;

			void _addRefVID(int vid);
			void _unRefVID(int vid);
			PolyF::Ptr& _addPoly(int i0, int i1, int i2);	//!< ポリゴン生成 & 頂点リストへの追加など
			void _remPoly(Idx3* p);					//!< ポリゴンのメモリ解放 & 頂点リストからの解除など
			Idx3* _pickPoly();						//!< 候補ポリゴンから1つポリゴンを抜き出す (本当は何か基準が要る)
			void _remTopPoly();						//!< 前回選ばれた候補を削除
			/*! @param[out] dstI 消去されるポリゴンリスト
				@param[out] dst 新たに面を張るエッジリスト
				@param[in] vert 基準頂点
				@param[in] poly これから処理するポリゴン */
			void _enumEdge(PEVec& dst, IdxP3Set& dstI, const Vec3& vert, Idx3* poly);		//!< ポリゴン張りなおす際にどのエッジが境界か調べる
			void _check();

			public:
				CnvP(const int (&initial)[4], const Vec3List& vtx);	//!< 必ず4面体から探索を始める
				bool quickHull();
				Idx3List getResult(Vec3List& dstV) const;
		};
		//! 凸包
		class ConvexP : public ITagP<ConvexP> {
			Vec3List	_vtx;
			enum RFLAG {
				RFL_GCENTER = 0x01,
				RFL_POLYFACE = 0x02
			};
			mutable uint32_t	_rflg;
			mutable	Idx3List	_pface;		//! 凸包ポリゴンインデックス
			mutable Vec3		_vGCenter;	//!< 重心座標

			void _init(int) {}
			template <class... Args>
			void _init(int idx, const Vec3& v, Args&&... args) {
				_vtx[idx] = v;
				_init(idx+1, std::forward<Args>(args)...);
			}

			public:
				ConvexP();
				ConvexP(const Vec3* src, int n);
				ConvexP(ConvexP&& c);
				ConvexP(Vec3List&& src);

				template <class... Args>
				ConvexP(const Vec3& v, Args&&... args): _vtx(sizeof...(args)+1) {
					_init(0, v, std::forward<Args>(args)...);
				}
				void addVtx(const Vec3& v);
				const Vec3List& getVtxArray() const;
				const Vec3& getVtx(int n) const;
				int getNVtx() const;
				void popVtx();
				//! quickHull法を用いて凸多面体にする
				/*! ついでにポリゴンインデックスをキャッシュ */
				bool quickHull();
				bool hasPoint(const Vec3& p);
				void swap(ConvexP& c) noexcept;

				static Vec3 DualTransform(const Plane& plane);
				static Plane DualTransform(const Vec3& p);
				static Vec3List DualTransform(const PlaneList& plL);
				Vec3List exportDualTransform();
				void dualTransform(Vec3List& dst, const Vec3& dir, const Vec3& pos);
				ConvexP operator + (const Vec3& ofs) const;
				ConvexP operator - (const Vec3& ofs) const;
				ConvexP& operator += (const Vec3& ofs);
				ConvexP& operator -= (const Vec3& ofs);
				const Idx3List& getPolyFace();

				// ---- cacheable functions ----
				const Vec3& bs_getGCenter() const;
				const Vec3& bs_getCenter() const;
				Sphere bs_getBVolume() const;
				float bs_getArea() const;
				Mat33 bs_getInertia() const;
				AABB bs_getAABB() const;
				// -----------------------------
				Vec3 support(const Vec3& dir) const;

				// すべてGJK法で衝突判定するのでここには判定関数を記述しない
				spn::none_t hit(...) const;
		};
		//! 円柱
		struct Cylinder : ITagP<Cylinder>, Capsule {
			// ---- cacheable functions ----
			Vec3 bs_getGCenter() const;
			Vec3 bs_getCenter() const;
			Sphere bs_getBVolume() const;
			float bs_getArea() const;
			Mat33 bs_getInertia() const;
			// -----------------------------
			Vec3 support(const Vec3& dir) const;
			Cylinder operator * (const AMat43& m) const;

			// 底面と上面を計算し平面を返す
			Vec3x2 getBeginPlaneV() const;
			Vec3x2 getEndPlaneV() const;
			Plane getBeginPlane() const;
			Plane getEndPlane() const;

			// 凸ポリゴンを円柱ローカル座標に変換し上面と底面でクリップ
			void translateLocal(ColCv& cp) const;
			Vec3 translateLocal(const Vec3& vec) const;

			spn::none_t hit(...) const;
			using ITagP<Cylinder>::GetCID;
		};

		//! GJK法による衝突判定(3D)
		/*! ヒットチェックのみ。衝突時は内部点を出力可 */
		class GSimplex {
			enum class Result {
				Default,
				Hit,
				NoHit
			};
			//! 頂点数は最大でも4
			constexpr static int NUM_VERT = 4;

			Vec3	_vtx[NUM_VERT];
			Vec3	_posB[NUM_VERT];
			int		_nVtx;
			Result	_hr;

			const IModel &_m0, &_m1;

			const static Vec3 cs_randVec[6];	//!< 2点以下でNOHITになった時の方向決めに使用
			const static int cs_index[4][4];	//!< 4面体インデックス (最後は使われてない番号)
			int			_retryCount;
			int			_nearIDX;		//!< 4面体のどの面を使ったか
			Vec3		_nearPosB;		//!< _nearIDXが負数の場合はこちらを使う
			Vec3		_nearP;			//!< nearestDir()の結果の座標

			//! minkowski差を, 重複してなければ追加 (不要な頂点の削除)
			bool _addVtx(const Vec3& vA, const Vec3& vB, int id);
			//! 次に探索すべき方向を求める
			/*! @return 算出された方向, 衝突と判定できるか(0=未定, 1=衝突, -1=非衝突), 4面体において使用されなかった頂点 */
			std::tuple<Vec3,Result,int> getNearestDir();
			//! 既存の頂点と重複を調べる
			bool _checkDuplication(const Vec3& v) const;
			void _clear();
			bool _gjkMethod();

			public:
				GSimplex(const IModel& m0, const IModel& m1);
				bool getResult() const;
				Vec3x2 getNearestPair() const;
				Vec3 getInterPoint() const;
				int getNVtx() const;
				const Vec3& getVtx(int n) const;
		};
		//! EPS法で用いる凸包クラス
		class EConvex : public ConvexP {
			static uint32_t _ComposeIDX(int id0, int id1, int id2);
			static std::tuple<int,int,int> _DecompIDX(uint32_t id);

			// 新しい点は末尾に追加していく
			using MapD_ID = std::multimap<float, Idx3>;
			using MapID_D = std::map<uint32_t, float>;	//!< 10bitずつ小さい順に上位ビットから並べる
			MapD_ID		_d_id;	// distance -> index
			MapID_D		_id_d;	// index -> distance

			//! 新たにFaceを追加 (2つのmapに登録)
			void _addPoly(int idx0, int idx1, int idx2);
			void _addPoly(const Idx3List& src);
			//! Faceを削除 (2つのmapから解除)
			void _delPoly(uint32_t id);
			//! MinkowskiSubで新たに凸包の頂点を探して追加 (ダブったら何もしない)
			bool _dividePoly(const IModel& m0, const IModel& m1, const Idx3& poly, const Vec3& dir);

			const IModel	&_m0, &_m1;
			Vec3	_dir;
			float	_dist;

			std::pair<Vec3,float> _epaMethod();

			public:
				EConvex(const GSimplex& gs, const IModel& m0, const IModel& m1);
				std::pair<Vec3, float> getAvoidVector() const;
		};
		//! ミンコフスキー差を求める
		Vec3 MinkowskiSub(const IModel& m0, const IModel& m1, const Vec3& dir);
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
			float d = ls0.dir.dot(ls0.dir) * ls1.dir.dot(ls1.dir) - spn::Square(ls0.dir.dot(ls1.dir));
			if(std::fabs(d) <= Point::NEAR_THRESHOLD) {
				// 2つの直線は平行
				return Vec3x2(ls0.pos, NearestPoint(ls1, ls0.pos, [](float f){return f;}));
			}
			Mat22 m(ls1.dir.dot(ls1.dir), ls0.dir.dot(ls1.dir),
					ls0.dir.dot(ls1.dir), ls0.dir.dot(ls0.dir));
			float fv[2] = {(ls1.pos - ls0.pos).dot(ls0.dir),
							(ls0.pos - ls1.pos).dot(ls1.dir)};
			float fvt[2] = {m.ma[0][0]*fv[0] + m.ma[0][1]*fv[1],
							m.ma[1][0]*fv[0] + m.ma[1][1]*fv[1]};
			fvt[0] *= 1/d;
			fvt[1] *= 1/d;

			return Vec3x2(Vec3(ls0.pos + ls0.dir*clip0(fvt[0])),
							Vec3(ls1.pos + ls1.dir*clip1(fvt[1])));
		}
		#undef DEF_INVALID_BSFUNCS

		struct Types {
			using CTGeo = ::boom::geo3d::CTGeo;
			using MMgr = ::boom::geo3d::ModelMgr;
			using IModel = ::boom::geo3d::IModel;
			using GJK = ::boom::geo3d::GSimplex;
			using Narrow = ::boom::Narrow<Types>;
		};
	}
}
