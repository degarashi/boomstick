//! 3D形状の定義
#pragma once
#include "geom.hpp"
#include "spinner/resmgr.hpp"
#include "spinner/assoc.hpp"
#include <boost/pool/object_pool.hpp>
#include <memory>

namespace boom {
	namespace geo2d {
		struct Point;
		struct Segment;
		struct AABB;
		struct Poly;
		struct Circle;
		struct Convex;
		using CTGeo = spn::CType<Convex, Circle, Poly, AABB, Segment, Point>;
		template <class T>
		using ITagP = ITagP_base<T, CTGeo>;

		using PointL = std::vector<Vec2>;
		using SegL = std::vector<Segment>;
		using Convex2 = std::pair<Convex, Convex>;
		struct Line;
		struct Ray;

		struct Circle : ITagP<Circle> {
			Vec2	vCenter;
			float	fRadius;

			Circle() = default;
			Circle(const Vec2& c, float r);
			// ---- cacheable functions ----
			float bs_getArea() const;
			float bs_getInertia() const;
			const Vec2& bs_getCenter() const;
			const Circle& bs_getBVolume() const;
			// -----------------------------
			Vec2 support(const Vec2& dir) const;
			bool isInner(const Vec2& pos) const;

			spn::none_t hit(...) const;
			bool hit(const Vec2& pt) const;
			bool hit(const Circle& c) const;

			void getArcPoints(PointL& dst, float ang0, float ang1, float deep) const;
			Circle operator * (const AMat32& m) const;
			Circle operator * (float s) const;
		};

		struct IModel : ::boom::IModelNode {
			// ---- cacheable functions ----
			virtual Vec2 im_getCenter() const = 0;
			virtual float im_getArea() const = 0;
			virtual float im_getInertia() const = 0;
			virtual Circle im_getBVolume() const = 0;
			// -----------------------------
			//! サポート射像
			/*! 均等でないスケーリングは対応しない、移動は後でオフセット、回転はdirを逆にすれば代用可
				・・との理由で行列変換後の物体に対する射像は無し */
			virtual Vec2 im_support(const Vec2& dir) const = 0;
			//! ある座標が図形の内部に入っているか
			virtual bool im_isInner(const Vec2& pos) const = 0;
			virtual uint32_t getCID() const = 0;

			virtual Vec2 toLocal(const Vec2& v) const { return v; }
			virtual Vec2 toLocalDir(const Vec2& v) const { return v; }
			virtual Vec2 toWorld(const Vec2& v) const { return v; }
			virtual Vec2 toWorldDir(const Vec2& v) const { return v; }
			virtual spn::Optional<const AMat32&> getToLocal() const { return spn::none; }
			virtual spn::Optional<const AMat32&> getToWorld() const { return spn::none; }
			const static AMat32 cs_idMat;
			const AMat32& getToLocalI() const;
			const AMat32& getToWorldI() const;
		};
		using UPModel = std::unique_ptr<IModel>;
		#define mgr_model2d (::boom::geo2d::ModelMgr::_ref())
		class ModelMgr : public spn::ResMgrA<UPModel, ModelMgr> {};
		DEF_AHANDLE(ModelMgr, Mdl, UPModel, UPModel)

		template <class T>
		struct Model : IModelP_base<T, IModel> {
			using IModelP_base<T, IModel>::IModelP_base;
			// 各種ルーチンの中継
			Circle im_getBVolume() const override { return T::bs_getBVolume(); }
			float im_getInertia() const override { return T::bs_getInertia(); }
			float im_getArea() const override { return T::bs_getArea(); }
			Vec2 im_getCenter() const override { return T::bs_getCenter(); }
			Vec2 im_support(const Vec2& dir) const override { return T::support(dir); }
			bool im_isInner(const Vec2& pos) const override { return T::isInner(pos); }
		};
		//! 座標変換ありのモデル基底
		/*! Modelに被せる形で使う */
		template <class T>
		class TModel : public ITagP<T>, public IModel {
			using VMdl = VModel<ModelMgr>;
			VMdl	_model;
			AMat32	_mToWorld,
					_mToLocal;
			void _calcInv() {
				_mToWorld.convertA33().inversion(reinterpret_cast<AMat33&>(_mToLocal));
			}
			public:
				TModel(const TModel&) = default;
				TModel(const IModel& mdl, const AMat32& mw): _model(mdl), _mToWorld(mw) { _calcInv(); }
				TModel(HMdl hMdl, const AMat32& mw): _model(hMdl), _mToWorld(mw) { _calcInv(); }

				Circle im_getBVolume() const override {
					return toWorld(_model.get().im_getBVolume()); }
				float im_getInertia() const override {
					return _model.get().im_getInertia() * _mToWorld; }
				float im_getArea() const override {
					return _model.get().im_getArea(); }
				Vec2 im_getCenter() const override {
					return toWorld(_model.get().im_getCenter()); }
				Vec2 im_support(const Vec2& dir) const override {
					return toWorldDir(_model.get().im_support(toLocalDir(dir))); }
				bool im_isInner(const Vec2& pos) const override {
					return _model.get().im_isInner(toLocal(pos)); }

				Vec2 toLocal(const Vec2& v) const override {
					return v.asVec3(1) * _mToLocal; }
				Vec2 toLocalDir(const Vec2& v) const override {
					return v.asVec3(0) * _mToLocal; }
				Vec2 toWorld(const Vec2& v) const override {
					return v.asVec3(1) * _mToWorld; }
				Vec2 toWorldDir(const Vec2& v) const override {
					return v.asVec3(0) * _mToWorld; }
				spn::Optional<const AMat32&> getToLocal() const override { return _mToLocal; }
				spn::Optional<const AMat32&> getToWorld() const override { return _mToWorld; }
				uint32_t getCID() const override { return ITagP<T>::GetCID(); }
				const void* getCore() const override { return _model.get().getCore(); }
		};

		using CircleM = Model<Circle>;

		using LNear = std::pair<Vec2, LinePos>;
		struct Point : Vec2, ITagP<Point> {
			constexpr static float NEAR_THRESHOLD = 1e-5f;

			// ---- cacheable functions ----
			const float& bs_getArea() const;
			const float& bs_getInertia() const;
			const Vec2& bs_getCenter() const;
			Circle bs_getBVolume() const;
			// -----------------------------
			using Vec2::Vec2;
			using Vec2::distance;
			float distance(const Segment& s) const;
			LNear nearest(const Segment& s) const;
			const Vec2& support(const Vec2& dir) const;

			spn::none_t hit(...) const;
			bool hit(const Point& p) const;
		};
		using PointM = Model<Point>;

		#define DEF_INVALID_BSFUNCS \
			const float& bs_getArea() const {INVOKE_ERROR} \
			const float& bs_getInertia() const {INVOKE_ERROR} \
			const Vec2& bs_getCenter() const {INVOKE_ERROR} \
			const Circle& bs_getBVolume() const {INVOKE_ERROR}
		//! 直線
		struct Line : ITagP<Line> {
			Vec2	pos, dir;

			Line() = default;
			Line(const Vec2& p, const Vec2& d);

			DEF_INVALID_BSFUNCS
			Vec2x2 nearest(const Line& st) const;
			Vec2 nearest(const Vec2& p) const;
			float distance(const Vec2& p) const;
			//! 点を線分上に置く
			Vec2 placeOnLine(const Vec2& p) const;
			//! 基準位置に対する方向ベクトルとの内積
			float posDot(const Vec2& p) const;
			Line operator * (const AMat32& m) const;

			spn::none_t hit(...) const;
		};
		//! 半直線
		struct Ray : ITagP<Ray> {
			Vec2	pos, dir;

			Ray() = default;
			Ray(const Vec2& p, const Vec2& d);

			DEF_INVALID_BSFUNCS
			const Line& asLine() const;
			Vec2x2 nearest(const Ray& r) const;
			Vec2 nearest(const Vec2& p) const;
			Ray operator * (const AMat32& m) const;

			spn::none_t hit(...) const;
		};
		//! 線分
		struct Segment : ITagP<Segment> {
			Vec2	from, to;

			Segment() = default;
			Segment(const Vec2& v0, const Vec2& v1);
			// ---- cacheable functions ----
			const float& bs_getArea() const {INVOKE_ERROR}
			const float& bs_getInertia() const {INVOKE_ERROR}
			Vec2 bs_getCenter() const;
			Circle bs_getBVolume() const;
			// -----------------------------
			Vec2 support(const Vec2& dir) const;

			float distance(const Segment& l) const;
			float length() const;
			float len_sq() const;
			spn::none_t hit(...) const;
			bool hit(const Segment& l) const;

			/*! \return Vec2(最寄り座標),LINEPOS(最寄り線分位置) */
			LNear nearest(const Vec2& p) const;
			//! 線分が交差する位置を調べる
			/*! \return first: 交差座標
						second: 交差しているライン位置 (ONLINE or NOHIT) */
			LNear crossPoint(const Segment& l) const;
			LNear crossPoint(const Line& l) const;
			float ratio(const Vec2& p) const;
			Vec2 getDir() const;
			Line asLine() const;
			bool online(const Vec2& p) const;

			Segment operator * (const AMat32& m) const;
		};
		//! AxisAlignedBox
		struct AABB : ITagP<AABB> {
			Vec2	minV, maxV;

			AABB() = default;
			AABB(const Vec2& min_v, const Vec2& max_v);
			// ---- cacheable functions ----
			float bs_getArea() const;
			float bs_getInertia() const;
			Vec2 bs_getCenter() const;
			Circle bs_getBVolume() const;
			// -----------------------------
			Vec2 support(const Vec2& dir) const;
			Vec2 nearest(const Vec2& pos) const;

			spn::none_t hit(...) const;
			bool hit(const Segment& l) const;
		};
		struct Poly : ITagP<Poly> {
			Vec2		point[3];

			Poly() = default;
			Poly(const Vec2& p0, const Vec2& p1, const Vec2& p2);
			// ---- cacheable functions ----
			float bs_getArea() const;
			float bs_getInertia() const;
			Vec2 bs_getCenter() const;
			Circle bs_getBVolume() const;
			// -----------------------------
			bool isInTriangle(const Vec2& p) const;
			std::pair<Vec2,int> nearest(const Vec2& p) const;

			Vec2 support(const Vec2& dir) const;
			void addOffset(const Vec2& ofs);
			static float CalcArea(const Vec2& p0, const Vec2& p1, const Vec2& p2);
			static float CalcArea(const Vec2& p0, const Vec2& p1);
			//! 鈍角を探す
			/*! \return 鈍角の番号 (負数は該当なし) */
			int getObtuseCorner() const;

			spn::none_t hit(...) const;
		};
		//! 多角形基本クラス
		struct Convex : ITagP<Convex> {
			using AreaL = std::vector<float>;
			struct AreaSum {
				float result;

				AreaSum(): result(0) {}
				void operator()(int /*n*/, const Vec2& p0, const Vec2& p1) {
					result += Poly::CalcArea(p0,p1);
				}
			};
			struct AreaList {
				AreaL	areaL;
				float	sum;

				AreaList(int n): areaL(n), sum(0) {}
				void operator()(int n, const Vec2& p0, const Vec2& p1) {
					float a = Poly::CalcArea(p0,p1);
					areaL[n] = a;
					sum += a;
				}
			};
			PointL	point;
			const static uint8_t cs_index[1<<8];

			Convex() = default;
			Convex(const Convex& c) = default;
			/*! \param[in] v 凸包と分かっている頂点 */
			Convex(std::initializer_list<Vec2> v);
			Convex(const PointL& pl);
			Convex(PointL&& pl);
			Convex(Convex&& c);
			Convex& operator = (Convex&& c);

			static Convex FromConcave(const PointL& src);
			// ---- cacheable functions ----
			float bs_getArea() const;
			Circle bs_getBVolume() const;
			// -----------------------------
			/*! 同時に求めると少し効率が良い為 */
			std::tuple<float,float,Vec2> area_inertia_center() const;
			Vec2 support(const Vec2& dir) const;
			bool isInner(const Vec2& pos) const;
			//! 2つに分割
			/*! \param[out] c0 線分の進行方向左側
				\param[out] c1 線分の進行方向右側 */
			Convex2 splitTwo(const Line& l) const;
			//! 2つに分割して左側を自身に適用
			Convex split(const Line& l);
			//! 2つに分割して右側は捨てる
			void splitThis(const Line& l);

			template <class CB>
			void iterate(CB cb) const {
				// 頂点数は3つ以上
				int nL = point.size();
				AssertP(Trap, nL > 2);

				// 先にブリッジの箇所を処理
				cb(nL-1, point.back(), point.front());
				for(int i=0 ; i<nL-1 ; i++)
					cb(i, point[i], point[i+1]);
			}
			void addOffset(const Vec2& ofs);

			//! 指定ポイントの内部的な領域IDと内外位置を取得
			/*! \return first=内外判定
						second=領域ID */
			std::pair<ConvexPos, int> checkPosition(const Vec2& pos) const;
			//! 内部的な通し番号における外郭ライン
			Segment getOuterSegment(int n) const;
			Line getOuterLine(int n) const;
			std::pair<bool,PointL> getOverlappingPoints(const Convex& mdl, const Vec2& inner) const;
			static Convex GetOverlappingConvex(const Convex& m0, const Convex& m1, const Vec2& inner);
			//! 凸包が直線と交差している箇所を2点計算
			std::tuple<bool,Vec2,Vec2> checkCrossingLine(const Line& l) const;
			Convex& operator *= (const AMat32& m);
			//! 頂点が時計回りになっているか
			bool checkCW() const;
			//! 頂点の並びを時計回りに修正
			void adjustLoop();
			int getNPoints() const;
			Vec2 getPoint(int n) const;

			spn::none_t hit(...) const;
			std::ostream& dbgPrint(std::ostream& os) const;
		};
		//! GJK法による衝突判定(2D)
		/*! ヒットチェックのみ。衝突時は内部点を出力可 */
		class GSimplex {
			protected:
				const IModel	&_m0, &_m1;
				Vec2	_vtx[3],		//!< 凸包を構成する頂点
						_posB[3],		//!< vtxを求める時に使ったB側のインデックス
						_inner;			//!< 内部点
				bool	_bHit;
				int		_nVtx;
			private:
				void _minkowskiSub(const Vec2& dir, int n);
				bool _gjkMethod();
			public:
				//! 初期化 = GJKによる判定(ヒットチェックだけ)
				GSimplex(const IModel& m0, const IModel& m1);
				bool getResult() const;
				//! 衝突時: 内部点を取得
				const Vec2& getInner() const;
		};
		//! GJKで最近傍対を求める
		/*! 常に頂点リストを時計回りに保つ */
		class GEpa : public GSimplex {
			constexpr static int MAX_VERT = 0x100;
			using VPool = boost::object_pool<Vec2x2>;
			static thread_local VPool tls_vPool;
			using VList = std::array<Vec2x2*, MAX_VERT>;

			VList	_vl;
			size_t	_szVl;
			Vec2x2* _allocVert(int n=-1);
			int _getIndex(const Vec2x2* vp) const;

			//! 点と辺の両対応
			struct LmLen {
				float			dist;
				Vec2			dir;
				const Vec2x2	*vtx[2];	//!< index[1]==nullptrの時は単一の頂点を表す

				bool operator < (const LmLen& len) const {
					return dist < len.dist;
				}
			};
			//! 最短距離リスト
			spn::AssocVec<LmLen>	_asv;
			void _printASV() const;

			union {
				Vec2	_pvec;
				Vec2x2	_nvec;
			};
			bool _addAsv(const Vec2x2& v0, const Vec2x2& v1);

			/*! \param[in] n 計算した頂点の挿入先インデックス */
			const Vec2x2& _minkowskiSub(const Vec2& dir, int n=-1);
			//! Hit時の脱出ベクトル
			/*! 最低でも3頂点以上持っている前提 */
			void _epaMethodOnHit();
			//! NoHit時の最短ベクトル
			void _epaMethodNoHit();

			//! NoHit: 頂点が2つしか無い時の補正
			void _recover2NoHit();
			//! Hit: 頂点が2つしか無い時の補正
			void _recover2OnHit();
			//! 頂点の並び順を時計回りに修正
			void _adjustLoop3();
			//! 頂点リスト(3つ以上)から最短距離リストを生成
			void _geneASV();
			void _clear();

			public:
				GEpa(const IModel& m0, const IModel& m1);
				~GEpa();
				/*! 非衝突時に有効
					\return A側の最近傍点, B側の最近傍点 */
				Vec2x2 getNearestPair() const;
				//! 衝突を回避するための最短移動ベクトル(A側)
				const Vec2& getPVector() const;
		};
		//! ミンコフスキー差を求める
		Vec2 MinkowskiSub(const IModel& m0, const IModel& m1, const Vec2& dir);
		//! DualTransform (point2D -> line2D)
		Line Dual(const Vec2& v);
		//! DualTransform (line2D -> point2D)
		Vec2 Dual(const Line& ls);
		Vec2 Dual2(const Vec2& v0, const Vec2& v1);
		//! DualTransform (point3D -> plane)
		Plane Dual(const Vec3& v);
		//! DualTransform (plane -> point3D)
		Vec3 Dual(const Plane& plane);

		template <class CLIP>
		inline Vec2 NearestPoint(const Line& ls, const Vec2& p, CLIP clip) {
			Vec2 toP = p - ls.pos;
			float d = ls.dir.dot(toP);
			return ls.pos + ls.dir * clip(d);
		}
		template <class CLIP>
		inline Vec2x2 NearestPoint(const Line& ls0, const Line& ls1, CLIP clip0, CLIP clip1) {
			float st0d = ls0.dir.len_sq(),
					st1d = ls1.dir.len_sq(),
					st01d = ls0.dir.dot(ls1.dir);
			float d = st0d * st1d - spn::Square(st01d);
			if(std::fabs(d) < 1e-5f) {
				// 2つの直線は平行
				return Vec2x2(ls0.pos, NearestPoint(ls1, ls0.pos, [](float f){return f;}));
			}
			Mat22 m0(st1d, st01d,
					st01d, st0d);
			Vec2	m1((ls1.pos - ls0.pos).dot(ls0.dir),
						(ls0.pos - ls1.pos).dot(ls1.dir));
			m1 = m0 * m1;
			m1 *= spn::Rcp22Bit(d);
			return Vec2x2(ls0.pos + ls0.dir * clip0(m1.x),
							ls1.pos + ls1.dir * clip1(m1.y));
		}
		#undef DEF_INVALID_BSFUNCS

		struct Types {
			using CTGeo = ::boom::geo2d::CTGeo;
			using MMgr = ::boom::geo2d::ModelMgr;
			using IModel = ::boom::geo2d::IModel;
			using GJK = ::boom::geo2d::GSimplex;
			using Narrow = ::boom::Narrow<Types>;
		};
	}
}
