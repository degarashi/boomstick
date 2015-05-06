//! 2D形状の定義
#pragma once
#include "geom.hpp"
#include "spinner/resmgr.hpp"
#include "spinner/assoc.hpp"
#include <boost/pool/object_pool.hpp>
#include <memory>
#include "spinner/structure/treenode.hpp"
#include "spinner/rflag.hpp"
#include "boundary.hpp"

namespace boom {
	namespace geo2d {
		using Int_OP = spn::Optional<int>;
		struct Point;
		struct Segment;
		struct AABB;
		struct Poly;
		struct Circle;
		struct Convex;
		struct Ray;
		struct Line;
		using CTGeo = spn::CType<Convex, Circle, Poly, AABB, Segment, Ray, Line, Point>;
		template <class T>
		using ITagP = ITagP_base<T, CTGeo>;

		using PointL = std::vector<Vec2>;
		using SegL = std::vector<Segment>;
		using Convex2 = std::pair<Convex, Convex>;
		struct Line;
		struct Ray;

		struct IModel;
		struct Circle : ITagP<Circle> {
			Vec2	vCenter;
			float	fRadius;

			Circle() = default;
			Circle(const Vec2& c, float r);

			// -----------------------------
			float bs_getArea() const;
			float bs_getInertia() const;
			const Vec2& bs_getCenter() const;
			const Circle& bs_getBVolume() const;
			// -----------------------------
			Vec2 support(const Vec2& dir) const;

			spn::none_t hit(...) const;
			bool hit(const Vec2& p, float t=NEAR_THRESHOLD) const;
			bool hit(const Segment& s, float t=NEAR_THRESHOLD) const;
			bool hit(const Circle& c, float t=NEAR_THRESHOLD) const;

			void getArcPoints(PointL& dst, float ang0, float ang1, float deep) const;
			Circle operator * (const AMat32& m) const;
			Circle operator * (float s) const;
			Circle& operator += (const Vec2& ofs);
			void distend(float width, float mindist);

			// ---- for MakeBoundary ----
			void setBoundary(const IModel* p);
			void appendBoundary(const IModel* p);

			friend std::ostream& operator << (std::ostream& os, const Circle& c);
		};
		//! IModelインタフェースの子ノードイテレータ
		struct MdlItr {
			TfBase_SP	_sp;

			MdlItr() = default;
			MdlItr(const TfBase_SP& sp);

			MdlItr& operator ++ ();
			bool operator == (const MdlItr& m) const;
			bool operator != (const MdlItr& m) const;
			explicit operator bool () const;
			const TfBase* get() const;
		};
		struct IModel : ::boom::IModelNode {
			// -----------------------------
			virtual Vec2 im_getCenter() const = 0;
			virtual float im_getArea() const = 0;
			virtual float im_getInertia() const = 0;
			virtual Circle im_getBVolume() const = 0;
			// -----------------------------
			virtual void im_transform[[noreturn]](void* dst, const AMat32& m) const = 0;
			//! サポート射像
			/*! 均等でないスケーリングは対応しない、移動は後でオフセット、回転はdirを逆にすれば代用可
				・・との理由で行列変換後の物体に対する射像は無し */
			virtual Vec2 im_support(const Vec2& dir) const = 0;
			//! 図形と点の判定
			virtual bool im_hitPoint(const Vec2& p, float t=NEAR_THRESHOLD) const = 0;
			virtual uint32_t getCID() const = 0;
			virtual MdlItr getInner() const;
			virtual Model_SP im_clone() const = 0;

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
		#define mgr_model2d (::boom::geo2d::ModelMgr::_ref())
		class ModelMgr : public spn::ResMgrA<Model_SP, ModelMgr> {};

		//! TにIModelインタフェースを付加
		template <class T>
		struct Model : IModelP_base<T, IModel> {
			using base_t = IModelP_base<T, IModel>;
			using base_t::base_t;

			Model() = default;
			Model(const T& t): base_t(t) {}
			Model(T&& t): base_t(std::move(t)) {}
			// 各種ルーチンの中継
			Model_SP im_clone() const override {
				return std::make_shared<Model<T>>(static_cast<const T&>(*this));
			}
			void im_transform(void* dst, const AMat32& m) const override { *reinterpret_cast<T*>(dst) = *this * m; }
			Circle im_getBVolume() const override { return T::bs_getBVolume(); }
			float im_getInertia() const override { return T::bs_getInertia(); }
			float im_getArea() const override { return T::bs_getArea(); }
			Vec2 im_getCenter() const override { return T::bs_getCenter(); }
			Vec2 im_support(const Vec2& dir) const override { return T::support(dir); }
			bool im_hitPoint(const Vec2& pos, float t=NEAR_THRESHOLD) const override { return T::hit(pos, t); }
			std::ostream& print(std::ostream& os) const override { return os << static_cast<const T&>(*this); }
		};
		template <class T>
		std::ostream& operator << (std::ostream& os, const Model<T>& m) {
			return m.print(os);
		}
		using CircleM = Model<Circle>;

		class TfBase : public spn::TreeNode<TfBase>,
						virtual public IModel
		{
			private:
				friend class spn::TreeNode<TfBase>;
				using base_t = spn::TreeNode<TfBase>;
			protected:
				void* _getUserData(void*, std::true_type);
				void* _getUserData(void* udata, std::false_type);
				virtual void onChildAdded(const SP& node);
				virtual void onChildRemove(const SP& node);
				virtual void onParentChange(const SP& from, const SP& to);
			public:
				//! 子ノードの取得
				MdlItr getInner() const override;
				bool hasInner() const override;
				virtual bool isLeaf() const { return false; }
				virtual spn::Pose2D& tf_refPose[[noreturn]]();
				virtual const spn::Pose2D& tf_getPose() const;
				Model_SP im_clone() const override;
				virtual TfBase_SP clone() const = 0;
		};
		template <class Boundary, class Ud>
		class TfNode_base : public TfBase,
							public Model<Boundary>
		{
			protected:
				using model_t = Model<Boundary>;
				Ud					_udata;
			public:
				using model_t::model_t;
				Model_SP im_clone() const override {
					return TfBase::im_clone();
				}
				/*! ユーザーデータがvoidの時は親ノードのデータを返す */
				void* getUserData() override {
					return _getUserData(&_udata, std::is_same<spn::none_t, Ud>());
				}
		};
		template <class Boundary, class Ud=spn::none_t>
		class TfNode_Static : public TfNode_base<Boundary, Ud> {
			private:
				using base_t = TfNode_base<Boundary, Ud>;
				using SP = typename base_t::SP;
				using Time_t = typename base_t::Time_t;
				mutable bool	_bChanged = true,
								_bValid;
			protected:
				void onChildAdded(const TfBase::SP& node) override {
					_bChanged = true;
					base_t::onChildAdded(node);
				}
				void onChildRemove(const TfBase::SP& node) override {
					_bChanged = true;
					base_t::onChildRemove(node);
				}
			public:
				SP clone() const override {
					return std::make_shared<TfNode_Static>(*this);
				}
				bool imn_refresh(Time_t t) const override {
					if(_bChanged) {
						_bChanged = false;
						std::vector<const IModel*> pm;
						// 直下の子ノードをリストアップ
						this->template iterateDepthFirst<false>([&pm](auto& node, int depth){
							if(depth == 0)
								return TfBase::Iterate::StepIn;
							pm.push_back(&node);
							return TfBase::Iterate::Next;
						});
						if((_bValid = !pm.empty())) {
							// 境界ボリュームの更新
							auto* core = const_cast<Boundary*>(reinterpret_cast<const Boundary*>(base_t::getCore()));
							auto op = MakeBoundary<Boundary>(&pm[0], pm.size(), t);
							if((_bValid = op))
								*core = *op;
						}
						// 子ノードが無い場合は常にヒットしない
					}
					return _bValid;
				}
		};
		template <class Boundary, class Ud=spn::none_t>
		class TfNode_Dynamic : public TfNode_base<Boundary, Ud> {
			private:
				using base_t = TfNode_base<Boundary, Ud>;
				using Time_t = typename base_t::Time_t;
				using SP = typename base_t::SP;
				using Time_OP = spn::Optional<Time_t>;
				mutable Time_OP		_opTime = 0;
				mutable bool		_bValid;
			protected:
				void onChildAdded(const SP& node) override {
					_opTime = spn::none;
					base_t::onChildAdded(node);
				}
				void onChildRemove(const SP& node) override {
					_opTime = spn::none;
					base_t::onChildRemove(node);
				}
			public:
				SP clone() const override {
					return std::make_shared<TfNode_Dynamic>(*this);
				}
				bool imn_refresh(Time_t t) const override {
					// 更新時刻が古い時のみ処理を行う
					if(!_opTime || *_opTime < t) {
						_opTime = t;

						std::vector<const IModel*> pm;
						// 子ノードをリストアップ
						this->template iterateDepthFirst<false>([&pm](auto& node, int depth){
							if(depth == 0)
								return TfBase::Iterate::StepIn;
							pm.push_back(&node);
							return TfBase::Iterate::Next;
						});
						if((_bValid = !pm.empty())) {
							// 境界ボリュームの更新
							auto* core = const_cast<Boundary*>(reinterpret_cast<const Boundary*>(base_t::getCore()));
							auto op = MakeBoundary<Boundary>(&pm[0], pm.size(), t);
							if((_bValid = op))
								*core = *op;
						}
						// 子ノードが無い場合は常にヒットしない
					}
					return _bValid;
				}
		};

		#define mgr_tf2d (::boom::geo2d::TfMgr::_ref())
		class TfMgr : public spn::ResMgrA<TfBase_SP, TfMgr> {};

		using LNear = std::pair<Vec2, LinePos>;
		struct Point : Vec2, ITagP<Point> {
			// -----------------------------
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
			Point operator * (const AMat32& m) const;

			spn::none_t hit(...) const;
			bool hit(const Vec2& p, float t=NEAR_THRESHOLD) const;
			friend std::ostream& operator << (std::ostream& os, const Point& p);
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

			// -----------------------------
			float bs_getArea() const;
			float bs_getInertia() const;
			const Vec2& bs_getCenter() const;
			const Circle& bs_getBVolume() const;
			// -----------------------------

			Vec2x2 nearest(const Line& st) const;
			Vec2 nearest(const Vec2& p) const;
			float distance(const Vec2& p) const;
			//! 点を線分上に置く
			Vec2 placeOnLine(const Vec2& p) const;
			//! 基準位置に対する方向ベクトルとの内積
			float posDot(const Vec2& p) const;
			Line operator * (const AMat32& m) const;
			LineDivision checkSide(const Vec2& p, float t=DOT_THRESHOLD) const;

			Vec2 support(const Vec2& dir) const;
			spn::none_t hit(...) const;
			bool hit(const Vec2& p, float t=NEAR_THRESHOLD) const;
			friend std::ostream& operator << (std::ostream& os, const Line& l);
		};
		using LineM = Model<Line>;
		//! 半直線
		struct Ray : ITagP<Ray> {
			Vec2	pos, dir;

			Ray() = default;
			Ray(const Vec2& p, const Vec2& d);

			// -----------------------------
			float bs_getArea() const;
			float bs_getInertia() const;
			const Vec2& bs_getCenter() const;
			const Circle& bs_getBVolume() const;
			// -----------------------------

			const Line& asLine() const;
			Vec2x2 nearest(const Ray& r) const;
			Vec2 nearest(const Vec2& p) const;
			Ray operator * (const AMat32& m) const;

			Vec2 support(const Vec2& dir) const;
			spn::none_t hit(...) const;
			bool hit(const Vec2& p, float t=NEAR_THRESHOLD) const;
			friend std::ostream& operator << (std::ostream& os, const Ray& r);
		};
		using RayM = Model<Ray>;
		//! 線分
		struct Segment : ITagP<Segment> {
			Vec2	from, to;

			Segment() = default;
			Segment(const Vec2& v0, const Vec2& v1);
			// -----------------------------
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
			bool hit(const Vec2& p, float t=NEAR_THRESHOLD) const;
			bool hit(const Segment& s, float t=NEAR_THRESHOLD) const;

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
			friend std::ostream& operator << (std::ostream& os, const Segment& s);
		};
		using SegmentM = Model<Segment>;
		//! AxisAlignedBox
		struct AABB : ITagP<AABB> {
			Vec2	minV, maxV;

			int _getAreaNumX(float p) const;
			int _getAreaNumY(float p) const;
			void _makeSegmentX(Segment& s, int num) const;
			void _makeSegmentY(Segment& s, int num) const;
			bool _checkHitX(const Segment& s0, int from, int to) const;
			bool _checkHitY(const Segment& s0, int from, int to) const;
			std::pair<int,int> _getAreaNum(const Vec2& p) const;

			AABB() = default;
			AABB(const Vec2& min_v, const Vec2& max_v);
			// -----------------------------
			float bs_getArea() const;
			float bs_getInertia() const;
			Vec2 bs_getCenter() const;
			Circle bs_getBVolume() const;
			// -----------------------------
			Vec2 support(const Vec2& dir) const;
			Vec2 nearest(const Vec2& pos) const;

			spn::none_t hit(...) const;
			bool hit(const Vec2& p, float t=NEAR_THRESHOLD) const;
			bool hit(const Segment& l, float t=NEAR_THRESHOLD) const;
			bool hit(const AABB& ab, float t=NEAR_THRESHOLD) const;
			AABB operator * (const AMat32& m) const;
			AABB& operator += (const Vec2& ofs);
			void distend(float width, float mindist);

			// ---- for MakeBoundary ----
			void setBoundary(const IModel* p);
			void appendBoundary(const IModel* p);

			friend std::ostream& operator << (std::ostream& os, const AABB& a);
		};
		using AABBM = Model<AABB>;
		struct Poly : ITagP<Poly> {
			Vec2		point[3];
			bool _isInTriangle(const Vec2& p, float threshold) const;

			Poly() = default;
			Poly(const Vec2& p0, const Vec2& p1, const Vec2& p2);
			// -----------------------------
			float bs_getArea() const;
			float bs_getInertia() const;
			Vec2 bs_getCenter() const;
			Circle bs_getBVolume() const;
			// -----------------------------
			//! 座標が三角形の内部に含まれるかを判定
			/*! ポリゴンは時計回りを想定
				辺上は含まない */
			bool isInTriangle(const Vec2& p) const;
			std::pair<Vec2,int> nearest(const Vec2& p) const;

			Vec2 support(const Vec2& dir) const;
			void addOffset(const Vec2& ofs);
			static float CalcArea(const Vec2& p0, const Vec2& p1, const Vec2& p2);
			static float CalcArea(const Vec2& p0, const Vec2& p1);
			//! 鈍角を探す
			/*! \return 鈍角の番号 (負数は該当なし) */
			int getObtuseCorner() const;
			//! 頂点の並びが時計回りかを判定
			bool isCW() const;
			//! 頂点の並びを反転
			void invert();
			LineDivision checkSide(const Line& l, float t=DOT_THRESHOLD) const;

			spn::none_t hit(...) const;
			//! 座標が三角形と衝突するか判定
			/*! isInTriangleとは異なり辺上もHitとみなす */
			bool hit(const Vec2& p, float t=NEAR_THRESHOLD) const;
			bool hit(const Poly& p, float t=NEAR_THRESHOLD) const;
			Poly operator * (const AMat32& m) const;
			Poly& operator += (const Vec2& ofs);
			void distend(float width, float mindist);
			friend std::ostream& operator << (std::ostream& os, const Poly& p);
		};
		using PolyM = Model<Poly>;
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
			// -----------------------------
			float bs_getArea() const;
			Circle bs_getBVolume() const;
			Vec2 bs_getCenter() const;
			float bs_getInertia() const;
			// -----------------------------
			/*! 同時に求めると少し効率が良い為 */
			std::tuple<float,float,Vec2> area_inertia_center() const;
			Vec2 support(const Vec2& dir) const;
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
			std::pair<ConvexPos, int> checkPosition(const Vec2& pos, float threshold=DOT_THRESHOLD) const;
			//! 内部的な通し番号における外郭ライン
			Segment getOuterSegment(int n) const;
			Line getOuterLine(int n) const;
			std::pair<bool,PointL> getOverlappingPoints(const Convex& mdl, const Vec2& inner) const;
			static Convex GetOverlappingConvex(const Convex& m0, const Convex& m1, const Vec2& inner);
			//! 凸包が直線と交差している箇所を2点計算
			std::tuple<bool,Vec2,Vec2> checkCrossingLine(const Line& l) const;
			Convex operator * (const AMat32& m) const;
			Convex& operator *= (const AMat32& m);
			Convex& operator += (const Vec2& ofs);
			void distend(float width, float mindist);
			//! 頂点が時計回りになっているか
			bool checkCW() const;
			//! 頂点の並びを時計回りに修正
			void adjustLoop();
			int getNPoints() const;
			Vec2 getPoint(int n) const;

			spn::none_t hit(...) const;
			bool hit(const Vec2& p, float t=NEAR_THRESHOLD) const;
			friend std::ostream& operator << (std::ostream& os, const Convex& c);
		};
		using ConvexM = Model<Convex>;
		//! GJK法による衝突判定(2D)
		/*! ヒットチェックのみ。衝突時は内部点を出力可 */
		class GSimplex {
			protected:
				const IModel	&_m0, &_m1;
				Poly	_poly;			//!< 凸包を構成するポリゴン
				Vec2	_posA[3],		//!< vtx(A-B)を求める時に使ったA側の座標
						_inner;			//!< 内部点
				bool	_bHit;			//!< 衝突の有無
				int		_nVtx;			//!< 使用された頂点の数(min=1, max=3)
			private:
				void _minkowskiSub(const Vec2& dir, int n);
				void _gjkMethod();
				void _setAsHit(int nv, const Vec2& inner);
				void _setAsNotHit(int nv);
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
			/*! ミンコフスキー差の凸形状を構成する最大頂点数(VList)
				必要ならvectorでやっても良いがとりあえず決め打ち */
			constexpr static int MAX_VERT = 0x100;
			using VPool = boost::object_pool<Vec2x2>;
			static thread_local VPool tls_vPool;
			using VList = std::array<Vec2x2*, MAX_VERT>;

			VList	_vl;
			size_t	_szVl;
			//! 新しく頂点メモリを確保
			/*! \param n 格納先のインデックス (負数なら末尾) */
			Vec2x2* _allocVert(int n=-1);
			//! 頂点ポインタからのインデックス検索
			/*! \param vp 検索する頂点ポインタ */
			Int_OP _getIndex(const Vec2x2* vp) const;

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
			//! デバッグ用
			void _printASV(std::ostream& os) const;

			union {
				Vec2x2	_pvec;
				Vec2x2	_nvec;
			};
			//! v0.firstとv1.firstからなる線分候補をリストに追加
			/*!	\return 最近傍点が原点と重なっていれば0x02, 線分候補が追加されれば0x01, それ以外は0x00 */
			int _addAsv(const Vec2x2& v0, const Vec2x2& v1);

			//! 指定方向へのミンコフスキー差
			/*! \param[in] n 計算した頂点の挿入先インデックス */
			const Vec2x2& _minkowskiSub(const Vec2& dir, int n=-1);
			//! Hit時の脱出ベクトル
			/*! 最低でも3頂点以上持っている前提 */
			void _epaMethodOnHit(float threshold);
			//! NoHit時の最短ベクトル
			void _epaMethodNoHit(float threshold);

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
				GEpa(const IModel& m0, const IModel& m1, float threshold=NEAR_THRESHOLD);
				~GEpa();
				/*! 非衝突時に有効
					\return A側の最近傍点, B側の最近傍点 */
				const Vec2x2& getNearestPair() const;
				/*! 衝突時にそれを回避するための最短移動ベクトル(A側)
					\return first=A側の最深点 second=A側の回避ベクトル */
				const Vec2x2& getPVector() const;
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
		bool IsCrossing(const Line& ls0, const Line& ls1, float len0, float len1, float t=NEAR_THRESHOLD);
		template <class CLIP0, class CLIP1>
		inline Vec2x2 NearestPoint(const Line& ls0, const Line& ls1, CLIP0 clip0, CLIP1 clip1) {
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
