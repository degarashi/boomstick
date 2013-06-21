#pragma once
#include "spinner/matrix.hpp"
#include "spinner/type.hpp"
#include "spinner/misc.hpp"
#include "spinner/assoc.hpp"
#include "spinner/pose.hpp"
#include "spinner/plane.hpp"
#include <cassert>
#include <vector>
#include <memory>
#include <unordered_map>
#include <boost/optional.hpp>

namespace boom {
	using spn::AVec2;
	using spn::Vec2;
	using spn::Vec3;
	using spn::Plane;
	using spn::AMat32;
	using spn::AMat33;
	using spn::_sseRcp22Bit;
	using spn::CType;
	using Float2 = std::pair<float,float>;
	using Vec2x2 = std::pair<Vec2,Vec2>;

	//! 90度回転行列(2D)
	extern const spn::AMat22 cs_mRot90[2];
	//! v0とv1で表される面積の2倍
	float Area_x2(const Vec2& v0, const Vec2& v1);
	float Area_x2(const Vec3& v0, const Vec3& v1);
	//! 三角形(v0,v1,v2)のvtに対する重心比率
	/*! \return first: (v1-v0)の比率 <br>
				second: (v2-v0)の比率 */
	template <class T>
	Float2 TriangleRatio(const T& v0, const T& v1, const T& v2, const T& vt) {
		Float2 ret;
		T toV1(v1-v0),
			toV2(v2-v0),
			toVT(vt-v0);
		float invarea = _sseRcp22Bit(Area_x2(toV1, toV2));
		ret.first = Area_x2(toV1, toVT) * invarea;		// 横比率
		ret.second = Area_x2(toVT, toV2) * invarea;		// 縦比率
		return ret;
	}
	//! 三角形(v0,v1,v2)のvtに対する重心比率をユーザー定義変数(f0,f1,f2)に適用した物を算出
	/*! \param[in] v0 三角形座標0
		\param[in] f0 ユーザー定義変数
		\return 補間された値 */
	template <class T>
	T TriangleLerp(const Vec2& v0, const Vec2& v1, const Vec2& v2, const Vec2& vt,
					const T& f0, const T& f1, const T& f2)
	{
		auto res = TriangleRatio(v0,v1,v2,vt);
		auto toT01(f0.l_intp(f1, res.second)),
				toT21(f2.l_intp(f1, res.second));
		return toT01.l_intp(toT21, res.first);
	}
	template <class T>
	T LineLerp(const Vec2& v0, const Vec2& v1, const Vec2& vt,
				const T& f0, const T& f1)
	{
		auto line(v1-v0),
			line2(vt-v0);
		float r = line2.length() * _sseRcp22Bit(line.length());
		return spn::Lerp(f0, f1, r);
	}

	namespace geo2d {
		struct PointCore;
		struct LineCore;
		struct PolyCore;
		struct CircleCore;
		struct ConvexCore;
		using CTGeo = CType<PointCore, LineCore, PolyCore, CircleCore, ConvexCore>;

		using PointL = std::vector<Vec2>;
		using LineL = std::vector<LineCore>;
		//! ペナルティ法における抗力と摩擦力
		struct RForce {
			struct F {
				Vec2	linear;
				float	torque;

				F& operator += (const F& f);
			};
			F	sdump,	//!< スプリング & ダンパによる力
				fricD;	//!< 動摩擦力

			RForce& operator += (const RForce& rf);
		};

		struct IModel {
			using sptr = std::shared_ptr<IModel>;
			using csptr = const sptr&;

			//! サポート射像
			/*! 均等でないスケーリングは対応しない、移動は後でオフセット、回転はdirを逆にすれば代用可
				・・との理由で行列変換後の物体に対する射像は無し */
			virtual Vec2 support(const Vec2& dir) const = 0;
			virtual Vec2 center() const = 0;
			//! 形状毎に一意なコリジョン管理IDを取得
			virtual uint32_t getCID() const = 0;
			//! ある座標が図形の内部に入っているか
			virtual bool isInner(const Vec2& pos) const { return false; }
			//! 外郭を構成する頂点で、mdlにめり込んでいる物を抽出
			/*! \return 時計回りでmdlにめり込んでいる頂点リスト。前後も含む */
			virtual PointL getOverlappingPoints(const IModel& mdl, const Vec2& inner) const {
				throw std::runtime_error("not supported function"); }
			//! 境界ボリューム(円)
//			virtual CircleCore getBBCircle() const = 0;
		};
		template <class T>
		struct IModelP : IModel {
			static uint32_t GetCID() { return CTGeo::Find<T>::result; }
			virtual uint32_t getCID() const override { return GetCID(); }
		};
		#define DEF_IMODEL_FUNCS \
			Vec2 support(const Vec2& dir) const override; \
			Vec2 center() const override;
//			CircleCore getBBCircle() const override;

		enum class LINEPOS {
			BEGIN,
			END,
			ONLINE,
			NOTHIT
		};
		using LNear = std::pair<Vec2,LINEPOS>;

		struct PointCore : Vec2 {
			constexpr static float NEAR_THRESHOLD = 1e-5f;

			using Vec2::Vec2;
			using Vec2::distance;
			float distance(const LineCore& l) const;
			LNear nearest(const LineCore& l) const;
			bool hit(const PointCore& p) const;
		};
		struct Point : PointCore, IModelP<PointCore> {
			using PointCore::PointCore;
			DEF_IMODEL_FUNCS
		};
		//! AxisAlignedBox
		struct BoxCore {
			Vec2	minV, maxV;

			BoxCore() = default;
			BoxCore(const Vec2& min_v, const Vec2& max_v);

			CircleCore bbCircle() const;
			Vec2 support(const Vec2& dir) const;
			Vec2 nearest(const Vec2& pos) const;
		};
		struct Box : BoxCore, IModelP<BoxCore> {
			using BoxCore::BoxCore;
			DEF_IMODEL_FUNCS
		};
		class BoxModel : public IModelP<BoxCore> {
			BoxCore		_box;

			public:
				BoxModel() = default;
				BoxModel(const BoxCore& b);
		};

		//! 直線
		struct StLineCore {
			Vec2	pos, dir;

			StLineCore() = default;
			StLineCore(const Vec2& p, const Vec2& d);

			Vec2x2 nearest(const StLineCore& st) const;
			Vec2 nearest(const Vec2& p) const;
			float distance(const Vec2& p) const;
			//! 点を線分上に置く
			Vec2 placeOnLine(const Vec2& p) const;
			//! 基準位置に対する方向ベクトルとの内積
			float posDot(const Vec2& p) const;
		};
		//! 半直線
		struct RayCore {
			Vec2	pos, dir;

			RayCore() = default;
			RayCore(const Vec2& p, const Vec2& d);

			const StLineCore& asStLine() const;
			Vec2x2 nearest(const RayCore& r) const;
			Vec2 nearest(const Vec2& p) const;
		};

		//! 線分
		struct LineCore {
			Vec2	point[2];

			LineCore() = default;
			LineCore(const Vec2& v0, const Vec2& v1);

			float distance(const LineCore& l) const;
			float length() const;
			float len_sq() const;
			bool hit(const LineCore& l) const;

			/*! \return Vec2(最寄り座標),LINEPOS(最寄り線分位置) */
			LNear nearest(const Vec2& p) const;
			//! 線分が交差する位置を調べる
			/*! \return first: 交差座標
						second: 交差しているライン位置 (ONLINE or NOHIT) */
			LNear crossPoint(const LineCore& l) const;
			LNear crossPoint(const StLineCore& l) const;
			float ratio(const Vec2& p) const;
			Vec2 support(const Vec2& dir) const;
			StLineCore toStLine() const;
			bool online(const Vec2& p) const;
		};
		struct Line : LineCore, IModelP<LineCore> {
			using LineCore::LineCore;
			DEF_IMODEL_FUNCS
		};

		//! 円の共通クラス
		/*! 原点中心とした円
			共通クラスには仮想関数テーブル分のポインタを含めたくない為に少々面倒臭い構造になっている
			面積などを求める比較的重い計算メソッドは用意されるがキャッシュなどはしない */
		struct CircleCore {
			Vec2	center;
			float	radius;

			CircleCore() = default;
			CircleCore(const Vec2& c, float r);

			float area() const;
			float inertia() const;			//!< 重心を軸とした慣性モーメント

			Vec2 support(const Vec2& dir) const;
			bool hit(const Vec2& pt) const;
			// radius
			void getArcPoints(PointL& dst, float ang0, float ang1, float deep) const;
		};
		//! 円の基本クラス
		/*! 共通クラスを当たり判定対応にラップした物 */
		struct Circle : CircleCore, IModelP<CircleCore> {
			using CircleCore::CircleCore;
			DEF_IMODEL_FUNCS
			bool isInner(const Vec2& pos) const override;
		};
		//! 円の剛体用クラス
		/*! キャッシュを管理する関係で形状へのアクセスは関数を通す仕様 元クラスはcompososition */
		class CircleModel : public IModelP<CircleCore> {
			CircleCore			_circle;
			mutable float		_area,
								_inertia;
			mutable uint32_t	_rflag;

			public:
				CircleModel();
				explicit CircleModel(const CircleCore& c);
				const CircleCore& getCore() const;

				float area() const;

				const Vec2& getCenter() const;
				float getRadius() const;
				void setCenter(const Vec2& v);
				void setRadius(float r);

				DEF_IMODEL_FUNCS
				bool isInner(const Vec2& pos) const override;
		};

		struct PolyCore {
			Vec2		point[3];

			PolyCore() = default;
			PolyCore(const Vec2& p0, const Vec2& p1, const Vec2& p2);

			float area() const;
			float inertia() const;
			Vec2 center() const;
			bool isInTriangle(const Vec2& p) const;
			std::pair<Vec2,int> nearest(const Vec2& p) const;

			Vec2 support(const Vec2& dir) const;
			void addOffset(const Vec2& ofs);
			static float CalcArea(const Vec2& p0, const Vec2& p1, const Vec2& p2);
			static float CalcArea(const Vec2& p0, const Vec2& p1);
		};
		struct Poly : PolyCore, IModelP<PolyCore> {
			using PolyCore::PolyCore;
			DEF_IMODEL_FUNCS
			bool isInner(const Vec2& pos) const override;
		};
		class PolyModel : IModel {
			PolyCore		_poly;
			mutable float	_area,
							_inertia;
			mutable AVec2	_center;

			enum REFLAG {
				RFL_CENTER = 0x01,
				RFL_AREA = 0x02,
				RFL_INERTIA = 0x04,
				RFL_ALL = 0xff
			};
			mutable uint32_t	_rflag;

			public:
				PolyModel();
				PolyModel(const Vec2& p0, const Vec2& p1, const Vec2& p2);
				PolyModel(const PolyModel& p) = default;

				float getInertia() const;
				float getArea() const;
				const AVec2& getCenter() const;

				void setPoint(int n, const Vec2& v);
				void addOffset(const Vec2& ofs);

				DEF_IMODEL_FUNCS
				bool isInner(const Vec2& pos) const override;
		};

		//! 多角形基本クラス
		struct ConvexCore {
			using AreaL = std::vector<float>;
			struct AreaSum {
				float result;

				AreaSum(): result(0) {}
				void operator()(int n, const Vec2& p0, const Vec2& p1) {
					result += PolyCore::CalcArea(p0,p1);
				}
			};
			struct AreaList {
				AreaL	areaL;
				float	sum;

				AreaList(int n): areaL(n), sum(0) {}
				void operator()(int n, const Vec2& p0, const Vec2& p1) {
					float a = PolyCore::CalcArea(p0,p1);
					areaL[n] = a;
					sum += a;
				}
			};
			PointL	point;
			const static uint8_t cs_index[1<<8];

			ConvexCore() = default;
			/*! \param[in] v 凸包と分かっている頂点 */
			ConvexCore(std::initializer_list<Vec2> v);
			ConvexCore(const PointL& pl);
			ConvexCore(PointL&& pl);

			static ConvexCore FromConcave(const PointL& src);

			float area() const;
			float inertia() const;
			Vec2 center() const;
			/*! 同時に求めると少し効率が良い */
			std::tuple<float,float,Vec2> area_inertia_center() const;
			Vec2 support(const Vec2& dir) const;
			bool isInner(const Vec2& pos) const;
			//! 2つに分割
			/*! \param[out] c0 線分の進行方向左側
				\param[out] c1 線分の進行方向右側 */
			std::pair<ConvexCore, ConvexCore> splitTwo(const StLineCore& l) const;
			//! 2つに分割して左側を自身に適用
			ConvexCore split(const StLineCore& l);
			//! 2つに分割して右側は捨てる
			void splitThis(const StLineCore& l);

			template <class CB>
			void iterate(CB cb) const {
				// 頂点数は3つ以上
				int nL = point.size();
				assert(nL > 2);

				// 先にブリッジの箇所を処理
				cb(nL-1, point.back(), point.front());
				for(int i=0 ; i<nL-1 ; i++)
					cb(i, point[i], point[i+1]);
			}
			void addOffset(const Vec2& ofs);

			enum class POSITION {
				INNER,
				ONLINE,
				OUTER
			};
			using CPos = std::pair<POSITION, int>;
			//! 指定ポイントの内部的な領域IDと内外位置を取得
			/*! \return first=内外判定
						second=領域ID */
			CPos checkPosition(const Vec2& pos) const;
			//! 内部的な通し番号における外郭ライン
			LineCore getOuterLine(int n) const;
			PointL getOverlappingPoints(const IModel& mdl, const Vec2& inner) const;
		};
		struct Convex : ConvexCore, IModelP<ConvexCore> {
			using ConvexCore::ConvexCore;
			DEF_IMODEL_FUNCS
			bool isInner(const Vec2& pos) const override;
			//! 凸包が重なっている領域を求める
			/*! 既に重なっている事が分かっている前提
				\param[in] cnv 他方の物体
				\param[in] inner 重複領域の内部点 */
			Convex getOverlappingConvex(const Convex& cnv, const Vec2& inner) const;
			PointL getOverlappingPoints(const IModel &mdl, const Vec2& inner) const override;
		};
		class ConvexModel : public IModelP<ConvexCore> {
			private:
				ConvexCore		_convex;
				mutable AVec2	_center;
				mutable float	_area,
								_inertia;
				enum REFLAG {
					RFL_CENTER_INERTIA = 0x01,
					RFL_ALL = 0xff
				};
				mutable uint32_t _rflag;
				void _refreshCInertia() const;

			public:
				ConvexModel();
				ConvexModel(const PointL& pl);
				ConvexModel(PointL&& pl);

				float getArea() const;
				float getInertia() const;
				const AVec2& getCenter() const;
				const PointL& getPoint() const;
				PointL& refPoint();

				void addOffset(const Vec2& ofs);
				DEF_IMODEL_FUNCS
				bool isInner(const Vec2& pos) const override;
				PointL getOverlappingPoints(const IModel &mdl, const Vec2& inner) const override;
		};

		//! 剛体制御用
		class RPose : public spn::Pose2D {
			Vec2			_vel,
							_acc;
			GAP_MATRIX_DEF(mutable, _finalInv, 3,3,
			   ((float _rotVel))
			   ((float _rotAcc))
			   ((mutable uint32_t _invAccum))	//!< 逆行列キャッシュを作った時のカウンタ値
			)
			//! Pose2Dのaccumカウンタとは別で速度に関する変数が書き換えられた際にインクリメント
			uint32_t		_velAccum;
			protected:
				void _setAsChanged();
			public:
				struct Value : spn::Pose2D::Value {
					Vec2	&vel, &acc;
					float	&rotVel, &rotAcc;

					Value(RPose& rp);
					~Value();
				};
				struct TValue : spn::Pose2D::TValue {
					Vec2	vel, acc;
					float	rotVel, rotAcc;

					TValue() = default;
					TValue(const Value& v):
						vel(v.vel), acc(v.acc), rotVel(v.rotVel), rotAcc(v.rotAcc) {}
				};

				RPose();
				RPose(const RPose& rp);
				void identity();

				// ---- setter / getter ----
				void setVelocity(const Vec2& v);
				const Vec2& getVelocity() const;
				void setRotVel(float m);
				float getRotVel() const;
				void setAccel(const Vec2& v);
				const Vec2& getAccel() const;
				void setRotAccel(float a);
				float getRotAccel() const;

				//! ある地点の移動方向と速度を調べる
				/*! \param[in] at 速度を調べるワールド座標 */
				Vec2 getVelocAt(const Vec2& at) const;
				//! ある地点の線形速度と回転ベクトルを別々に取得
				/*! \param[in] at 速度を調べるワールド座標 */
				Vec2x2 getVelocities(const Vec2& at) const;
				Vec2 toLocal(const Vec2& wpos) const;
				Vec2 toLocalDir(const Vec2& wdir) const;
				Vec2 toWorld(const Vec2& lpos) const;
				Vec2 toWorldDir(const Vec2& ldir) const;
				uint32_t getVelocityAccum() const;

				const AMat33& getFinalInv() const;
				RPose lerp(const RPose& p1, float t) const;
				Value refValue();
		};

		//! 剛体ラッパ (形状 + 姿勢)
		class Rigid : public IModel {
			public:
				using sptr = std::shared_ptr<Rigid>;
				using csptr = const sptr&;
			private:
				IModel::sptr	_spModel;	//!< 形状
				RPose			_pose;		//!< 姿勢

			public:
				Rigid(IModel::csptr sp);
				Rigid(IModel::csptr sp, const RPose& pose);

				// --- シミュレーションに関する関数など ---
				Vec2 support(const Vec2& dir) const override;
				Vec2 center() const override;
				uint32_t getCID() const override;
				bool isInner(const Vec2& pos) const override;

				RPose& refPose();
				const RPose& getPose() const;
				const IModel& getModel() const;
		};
		//! 位置と速度を与えた時にかかる加速度(抵抗力)を計算
		class IResist : public std::enable_shared_from_this<IResist> {
			public:
				using sptr = std::shared_ptr<IResist>;
				using csptr = const sptr&;

			protected:
				sptr	_spNext;		//!< 次の抵抗力計算インタフェース (null可)
				void _callNext(RForce::F& acc, int index, const RPose::Value& pose) const;

			public:
				virtual ~IResist() {}
				void addNext(csptr sp);
				//! 抵抗力計算の起点
				RForce::F resist(int index, const RPose::Value& pose) const;
				//! 抵抗力を計算
				/*! \param[in] r 衝突判定クエリの為の剛体マネージャ
					\param[in] pose 現在の姿勢
					\param[inout] accum 抵抗ベクトルの出力先(加算) */
				virtual void resist(RForce::F& acc, int index, const RPose::Value& pose) const;
				//! 内部で何か変数を蓄える場合に備えて要素数を伝える
				virtual void setNumOfRigid(int n) {}
		};
		using RPoseOPT = boost::optional<RPose::Value>;
		using TValueA = std::vector<RPose::TValue>;
		//! 積分インタフェース
		struct IItg {
			using sptr = std::shared_ptr<IItg>;
			using csptr = const sptr&;

			//! 剛体の姿勢を積分計算
			/*! Resistインタフェースは線形リストでつなげる
				\param[inout] st 現在の姿勢
				\param[in] pres 抵抗力計算インタフェース
				\param[in] dt 前回フレームからの時間(sec) */
			virtual void advance(int pass, RPoseOPT* value, int nRigid, const IResist* pres, float dt) = 0;
			//! 1回分の時間を進めるのに必要なステップ数
			virtual int numOfIteration() const = 0;

			virtual void beginIteration(int n) {}
			virtual void endIteration() {}
		};
		// ---------------------- 積分アルゴリズム ----------------------
		namespace itg {
			#define DEF_IITG_FUNCS \
				void advance(int pass, RPoseOPT* value, int nRigid, const IResist* pres, float dt) override; \
				int numOfIteration() const override;
			//! オイラー法
			struct Eular : IItg {
				DEF_IITG_FUNCS
			};
			//! 改良版オイラー法
			class ImpEular : public IItg {
				TValueA		_tvalue;
				public:
					DEF_IITG_FUNCS
					void beginIteration(int n) override;
					void endIteration() override;
			};
			//! 4次のルンゲ・クッタ法
			class RK4 : public IItg {
				// 計算途中で必要になる一時変数 [nRigid * 4]
				// Rigid0[0,1,2,3], Rigid1[0,1,2,3] ... と続ける
				TValueA		_tvalue;
				public:
					DEF_IITG_FUNCS
					void beginIteration(int n) override;
					void endIteration() override;
			};
		}
		class RigidMgr;
		// ---------------------- 抵抗力計算 ----------------------
		namespace resist {
			//! 衝突による加速度
			class Impact : public IResist {
				RigidMgr&	_rmgr;
				// 必ず剛体は0番から順に走査するので内部でカウンタを持つ
				std::vector<RForce::F>	_state;

				public:
					Impact(RigidMgr& rm);
					void resist(RForce::F& acc, int index, const RPose::Value &pose) const override;
					void setNumOfRigid(int n) override;
			};
			//! 重力による加速度
			class Gravity : public IResist {
				Vec2	_grav;
				public:
					Gravity(const Vec2& v);
					void resist(RForce::F& acc, int index, const RPose::Value &pose) const override;
			};
		}

		//! コリジョン判定の結果を参照しやすい形で格納
		struct ColResult {
			using PtrArray = std::vector<Rigid*>;
			using CursorArray = std::unordered_map<int, std::pair<uint16_t,uint16_t>>;
			PtrArray		array;
			CursorArray		cursor;

			ColResult() = default;
			ColResult(ColResult&& cr);
			ColResult& operator = (ColResult&& cr);
			ColResult& operator = (const ColResult& cr) = default;
		};
		//! 剛体マネージャ
		class RigidMgr {
			using RList = std::vector<Rigid::sptr>;
			RList			_rlist;
			IItg::sptr		_itg;	//!< 適用する積分アルゴリズム
			IResist::sptr	_res;	//!< 重力などの力

			public:
				RigidMgr(IItg::csptr itg);

				void add(Rigid::csptr sp);
				void add(IResist::csptr sp);
				void simulate(float dt);

				int getNRigid() const;
				Rigid::csptr getRigid(int n) const;

				ColResult checkCollision() const;	// コリジョンチェックして内部変数に格納
		};
		//! ヒットチェックだけ
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
				const Vec2& getInner() const;
		};
		/*! 常に頂点リストを時計回りに保つ */
		class GEpa : public GSimplex {
			std::vector<Vec2x2*>	_vl;
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
			void _printASV();

			union {
				Vec2	_pvec;
				Vec2x2	_nvec;
			};
			void _addAsv(const Vec2& v0, const Vec2& v1, const Vec2x2* (&vtx)[2]);

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
			void _adjustLoop();
			//! 頂点リスト(3つ以上)から最短距離リストを生成
			void _geneASV();

			public:
				GEpa(const IModel& m0, const IModel& m1);
				Vec2x2 getNearestPair() const;
				const Vec2& getPVector() const;
		};

		//! GJK法による衝突判定
		/*! \param[out] dst 内部点の出力先ポインタ(null可)
			\return 衝突していればtrue */
		bool GJKMethod(const IModel& m0, const IModel& m1, Vec2* dst=nullptr);
		//! GJKで最近傍対を求める
		/*! \return 衝突の有無(bool)
					nohit: m0側の最近傍点(Vec2), hit: m0側最深点
					nohit: m1側の最近傍点(Vec2), hit: 最短移動ベクトル */
		std::tuple<bool,Vec2,Vec2> GJKPair(const IModel& m0, const IModel& m1);
		//! ミンコフスキー差を求める
		Vec2 MinkowskiSub(const IModel& m0, const IModel& m1, const Vec2& dir);
		//! 固有のアルゴリズムによる衝突判定
		bool HitCheck(const IModel& mdl0, const IModel& mdl1);
		//! 固有のアルゴリズムでmdlFromのmdlToに対する最深点を算出
		Vec2 HitPos(const IModel& mdlFrom, const IModel& mdlTo);
		//! ペナルティ法における抗力計算
		/*! 重なり領域だけを使うことは無くて、抗力の算出とセットなので内部で積分計算 */
		RForce CalcForce(const Rigid& r0, const Rigid& r1, const Vec2& inner, const StLineCore& div);

		//! DualTransform (point2D -> line2D)
		StLineCore Dual(const Vec2& v);
		//! DualTransform (line2D -> point2D)
		Vec2 Dual(const StLineCore& l);
		//! DualTransform (point3D -> plane)
		Plane Dual(const Vec3& v);
		//! DualTransform (plane -> point3D)
		Vec3 Dual(const Plane& plane);
	}
}
