#pragma once
#include "spinner/matrix.hpp"
#include "spinner/type.hpp"
#include "spinner/misc.hpp"
#include "spinner/assoc.hpp"
#include "spinner/pose.hpp"
#include <cassert>
#include <vector>
#include <memory>

namespace boom {
	using spn::AVec2;
	using spn::Vec2;
	using spn::Vec3;
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
		};

		struct PointCore;
		struct LineCore;
		struct PolyCore;
		struct CircleCore;
		struct ConvexCore;
		using CTGeo = CType<PointCore, LineCore, PolyCore, CircleCore, ConvexCore>;
		template <class T>
		struct IModelP : IModel {
			virtual uint32_t getCID() const override { return CTGeo::Find<T>::result; }
		};

		struct PointCore : Vec2 {
			using Vec2::Vec2;
			float distance(const LineCore& l) const;
			std::pair<Vec2,int> nearest(const LineCore& l) const;
			bool hit(const PointCore& p) const;
		};
		struct Point : PointCore, IModelP<PointCore> {
			using PointCore::PointCore;
			Vec2 support(const Vec2& dir) const override;
			Vec2 center() const override;
		};

		//! 直線
		struct StLineCore {
			Vec2	pos, dir;

			StLineCore() = default;
			StLineCore(const Vec2& p, const Vec2& d);

			Vec2x2 nearest(const StLineCore& st) const;
			Vec2 nearest(const Vec2& p) const;
			float distance(const Vec2& p) const;
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
			enum LINEPOS {
				BEGIN,
				END,
				ONLINE
			};

			Vec2	point[2];

			LineCore() = default;
			LineCore(const Vec2& v0, const Vec2& v1);

			float distance(const LineCore& l) const;
			float length() const;
			float len_sq() const;
			bool hit(const LineCore& l) const;

			using LNear = std::pair<Vec2,LINEPOS>;
			/*! \return Vec2(最寄り座標),LINEPOS(最寄り線分位置) */
			LNear nearest(const Vec2& p) const;
			float ratio(const Vec2& p) const;
			Vec2 support(const Vec2& dir) const;
			StLineCore toStLine() const;
			bool online(const Vec2& p) const;
		};
		struct Line : LineCore, IModelP<LineCore> {
			using LineCore::LineCore;
			Vec2 support(const Vec2& dir) const override;
			Vec2 center() const override;
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
		};
		//! 円の基本クラス
		/*! 共通クラスを当たり判定対応にラップした物 */
		struct Circle : CircleCore, IModelP<CircleCore> {
			using CircleCore::CircleCore;
			Vec2 support(const Vec2& dir) const override;
			Vec2 center() const override;
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

				Vec2 support(const Vec2& dir) const override;
				Vec2 center() const override;
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
			Vec2 support(const Vec2& dir) const override;
			Vec2 center() const override;
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

				Vec2 support(const Vec2& dir) const override;
				Vec2 center() const override;
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

			using PointL = std::vector<Vec2>;
			PointL	point;

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
		};
		struct Convex : ConvexCore, IModelP<ConvexCore> {
			using ConvexCore::ConvexCore;
			Vec2 support(const Vec2& dir) const override;
			Vec2 center() const override;
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
				ConvexModel(const ConvexCore::PointL& pl);
				ConvexModel(ConvexCore::PointL&& pl);

				float getArea() const;
				float getInertia() const;
				const AVec2& getCenter() const;
				const ConvexCore::PointL& getPoint() const;
				ConvexCore::PointL& refPoint();

				void addOffset(const Vec2& ofs);
				Vec2 support(const Vec2& dir) const override;
				Vec2 center() const override;
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

				RPose& refPose();
				const RPose& getPose() const;
		};
		//! 位置と速度を与えた時にかかる加速度(抵抗力)を計算
		class IResist : public std::enable_shared_from_this<IResist> {
			public:
				using sptr = std::shared_ptr<IResist>;
				using csptr = const sptr&;
				struct Accel {
					Vec2	linear;
					float	rot;
				};

			protected:
				sptr	_spNext;		//!< 次の抵抗力計算インタフェース (null可)
				void _callNext(Accel& acc, const RPose::Value& pose) const;

			public:
				virtual ~IResist() {}
				void addNext(csptr sp);
				//! 抵抗力計算の起点
				Accel resist(const RPose::Value& pose) const;
				//! 抵抗力を計算
				/*! \param[in] r 衝突判定クエリの為の剛体マネージャ
					\param[in] pose 現在の姿勢
					\param[inout] accum 抵抗ベクトルの出力先(加算) */
				virtual void resist(Accel& acc, const RPose::Value& pose) const;
		};
		//! 積分インタフェース
		struct IItg {
			using sptr = std::shared_ptr<IItg>;
			using csptr = const sptr&;

			//! 剛体の姿勢を積分計算
			/*! Resistインタフェースは線形リストでつなげる
				\param[inout] st 現在の姿勢
				\param[in] pres 抵抗力計算インタフェース
				\param[in] r 衝突判定クエリの為の剛体マネージャ for IResist
				\param[in] dt 前回フレームからの時間(sec) */
			virtual void advance(RPose::Value& st, const IResist* pres, float dt) = 0;
		};
		//! オイラー積分
		struct Itg_Eular : IItg {
			void advance(RPose::Value& st, const IResist* pres, float dt) override;
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
	};
}
