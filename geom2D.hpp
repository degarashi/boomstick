#pragma once
#include "spinner/matrix.hpp"
#include "spinner/type.hpp"
#include "spinner/misc.hpp"
#include "spinner/assoc.hpp"
#include "spinner/pose.hpp"
#include "spinner/plane.hpp"
#include "spinner/resmgr.hpp"
#include "collision.hpp"
#include "cache2D.hpp"
#include <cassert>
#include <vector>
#include <memory>

namespace boom {
	using spn::AVec2;
	using spn::Vec2;
	using spn::Vec3;
	using spn::Plane;
	using spn::AMat32;
	using spn::Mat32;
	using spn::Mat33;
	using spn::AMat33;
	using spn::_sseRcp22Bit;
	using spn::CType;
	using spn::CramersRule;
	using spn::CramerDet;
	using Float2 = std::pair<float,float>;
	using int2 = std::pair<int,int>;
	using int2x2 = std::pair<int2,int2>;
	using uint2 = std::pair<uint32_t,uint32_t>;
	using Vec2x2 = std::pair<Vec2,Vec2>;

	//! 90度回転行列(2D)
	extern const spn::AMat22 cs_mRot90[2];
	//! v0とv1で表される面積の2倍
	float Area_x2(const Vec2& v0, const Vec2& v1);
	float Area_x2(const Vec3& v0, const Vec3& v1);
	//! 三角形(v0,v1,v2)のvtに対する重心比率
	/*! \return first: (v1-v0)の比率 <br>
				second: (v2-v0)の比率 */
	inline Vec2 TriangleRatio(const Vec2& v0, const Vec2& v1, const Vec2& v2, const Vec2& vt) {
		Float2 ret;
		Vec2 	toV1(v1-v0),
				toV2(v2-v0),
				toVT(vt-v0);
		float det = spn::CramerDet(toV1, toV2);
		return spn::CramersRule(toV1, toV2, toVT, _sseRcp22Bit(det));
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
		auto toT01 = f1-f0,
			toT02 = f2-f0;
		return (toT01 * res.x) + (toT02 * res.y) + f0;
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
		struct BoxCore;
		struct PolyCore;
		struct CircleCore;
		struct ConvexCore;
		using CTGeo = CType<PointCore, LineCore, BoxCore, PolyCore, CircleCore, ConvexCore>;

		using PointL = std::vector<Vec2>;
		using LineL = std::vector<LineCore>;
		struct StLineCore;
		struct Convex;
		// class Rigid;
		using Convex2 = std::pair<Convex,Convex>;
		using ConvexCore2 = std::pair<ConvexCore, ConvexCore>;
		struct IModel {
			enum class POSITION {
				INNER,
				ONLINE,
				OUTER
			};
			using CPos = std::pair<POSITION, int>;
			using PosL = std::pair<bool, PointL>;

			//! サポート射像
			/*! 均等でないスケーリングは対応しない、移動は後でオフセット、回転はdirを逆にすれば代用可
				・・との理由で行列変換後の物体に対する射像は無し */
			virtual Vec2 support(const Vec2& dir) const = 0;
			virtual Vec2 getCenter() const = 0;
			//! 形状毎に一意なコリジョン管理IDを取得
			virtual uint32_t getCID() const = 0;
			//! ある座標が図形の内部に入っているか
			virtual bool isInner(const Vec2& /*pos*/) const { return false; }

			#define INVOKE_ERROR throw std::runtime_error(std::string("not supported function: ") + __func__);
			//! 外郭を構成する頂点で、mdlにめり込んでいる物を抽出
			/*! \return 時計回りでmdlにめり込んでいる頂点リスト。前後も含む */
			virtual PosL getOverlappingPoints(const IModel& /*mdl*/, const Vec2& /*inner*/) const { INVOKE_ERROR }
			//! 境界ボリューム(円)
			virtual CircleCore getBCircle() const;

			virtual float getArea(bool /*bInv*/) const { INVOKE_ERROR }
			virtual float getInertia(bool /*bInv*/) const { INVOKE_ERROR }

			#define DEF_CASTFUNC(typ) virtual boost::optional<typ&> castAs##typ() { return boost::none; } \
				virtual boost::optional<const typ&> castAs##typ() const { auto ref = const_cast<IModel*>(this)->castAs##typ(); \
					if(ref) return *ref; return boost::none; }
			// DEF_CASTFUNC(Rigid)

			// ---------- Convex専用 ----------
			//! 物体を構成する頂点数を取得
			/*! Circle, Lineなど物によってはサポートしていない */
			virtual int getNPoints() const { INVOKE_ERROR }
			virtual Vec2 getPoint(int /*n*/) const { INVOKE_ERROR }
			virtual CPos checkPosition(const Vec2& /*pos*/) const { INVOKE_ERROR }
			virtual Convex2 splitTwo(const StLineCore& line) const;
			virtual std::ostream& dbgPrint(std::ostream& /*os*/) const { INVOKE_ERROR }
			friend std::ostream& operator << (std::ostream& os, const IModel& mdl);
		};
		#define mgr_model	ModelMgr::_ref()
		class ModelMgr : public spn::ResMgrA<std::unique_ptr<IModel>, ModelMgr> {};
		DEF_HANDLE(ModelMgr, Mdl, std::unique_ptr<IModel>)

		template <class T>
		struct IModelP : IModel {
			static uint32_t GetCID() { return CTGeo::Find<T>::result; }
			virtual uint32_t getCID() const override { return GetCID(); }
		};
		#define DEF_IMODEL_FUNCS \
			Vec2 support(const Vec2& dir) const override; \
			Vec2 getCenter() const override;
		#define DEF_IMODEL_MASS \
			float getArea(bool bInv=false) const override; \
			float getInertia(bool bInv=false) const override; \
			CircleCore getBCircle() const override;

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
			CircleCore getBCircle() const override;
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
			CircleCore getBCircle() const;
		};
		struct Line : LineCore, IModelP<LineCore> {
			using LineCore::LineCore;
			DEF_IMODEL_FUNCS
			CircleCore getBCircle() const override;
		};

		//! 円の共通クラス
		/*! 原点中心とした円
			共通クラスには仮想関数テーブル分のポインタを含めたくない為に少々面倒臭い構造になっている
			面積などを求める比較的重い計算メソッドは用意されるがキャッシュなどはしない */
		struct CircleCore {
			using CT = CType<TagArea, TagInertia>;
			template <class T>
			using Wrap = CoreRaw<T>;

			float area() const;
			float inertia() const;				//!< 重心を軸とした慣性モーメント
			const Vec2& center() const;
			const CircleCore& bcircle() const;
			// --------------------------------------------
			Vec2	vCenter;
			float	fRadius;

			CircleCore() = default;
			CircleCore(const Vec2& c, float r);

			Vec2 support(const Vec2& dir) const;
			bool hit(const Vec2& pt) const;
			bool hit(const CircleCore& c) const;
			// radius
			void getArcPoints(PointL& dst, float ang0, float ang1, float deep) const;
			CircleCore operator * (const AMat32& m) const;
		};
		//! 円の基本クラス
		/*! 共通クラスを当たり判定対応にラップした物 */
		struct Circle : CircleCore, IModelP<CircleCore> {
			using CircleCore::CircleCore;
			DEF_IMODEL_FUNCS
			DEF_IMODEL_MASS
			bool isInner(const Vec2& pos) const override;
		};

		//! 円の剛体用クラス
		/*! キャッシュを管理する関係で形状へのアクセスは関数を通す仕様 元クラスはcompososition */
		class CircleModel : public IModelP<CircleCore>, public Cache<CircleCore> {
			public:
				CircleModel() = default;
				explicit CircleModel(const CircleCore& c);

				float getRadius() const;
				void setCenter(const Vec2& v);
				void setRadius(float r);

				DEF_IMODEL_FUNCS
				DEF_IMODEL_MASS
				bool isInner(const Vec2& pos) const override;
		};

		//! AxisAlignedBox
		struct BoxCore {
			DEF_CORE_FUNCS(TagArea, TagInertia, TagCenter, TagBCircle)
			// --------------------------------------------
			Vec2	minV, maxV;

			BoxCore() = default;
			BoxCore(const Vec2& min_v, const Vec2& max_v);

			Vec2 support(const Vec2& dir) const;
			Vec2 nearest(const Vec2& pos) const;
			bool hit(const LineCore& l) const;
		};
		struct Box : BoxCore, IModelP<BoxCore> {
			using BoxCore::BoxCore;
			DEF_IMODEL_FUNCS
			CircleCore getBCircle() const override;
		};
		class BoxModel : public IModelP<BoxCore>, public Cache<BoxCore> {
			BoxCore		_box;

			public:
				BoxModel() = default;
				BoxModel(const BoxCore& b);
				DEF_IMODEL_FUNCS
				DEF_IMODEL_MASS
		};
		struct PolyCore {
			DEF_CORE_FUNCS(TagArea, TagInertia, TagCenter, TagBCircle)
			// --------------------------------------------
			Vec2		point[3];

			PolyCore() = default;
			PolyCore(const Vec2& p0, const Vec2& p1, const Vec2& p2);

			bool isInTriangle(const Vec2& p) const;
			std::pair<Vec2,int> nearest(const Vec2& p) const;

			Vec2 support(const Vec2& dir) const;
			void addOffset(const Vec2& ofs);
			static float CalcArea(const Vec2& p0, const Vec2& p1, const Vec2& p2);
			static float CalcArea(const Vec2& p0, const Vec2& p1);
			//! 鈍角を探す
			/*! \return 鈍角の番号 (負数は該当なし) */
			int getObtuseCorner() const;
		};
		struct Poly : PolyCore, IModelP<PolyCore> {
			using PolyCore::PolyCore;
			DEF_IMODEL_FUNCS
			bool isInner(const Vec2& pos) const override;
			CircleCore getBCircle() const override;
		};
		class PolyModel : public IModelP<PolyCore>, public Cache<PolyCore>, public spn::CheckAlign<16,PolyModel> {
			public:
				PolyModel() = default;
				PolyModel(const Vec2& p0, const Vec2& p1, const Vec2& p2);
				PolyModel(const PolyModel& p) = default;

				void setPoint(int n, const Vec2& v);
				void addOffset(const Vec2& ofs);

				DEF_IMODEL_FUNCS
				DEF_IMODEL_MASS
				bool isInner(const Vec2& pos) const override;
		};

		//! 多角形基本クラス
		struct ConvexCore {
			using CT = CType<TagArea, TagInertia, TagCenter, TagBCircle>;
			template <class T>
			using Wrap = CoreWrap<T>;

			CircleCore bcircle() const;
			template <class INFO>
			void getInfo(INFO& info, TagArea) const {
				auto res = area_inertia_center();
				info.refCache(TagArea()) = std::get<0>(res);
				info.refCache(TagInertia()) = std::get<1>(res);
				info.refCache(TagCenter()) = std::get<2>(res);
			}
			template <class INFO>
			void getInfo(INFO& info, TagInertia) const {
				getInfo(info, TagArea()); }
			template <class INFO>
			void getInfo(INFO& info, TagCenter) const {
				getInfo(info, TagArea()); }

			// --------------------------------------------
			using AreaL = std::vector<float>;
			struct AreaSum {
				float result;

				AreaSum(): result(0) {}
				void operator()(int /*n*/, const Vec2& p0, const Vec2& p1) {
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
			ConvexCore(const ConvexCore& c) = default;
			/*! \param[in] v 凸包と分かっている頂点 */
			ConvexCore(std::initializer_list<Vec2> v);
			ConvexCore(const PointL& pl);
			ConvexCore(PointL&& pl);
			ConvexCore(ConvexCore&& c);
			ConvexCore& operator = (ConvexCore&& c);

			static ConvexCore FromConcave(const PointL& src);
			float area() const;

			/*! 同時に求めると少し効率が良い */
			std::tuple<float,float,Vec2> area_inertia_center() const;
			Vec2 support(const Vec2& dir) const;
			bool isInner(const Vec2& pos) const;
			//! 2つに分割
			/*! \param[out] c0 線分の進行方向左側
				\param[out] c1 線分の進行方向右側 */
			ConvexCore2 splitTwo(const StLineCore& l) const;
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

			//! 指定ポイントの内部的な領域IDと内外位置を取得
			/*! \return first=内外判定
						second=領域ID */
			IModel::CPos checkPosition(const Vec2& pos) const;
			//! 内部的な通し番号における外郭ライン
			LineCore getOuterLine(int n) const;
			StLineCore getOuterStLine(int n) const;
			IModel::PosL getOverlappingPoints(const IModel& mdl, const Vec2& inner) const;
			// 凸包が直線と交差している箇所を2点計算
			std::tuple<bool,Vec2,Vec2> checkCrossingLine(const StLineCore& l) const;
			// 直線がヒットしてるか判定後、始点と終点のチェック
			bool hit(const LineCore& l) const;
			ConvexCore& operator *= (const AMat32& m);

			std::ostream& dbgPrint(std::ostream& os) const;
		};
		struct Convex : ConvexCore, IModelP<ConvexCore> {
			using ConvexCore::ConvexCore;
			Convex() = default;
			Convex(ConvexCore&& c);
			DEF_IMODEL_FUNCS
			bool isInner(const Vec2& pos) const override;
			PosL getOverlappingPoints(const IModel &mdl, const Vec2& inner) const override;
			CircleCore getBCircle() const override;
			int getNPoints() const override;
			Vec2 getPoint(int n) const override;
			CPos checkPosition(const Vec2& pos) const override;
			Convex2 splitTwo(const StLineCore& l) const override;
			Convex& operator *= (const AMat32& m);
			std::ostream& dbgPrint(std::ostream& os) const override;
			//! 凸包が重なっている領域を求める
			/*! 既に重なっている事が分かっている前提
				\param[in] m0 モデルその1(Convexとみなす)
				\param[in] m1 モデルその2(Convexとみなす)
				\param[in] inner 重複領域の内部点 */
			static Convex GetOverlappingConvex(const IModel& m0, const IModel& m1, const Vec2& inner);
		};
		class ConvexModel : public IModelP<ConvexCore>, public Cache<ConvexCore>, public spn::CheckAlign<16,ConvexModel> {
			public:
				ConvexModel(const ConvexModel& m) = default;
				ConvexModel(std::initializer_list<Vec2> v);
				ConvexModel(const PointL& pl);
				ConvexModel(PointL&& pl);

				const PointL& getPoint() const;
				PointL& refPoint();

				void addOffset(const Vec2& ofs);
				DEF_IMODEL_FUNCS
				DEF_IMODEL_MASS
				bool isInner(const Vec2& pos) const override;
				PosL getOverlappingPoints(const IModel &mdl, const Vec2& inner) const override;
				int getNPoints() const override;
				Vec2 getPoint(int n) const override;
				CPos checkPosition(const Vec2& pos) const override;
				Convex2 splitTwo(const StLineCore& line) const override;
				std::ostream& dbgPrint(std::ostream& os) const override;
		};

		//! 剛体制御用
		class RPose : public spn::Pose2D {
			Vec2			_vel,
							_acc;
			GAP_MATRIX_DEF(mutable, _finalInv, 3,2,
			   ((float, _rotVel, float, _rotAcc))
			   ((uint32_t, _velAccum))			//!< Pose2Dのaccumカウンタとは別で速度に関する変数が書き換えられた際にインクリメント
			   ((mutable uint32_t, _invAccum))	//!< 逆行列キャッシュを作った時のカウンタ値
			)
			void _identitySingle();

			protected:
				void _setAsChanged();
			public:
				struct TValue;
				struct Value : spn::Pose2D::Value {
					Vec2	&vel, &acc;
					float	&rotVel, &rotAcc;

					Value(RPose& rp);
					~Value();
					Value(const Value& v);
					Value& operator = (const TValue& tv);
					void operator = (const Value& v) = delete;
				};
				struct TValue : spn::Pose2D::TValue {
					Vec2	vel, acc;
					float	rotVel, rotAcc;

					TValue() = default;
					TValue(const Value& v):
						vel(v.vel), acc(v.acc), rotVel(v.rotVel), rotAcc(v.rotAcc) {}
					TValue& operator = (const Value& v);
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

				const AMat32& getToLocal() const;
				RPose lerp(const RPose& p1, float t) const;
				Value refValue();
		};

		struct IResist;
		//! 座標変換ありのモデル基底
		/*! 既存のモデルに被せる形で使う */
		template <class MDL, class BASE>
		class TModel : public IModel, public BASE {
			MDL				_model;		//!< 形状

			const IModel& _getModel(const HLMdl& mdl) const;
			const IModel& _getModel(const IModel& mdl) const;

			public:
				TModel(const TModel& t) = default;
				TModel(const MDL& mdl);
				template <class... T>
				TModel(const MDL& mdl, const T&... args): BASE(args...), _model(mdl) {}
				const MDL& getModel() const;

				uint32_t getCID() const override;
				DEF_IMODEL_FUNCS
				DEF_IMODEL_MASS
				bool isInner(const Vec2& pos) const override;
				PosL getOverlappingPoints(const IModel& mdl, const Vec2& inner) const override;
				int getNPoints() const override;
				Vec2 getPoint(int n) const override;
				CPos checkPosition(const Vec2& pos) const override;
				Convex2 splitTwo(const StLineCore& line) const override;
				std::ostream& dbgPrint(std::ostream& os) const override;
		};

		class TR_Mat {
			AMat32	_mToLocal,
					_mToWorld;
			public:
				static struct tagInverse {} TagInverse;
				TR_Mat() = default;
				TR_Mat(const AMat32& m);
				TR_Mat(const TR_Mat& t) = default;
				TR_Mat(const TR_Mat& t, tagInverse);
				TR_Mat(const RPose& rp);
				TR_Mat(const RPose& rp, tagInverse);
				Vec2 toLocal(const Vec2& v) const;
				Vec2 toLocalDir(const Vec2& v) const;
				Vec2 toWorld(const Vec2& v) const;
				Vec2 toWorldDir(const Vec2& v) const;
				const AMat32& getToLocal() const;
				const AMat32& getToWorld() const;
		};
		template <class BASE>
		using TModelH = TModel<HLMdl, BASE>;
		template <class BASE>
		using TModelR = TModel<const IModel&, BASE>;

		//! narrow-phase collision check
		/*! IModelインタフェース同士をBCircleで判定後、GJKによる判定 */
		class NarrowC_Model {
			Vec2	_inner;

			public:
				using info_type = decltype(_inner);
				const Vec2& getInfo() { return _inner; }
				bool operator()(const IModel& mdl0, const IModel& mdl1);
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
