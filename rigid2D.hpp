#pragma once
#include "geom2D.hpp"
#include "spinner/withbool.hpp"

namespace boom {
	namespace geo2d {
		//! ペナルティ法における抗力と摩擦力
		struct RForce {
			struct F {
				Vec2	linear;
				float	torque;

				F() = default;
				F(float s): linear(s), torque(s) {}
				F& operator += (const F& f);
				F& operator *= (float s);
				F operator * (float s) const;
				friend std::ostream& operator << (std::ostream& os, const F& f);
			};
			F	sdump,	//!< スプリング & ダンパによる力
			fricD;	//!< 動摩擦力

			RForce() = default;
			RForce(float s): sdump(s), fricD(s) {}
			RForce& operator += (const RForce& rf);
			RForce& operator *= (float s);
			RForce operator * (float s) const;
			friend std::ostream& operator << (std::ostream& os, const RForce& f);
		};
		using CResult = ColResult<512, c_info::Pairs<c_ent::Sum<RForce>, uint32_t>>;
		//! 剛体ラッパ (形状 + 姿勢)
		class Rigid : public TModelH<RPose>, public spn::CheckAlign<16,Rigid> {
			public:
				constexpr static int NUM_RESIST = 4;
				using SPResist = std::shared_ptr<IResist>;
			private:
				SPResist	_resist[NUM_RESIST];	//!< 抵抗計算用
				uint32_t	_rID[NUM_RESIST];
			public:
				constexpr static uint32_t DEFAULT_ID = ~0;
				using TModelH<RPose>::TModelH;

				// --- シミュレーションに関する関数など ---
				RPose& refPose();
				const RPose& getPose() const;

				using CheckAlign<16,Rigid>::NewUF;
				void addR(const SPResist& sp, uint32_t id=DEFAULT_ID);
				//! 登録時に付加したIDのresistインタフェースを取得
				const SPResist& getR(uint32_t id);
				void remR(uint32_t id);
				void remRAll();
				//! 抵抗力を計算
				/*! \param[in] index 通し番号 */
				RForce::F resist(int index, const CResult& cr) const;
				//! 質点に対して力を加える
				void addForce(const Vec2& wpos, const Vec2& f);
				void addLinearForce(const Vec2& f);
				void applyForce(const RForce::F& f);
				static RForce::F CalcForce(const RPose& rp, const Vec2& wpos, const Vec2& f);
		};
		using UPRigid = decltype(Rigid::NewUF());

		#define mgr_rigid reinterpret_cast<RigidRes&>(ModelMgr::_ref())
		class RigidRes : public ModelMgr {
			public:
				// unique_ptrなのでポインタのコピーは不可。常にmoveで指定
				LHdl acquireModel(UPModel&& mdl) {
					return acquire(std::move(mdl));
				}
				AnotherLHandle<UPRigid> acquireRigid(UPRigid&& rig) {
					LHdl lh = acquire(std::move(rig));
					return Cast<UPRigid>(std::move(lh));
				}
		};
		DEF_HANDLE(RigidRes, Rig, UPRigid)

		struct IItg;
		//! 剛体マネージャ
		class RigidMgr {
			using SPItg = std::shared_ptr<IItg>;
			using NPID = uint16_t;
			struct NPair {
				Vec2	pos,				//!< 接触位置
						dir;				//!< 接触法線
				HRig	hAnother;			//!< 相手のハンドル
				int		nFrame;				//!< 連続ヒット数 (接触判定には関係しないが外部から物体の静止判定をする際に使用)
				NPID	nextID, prevID;		/*!< ~0の時は無効 */
			};
			using NP_Pool = spn::noseq_list<NPair, NPID>;
			NP_Pool			_npPool;

			struct ERig {
				using TPose = typename spn::Pose2D::TValue;

				HRig		hRig;
				TPose		prePose;		//!< 前の姿勢
				uint16_t	nHitPrev,		//!< 1フレーム前のHitエントリカウント
							nHitCur,		//!< 現在のHitエントリカウント
							nHitNew;		//!< 新しく追加したHitエントリの数
				NPID		hitID;			//!< Hitリストの先頭

				constexpr static NPID invalid = ~0;

				ERig(HRig hR);
				operator HRig () const;
				operator const Rigid& () const;
				operator Rigid& ();
				void savePose();

				void resetHitCount();
				//! エントリ数を比較してもう衝突していないエントリを削除
				void removeOld(NP_Pool& pool);
				//! 既存のエントリを探して先頭に持ってくる。エントリが無ければ作成
				spn::WithBool<NPair&> addHit(HRig hR, NP_Pool& pool);
			};
			using BroadC = BroadC_RoundRobin<ERig, const Rigid&>;
			using NarrowC = NarrowC_Model;
			BroadC			_broadC;
			SPItg			_itg;		//!< 適用する積分アルゴリズム
			void _checkCollision();		//! コリジョンチェックして内部変数に格納
			// NOTE: ひとまずは反発係数固定での実装
			RCoeff			_coeff = {0,0,0,0,0};	//!< 力積計算に使う係数
			CResult			_cresult;				//!< 力積をRigidから参照する為のクラス

			public:
				using iterator = typename BroadC::iterator;
				using id_type = typename BroadC::id_type;

				RigidMgr(const SPItg& itg);
				void setCoeff(const RCoeff& coeff);
				void simulate(float dt);
				//! broad/narrow phase collisionに通った物体のペアを受け取る関数。Inner-pointを受け取り力積計算をする
				/*!	BroadCが渡すオブジェクトに適合していればIModelだろうがなんだろうがOK
				 *					\param[in] id0 ObjectAの通し番号
				 *					\param[in] id1 ObjectBの通し番号
				 *					\param[in] mdl0 ObjectAのインタフェース
				 *					\param[in] mdl1 ObjectBのインタフェース */
				void operator ()(int id0, int id1, ERig& er0, ERig& er1, const Vec2& inner);

				id_type addA(HRig hRig);
				id_type addB(HRig hRig);
				void remA(id_type id);
				void remB(id_type id);
				typename BroadC::const_iterator cbeginA() const;
				typename BroadC::const_iterator cendA() const;
		};

		//! 位置と速度を与えた時にかかる加速度(抵抗力)を計算
		struct IResist : std::enable_shared_from_this<IResist> {
			using sptr = std::shared_ptr<IResist>;
			using csptr = const sptr&;

			//! 抵抗力を計算
			/*! \param[inout] acc 加速度
			*	\param[in] r 姿勢
			*	\param[in] index 自分のID
			*	\param[in] cr 当たり判定結果 */
			virtual void resist(RForce::F& acc, const Rigid& r, int index, const CResult& cr) const = 0;
		};

		using RItr = RigidMgr::iterator;
		using TValueA = std::vector<RPose::TValue>;
		//! 積分インタフェース
		struct IItg {
			using sptr = std::shared_ptr<IItg>;
			using csptr = const sptr&;

			//! 剛体の姿勢を積分計算
			/*! \param[inout] st 現在の姿勢
			*	\param[in] dt 前回フレームからの時間(sec) */
			virtual void advance(int pass, int offset, RItr itr, RItr itrE, const CResult& cr, float dt) = 0;
			//! 1回分の時間を進めるのに必要なステップ数
			virtual int numOfIteration() const = 0;

			virtual void beginIteration(int /*n*/) {}
			virtual void endIteration() {}
		};

		// ---------------------- 積分アルゴリズム ----------------------
		namespace itg {
			#define DEF_IITG_FUNCS \
			virtual void advance(int pass, int offset, RItr itr, RItr itrE, const CResult& cr, float dt) override; \
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
		// ---------------------- 抵抗力計算 ----------------------
		namespace resist {
			//! 衝突による加速度
			class Impact : public IResist {
				public:
					void resist(RForce::F& acc, const Rigid& r, int index, const CResult& cr) const override;
			};
			//! 重力による加速度
			class Gravity : public IResist {
				Vec2	_grav = Vec2(0,0);
				public:
					Gravity() = default;
					Gravity(const Vec2& v);
					void setGravity(const Vec2& v);
					void resist(RForce::F& acc, const Rigid& r, int index, const CResult& cr) const override;
			};
			//! 空気抵抗
			/*! 速度に比例した抵抗を掛ける */
			class Air : public IResist {
				float	_cLinear=0, _cRot=0;
				public:
					Air() = default;
					Air(float cLinear, float cRot);
					void setAir(float cLinear, float cRot);
					void resist(RForce::F& acc, const Rigid& r, int index, const CResult& cr) const override;
			};
			//! 点ジョイント
			class Point : public IResist {
				Vec2	_pos;
				public:
					Point() = default;
					Point(const Vec2& v);
					void setPos();
					void resist(RForce::F& acc, const Rigid& r, int index, const CResult& cr) const override;
			};
		}
		using DualRForce = DualValue<RForce>;
		//! ペナルティ法における抗力計算
		/*! 重なり領域だけを使うことは無くて、抗力の算出とセットなので内部で積分計算 */
		DualRForce CalcForce(const Rigid& r0, const Rigid& r1, const Vec2& inner, const RCoeff& coeff, const StLineCore& div);
	}
}
