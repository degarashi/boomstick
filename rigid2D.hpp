#pragma once
#include "geom2D.hpp"

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
		class Rigid;
		/*! NarrowC_ModelからInnerを受け取り, 力積計算 */
		class RigidCR : public ColResult<512, SharedEntry<AverageEntry<RForce>, uint32_t>> {
			// NOTE: ひとまずは反発係数固定での実装
			RCoeff		_coeff;		//!< 力積計算に使う係数
			public:
				using narrow_type = NarrowC_Model;
				RigidCR(const RCoeff& c);
				//! 当たり判定チェック(narrow-phase)をしてヒットしていたら力積計算
				/*! \param[in] id0 ObjectAの通し番号
				*					\param[in] id1 ObjectBの通し番号
				*					\param[in] mdl0 ObjectAのインタフェース
				*					\param[in] mdl1 ObjectBのインタフェース */
				void operator ()(int id0, int id1, const Rigid& r0, const Rigid& r1, const Vec2& inner);
				// BroadCが渡すオブジェクトに適合していればIModelだろうがなんだろうがOK
		};

		struct IItg;
		//! 剛体マネージャ
		class RigidMgr : public RigidCR {
			using SPRigid = std::shared_ptr<Rigid>;
			using SPItg = std::shared_ptr<IItg>;
			using BroadC = BroadC_RoundRobin<Rigid>;
			BroadC			_broadC;
			SPItg			_itg;		//!< 適用する積分アルゴリズム
			void _checkCollision();		//! コリジョンチェックして内部変数に格納

			public:
				using const_iterator = typename BroadC::const_iterator;
				using id_type = typename BroadC::id_type;

				RigidMgr(const SPItg& itg, const RCoeff& coeff);
				void simulate(float dt);

				id_type addA(const SPRigid& sp);
				id_type addB(const SPRigid& sp);
				void remA(id_type id);
				void remB(id_type id);
		};
		//! 剛体ラッパ (形状 + 姿勢)
		class Rigid : public TModelSP<RPose>, public spn::CheckAlign<16,Rigid> {
			public:
				constexpr static int NUM_RESIST = 4;
				using sptr = std::shared_ptr<Rigid>;
				using csptr = const sptr&;
				using SPResist = std::shared_ptr<IResist>;
			private:
				SPResist		_resist[NUM_RESIST];	//!< 抵抗計算用
			public:
				using TModelSP<RPose>::TModelSP;

				// --- シミュレーションに関する関数など ---
				RPose& refPose();
				const RPose& getPose() const;

				using CheckAlign<16,Rigid>::New;
				void addR(const SPResist& sp);
				//! 抵抗力を計算
				/*! \param[in] index 通し番号 */
				RForce::F resist(int index, const RigidCR& cr) const;
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
			virtual void resist(RForce::F& acc, const Rigid& r, int index, const RigidCR& cr) const = 0;
		};

		using RItr = RigidMgr::const_iterator;
		using TValueA = std::vector<RPose::TValue>;
		//! 積分インタフェース
		struct IItg {
			using sptr = std::shared_ptr<IItg>;
			using csptr = const sptr&;

			//! 剛体の姿勢を積分計算
			/*! \param[inout] st 現在の姿勢
			*	\param[in] dt 前回フレームからの時間(sec) */
			virtual void advance(int pass, RItr itr, RItr itrE, const RigidCR& cr, float dt) = 0;
			//! 1回分の時間を進めるのに必要なステップ数
			virtual int numOfIteration() const = 0;

			virtual void beginIteration(int n) {}
			virtual void endIteration() {}
		};

		// ---------------------- 積分アルゴリズム ----------------------
		namespace itg {
			#define DEF_IITG_FUNCS \
			virtual void advance(int pass, RItr itr, RItr itrE, const RigidCR& cr, float dt) override; \
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
					void resist(RForce::F& acc, const Rigid& r, int index, const RigidCR& cr) const override;
			};
			//! 重力による加速度
			class Gravity : public IResist {
				Vec2	_grav;
				public:
					Gravity(const Vec2& v);
					void resist(RForce::F& acc, const Rigid& r, int index, const RigidCR& cr) const override;
			};
		}

		//! ペナルティ法における抗力計算
		/*! 重なり領域だけを使うことは無くて、抗力の算出とセットなので内部で積分計算 */
		RForce CalcForce(const Rigid& r0, const Rigid& r1, const Vec2& inner, const RCoeff& coeff, const StLineCore& div);
	}
}
