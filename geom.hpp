//! 2D, 3D形状に共通するクラスの定義
#pragma once
#include "spinner/type.hpp"
#include "spinner/misc.hpp"
#include "spinner/matrix.hpp"
#include "spinner/plane.hpp"
#include "spinner/optional.hpp"
#include "spinner/pose.hpp"
#include "spinner/interpolatation.hpp"
#include "spinner/check_macro.hpp"
#include <boost/variant.hpp>
#include "handle.hpp"

namespace boom {
	constexpr static float DOT_THRESHOLD = 1e-4f,		//!< 内積による表裏判定基準
							NEAR_THRESHOLD = 1e-4f,		//!< 同じ座標と判断する基準
							ZEROVEC_LENGTH = 1e-5f,
							SQUARE_RATIO = 1e-2f,
							DOT_THRESHOLD_SQ = DOT_THRESHOLD * SQUARE_RATIO,
							NEAR_THRESHOLD_SQ = NEAR_THRESHOLD * SQUARE_RATIO,
							ZEROVEC_LENGTH_SQ = ZEROVEC_LENGTH * SQUARE_RATIO;

	using spn::Vec2; using spn::Vec3; using spn::Vec4;
	using spn::AVec2; using spn::AVec3; using spn::AVec4;
	using spn::Mat22; using spn::Mat32; using spn::Mat33; using spn::Mat43; using spn::Mat44;
	using spn::AMat22; using spn::AMat32; using spn::AMat33; using spn::AMat43; using spn::AMat44; using spn::Plane;
	using spn::Quat; using spn::AQuat;
	using spn::Plane; using spn::APlane;
	using spn::Pose2D; using spn::Pose3D;
	using spn::Bit;
	using Float2 = std::pair<float,float>;
	using Int2 = std::pair<int,int>;
	using Int2x2 = std::pair<Int2,Int2>;
	using Uint2 = std::pair<uint32_t, uint32_t>;
	using Vec2x2 = std::pair<Vec2, Vec2>;
	using Vec3x2 = std::pair<Vec3, Vec3>;
	using Vec3List = std::vector<Vec3>;
	using PlaneList = std::vector<Plane>;
	using IndexList = std::vector<uint16_t>;

	//! 90度回転行列(2D)
	extern const AMat22 cs_mRot90[2];
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
		return spn::CramersRule(toV1, toV2, toVT, spn::Rcp22Bit(det));
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
		float r = line2.length() * spn::Rcp22Bit(line.length());
		return spn::Lerp(f0, f1, r);
	}

	//! サポートされていない関数を読んだ時の実行時エラーを送出
	#define INVOKE_ERROR Assert(Trap, false, "not supported function: ", __func__) throw 0;
	//! IModelから指定の型にキャスト出来ればその参照を返す関数のデフォルト実装
	#define DEF_CASTFUNC(typ) virtual spn::Optional<typ&> castAs##typ() { return spn::none; } \
			virtual spn::Optional<const typ&> castAs##typ() const { auto ref = const_cast<IModel*>(this)->castAs##typ(); \
			if(ref) return *ref; return spn::none; }
	//! コリジョンID取得関数を付加
	/*! \tparam T	コリジョン構造体
		\tparam CT	コリジョン型リスト(CType) */
	template <class T, class CT>
	struct ITagP_base {
		static constexpr uint32_t GetCID() { return CT::template Find<T>::result; }
	};
	/*! \tparam T	コリジョン構造体 (ITagP_base)
		\tparam MDL	Modelインタフェース */
	template <class T, class MDL>
	struct IModelP_base : virtual MDL, T {
		using T::T;
		IModelP_base() = default;
		IModelP_base(T&& t): T(std::move(t)) {}
		IModelP_base(const T& t): T(t) {}
		const void* getCore() const override { return static_cast<const T*>(this); }
		uint32_t getCID() const override { return T::GetCID(); }
	};
	struct IModelNode {
		using Time_t = uint32_t;

		virtual ~IModelNode() {}
		virtual bool imn_refresh(Time_t t) const;
		virtual bool hasInner() const;
		//! モデルの実体へのポインタ
		virtual const void* getCore() const;
		//! 最寄りのユーザーデータを取得
		/*! このノードが持っていればそれを返し、無ければ親を遡って探す */
		virtual void* getUserData();
		virtual std::ostream& print(std::ostream& os) const = 0;
		friend std::ostream& operator << (std::ostream& os, const IModelNode& mdl);
	};

	//! IModelとHMdlの差異吸収
	template <class MMGR>
	class VModel {
		using HMdl = typename MMGR::SHdl;
		using HLMdl = typename MMGR::LHdl;
		using IModel = typename MMGR::data_type::element_type;
		using VMdl = boost::variant<const IModel&, HLMdl>;
		VMdl	_vMdl;

		struct Visitor : boost::static_visitor<const IModel&> {
			const IModel& operator()(const HLMdl& hl) const {
				return *hl.cref().get(); }
			const IModel& operator()(const IModel& m) const { return m; }
		};
		public:
			VModel(const IModel& mdl): _vMdl(mdl) {}
			VModel(HMdl hMdl): _vMdl(hMdl) {}
			const IModel& get() const {
				return boost::apply_visitor(Visitor(), _vMdl);
			}
	};

	//! ラインの位置を示す
	enum class LinePos {
		Begin,		//!< 始点
		End,		//!< 終点
		OnLine,		//!< ライン上
		NotHit
	};
	//! 凸包の位置を示す
	enum class ConvexPos {
		Inner,
		OnLine,
		Outer
	};
	//! 凸包の直線との位置関係
	enum LineDivision {
		OnLine = 0x00,
		Cw = 0x01,
		Ccw = 0x02,
		Bridge = 0x03
	};

	template <class... Ts>
	struct CCType : spn::CType<Ts...> {
		template <class CORE>
		struct CAnother {
			using result = CCType<decltype(Ts::template get<CORE>(std::declval<CORE>()))...>;
		};
	};

	#define DEF_CACHETAG(name, valueT, method)	struct name { \
		valueT value; \
		valueT operator ()() const; \
		template <class T> static auto get(const T& t) -> decltype(std::declval<T>().method()) { return t.method(); } };
	template <class T>
	struct GetSize_Ref0 {
		static constexpr int get() {
			return (std::is_reference<T>::value == 0) ? sizeof(T) : 0;
		}
	};
	template <class CT, class T>
	struct FlagGet {
		static constexpr int get() {
			return 1 << CT::template Find<T>::result; }
	};

	#define DEF_SHARE_F(z, data, elem)	decltype(std::declval<BOOST_PP_SEQ_ELEM(0,data)>().BOOST_PP_SEQ_ELEM(1,data)()) CalcIt(elem) const { return BOOST_PP_SEQ_ELEM(1,data)(); }
	#define DEF_SHARE(obj, method, seq)	BOOST_PP_SEQ_FOR_EACH(DEF_SHARE_F, (obj)(method), seq)
	template <class T>
	struct CacheBase : T {
		template <class TAG>
		auto CalcIt(TAG) const -> std::tuple<TAG> { return std::tuple<TAG>(TAG{TAG::get((const T&)*this)}); }
	};

	template <int N>
	struct Setter {
		template <class CACHE, class... Ts>
		static void set(CACHE& c, const std::tuple<Ts...>& tup) {
			using Typ = typename std::decay<decltype(std::get<N-1>(tup))>::type;
			c.template _getBuff<Typ>() = std::get<N-1>(tup).value;
			Setter<N-1>::set(c, tup);
		}
	};
	template <>
	struct Setter<0> {
		template <class CACHE, class... Ts>
		static void set(CACHE& /*c*/, const std::tuple<Ts...>& /*tup*/) {}
	};

	#define DEF_GETMETHOD(clazz, name, tag)	decltype(std::declval<clazz>().template _getInfo<tag>()) name() const { return clazz::template _getInfo<tag>(); }
	//! キャッシュ管理自動化クラス
	/*! \param CTTAG キャッシュタグリスト
		\param CORE 形状を定義する構造体 */
	template <class CTTAG, class CORE>
	class Cache {
		private:
			#define DEF_SETTER_FRIEND(z,n,dummy)	friend struct Setter<n>;
			BOOST_PP_REPEAT(8, DEF_SETTER_FRIEND, NOTHING)
			#undef DEF_SETTER_FRIEND

			template <class T>
			using T_OP = decltype(std::declval<T>()());

			//! キャッシュタグからフラグを算出
			template <class T>
			using FlagGet0 = FlagGet<CTTAG, T>;
			//! COREから実際に返される型
			using CTRET = typename CTTAG::template CAnother<CORE>::result;
			//! キャッシュ変数を格納するメモリ領域
			mutable uint8_t		_buff[CTRET::template SumN<CTRET::size, GetSize_Ref0>::result];
			//! キャッシュ管理フラグ (ビットが1なら有効)
			mutable uint32_t	_cacheFlag = 0;

			//! tupleで返された値を保存し、対応するフラグを立てる
			template <class... Ts>
			void _setInfo(const std::tuple<Ts...>& tup) const {
				Setter<sizeof...(Ts)>::set(*this, tup);
			}
			template <class T>
			T_OP<T>& _getBuff() const {
				uint8_t* ptr = (uint8_t*)_buff + CTTAG::template SumT<T>::result;
				return *reinterpret_cast<T_OP<T>*>(ptr);
			}
			//! キャッシュ対応の値
			template <class T>
			const T_OP<T>& _getInfo(std::false_type) const {
				using spn::Bit;
				// フラグが立って無ければ値を算出
				constexpr int flag = FlagGet0<T>::get();
				if(!Bit::Check(_cacheFlag, flag)) {
					Bit::Set(_cacheFlag, flag);
					_setInfo(_core.CalcIt(T()));
				}
				return _getBuff<T>();
			}
			//! キャッシュ非対応の値
			template <class T>
			const T_OP<T>& _getInfo(std::true_type) const {
				// そのまま返す
				_core.CalcIt(T());
			}

		protected:
			//! 元となる構造体
			CORE 		_core;
			//! 特定のフラグをリセットする
			template <class... Tags>
			void _invalidate() {
				auto val = spn::CType<Tags...>::template SumN<sizeof...(Tags), FlagGet0>::result;
				spn::Bit::Clear(_cacheFlag, val);
			}
		public:
			template <class T>
			const T_OP<T>& _getInfo() const {
// 				constexpr int idx = CTTAG::template Find<T>::result;
				return _getInfo<T>(std::is_reference<decltype(T::get(std::declval<CORE>()))>());
			}
			const CORE& getCoreRef() const { return _core; }
			CORE& getCoreRef() { return _core; }
	};
	DEF_CHECKMETHOD_OV(mdlhit, hit)
	using ColFunc = bool (*)(const void*, const void*);

	template <class Types>
	struct Narrow_Init0 {
		using CTGeo = typename Types::CTGeo;
		using GJK = typename Types::GJK;
		using IModel = typename Types::IModel;

		template <class T0, class T1>
		static bool Invoke(std::true_type, const IModel* t0, const IModel* t1) {
			// 専用アルゴリズムでの当たり判定
			return reinterpret_cast<const T0*>(t0->getCore())->hit(*reinterpret_cast<const T1*>(t1->getCore()));
		}
		template <class T0, class T1>
		static bool Invoke(std::false_type, const IModel* t0, const IModel* t1) {
			//GJKアルゴリズムでの当たり判定
			GJK gjk(*t0,*t1);
			return gjk.getResult();
		}
		template <class T0, class T1>
		static bool _CFRaw(const IModel* t0, const IModel* t1, std::true_type) {
			return Invoke<T0,T1>(CheckMethod_mdlhit<T0, const T1&>(), t0, t1); }
		template <class T0, class T1>
		static bool _CFRaw(const IModel* t0, const IModel* t1, std::false_type) {
			return Invoke<T1,T0>(CheckMethod_mdlhit<T1, const T0&>(), t1, t0); }
		//! リスト登録用関数ラッパ
		template <int N0, int N1, class B0, class B1,
					class = typename std::enable_if<spn::TType<B0,B1>::t_and::value>::type>
		static bool CFRaw(const void* m0, const void* m1) {
			using T0 = typename CTGeo::template At<N0>::type;
			using T1 = typename CTGeo::template At<N1>::type;
			return _CFRaw<T0,T1>(reinterpret_cast<const IModel*>(m0),
								 reinterpret_cast<const IModel*>(m1),
								 typename spn::NType<T0::GetCID(), T1::GetCID()>::less_eq());
		}
		template <int N0, int N1, class B0, class B1,
					class = int,
					class = typename std::enable_if<spn::TType<B0,B1>::t_nand::value>::type>
		static bool CFRaw(const void* /*m0*/, const void* /*m1*/) {
			return false;
		}
	};
	template <class Types, int M, int N>
	struct Narrow_InitB {
		static void Init(ColFunc*& dst) {
			using CTGeo = typename Types::CTGeo;
			*dst-- = Narrow_Init0<Types>::template CFRaw<M,N,
				typename spn::NType<M, CTGeo::size>::less, typename spn::NType<N, CTGeo::size>::less>;
			Narrow_InitB<Types, M, N-1>::Init(dst);
		}
	};
	template <class Types, int M>
	struct Narrow_InitB<Types, M, -1> {
		static void Init(ColFunc*& /*dst*/) {}
	};

	template <class Types, int WIDE_M, int M>
	struct Narrow_InitA {
		static void Init(ColFunc*& dst) {
			Narrow_InitB<Types, M, WIDE_M-1>::Init(dst);
			Narrow_InitA<Types, WIDE_M, M-1>::Init(dst);
		}
	};
	template <class Types, int WIDE_M>
	struct Narrow_InitA<Types, WIDE_M, -1> {
		static void Init(ColFunc*& /*dst*/) {}
	};

	template <class T, class IM>
	struct IModelSP : IM {
		const T& source;
		IModelSP(const T& src): source(src) {}
		using VEC = decltype(std::declval<IM>().im_getCenter());
		VEC im_support(const VEC& v) const override {
			return source.support(v);
		}
		const void* getCore() const override {
			return static_cast<const T*>(&source); }
	};
	//! Narrow-phase 判定関数群
	template <class Types>
	struct Narrow {
		using CTGeo = typename Types::CTGeo;
		using IModel = typename Types::IModel;
		using GJK = typename Types::GJK;
		using Time_t = typename IModelNode::Time_t;

		constexpr static int WideBits = spn::CTBit::MSB_N<CTGeo::size>::result + 1,
							ArraySize = 1<<(WideBits*2);
		static ColFunc cs_cfunc[ArraySize];

		//! 当たり判定を行う関数をリストにセットする
		static void Initialize() {
			constexpr int WideM = (1<<WideBits);
			ColFunc* cfp = cs_cfunc+ArraySize-1;
			Narrow_InitA<Types, WideM, WideM-1>::Init(cfp);
		}
		//! 当たり判定を行う関数ポインタを取得
		static ColFunc GetCFunc(int id0, int id1) {
			return cs_cfunc[(id0 << WideBits) | id1];
		}
		//! IModelインタフェースを介さず直接判定
		template <class T0, class T1,
			class = typename std::enable_if<!std::is_pointer<T0>::value>::type,
			class = typename std::enable_if<!std::is_pointer<T1>::value>::type>
		static bool Hit(const T0& t0, const T1& t1) {
			constexpr int id0 = CTGeo::template Find<T0>::result,
						id1 = CTGeo::template Find<T1>::result;
			IModelSP<T0, IModel> tmp0(t0);
			IModelSP<T1, IModel> tmp1(t1);
			return GetCFunc(id0, id1)(&tmp0, &tmp1);
		}
		//! 2つの物体(階層構造可)を当たり判定
		static bool Hit(const IModel* mdl0, const IModel* mdl1, Time_t t) {
			if(mdl0->imn_refresh(t) && mdl1->imn_refresh(t)) {
				if(GetCFunc(mdl0->getCID(), mdl1->getCID())(mdl0, mdl1)) {
					if(mdl0->hasInner() | mdl1->hasInner())
						return _HitL(mdl1, mdl0, false, t) ||
								_HitL(mdl0, mdl1, true, t);
					return true;
				}
			}
			return false;
		}
		//! mdl0を展開したものとmdl1を当たり判定
		/*! mdl1は既にrefreshをかけてある前提
			\param[in] mdl0 展開する方のインタフェース
			\param[in] mdl1 展開されない方のインタフェース
			\param[in] bSwap 左右を入れ替えて判定している場合はtrue */
		static bool _HitL(const IModel* mdl0, const IModel* mdl1, bool bSwap, Time_t t) {
			if(!mdl0->imn_refresh(t))
				return false;
			auto in = mdl0->getInner();
			if(in) {
				// TODO: ここで判定関数を取得してInnerに含まれる子オブジェクトは全て同じ型 という前提にする
				do {
					ColFunc cf = GetCFunc(in.get()->getCID(), mdl1->getCID());
					if(cf(in.get(), mdl1)) {
						if(_HitL(mdl1, in.get(), false, t))
							return true;
					}
				} while(++in);
			} else {
				if(bSwap)
					return true;
				return _HitL(mdl1, mdl0, true, t);
			}
			return false;
		}
	};
	template <class Types>
	ColFunc Narrow<Types>::cs_cfunc[ArraySize];
}
