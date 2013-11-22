//! 2D, 3D形状に共通するクラスの定義
#pragma once
#include "spinner/type.hpp"
#include "spinner/misc.hpp"
#include "spinner/matrix.hpp"
#include "spinner/plane.hpp"
#include "spinner/optional.hpp"

namespace boom {
	using spn::Vec2;
	using spn::Vec3;
	using spn::Mat33;
	using spn::Mat43;
	using spn::Mat44;
	using spn::AMat44;
	using spn::AMat43;
	using spn::Plane;
	using Vec2x2 = std::pair<Vec2, Vec2>;
	using Vec3x2 = std::pair<Vec3, Vec3>;
	//! サポートされていない関数を読んだ時の実行時エラーを送出
	#define INVOKE_ERROR Assert(Trap, false, "not supported function: ", __func__) throw 0;
	//! IModelから指定の型にキャスト出来ればその参照を返す関数のデフォルト実装
	#define DEF_CASTFUNC(typ) virtual spn::Optional<typ&> castAs##typ() { return spn::none; } \
			virtual spn::Optional<const typ&> castAs##typ() const { auto ref = const_cast<IModel*>(this)->castAs##typ(); \
			if(ref) return *ref; return spn::none; }
	template <class T, class CT>
	struct ITagP_base {
		static constexpr uint32_t GetCID() { return CT::template Find<T>::result; }
	};
	template <class T, class MDL>
	struct IModelP_base : MDL, T {
		using T::T;
		const void* getCore() const override { return static_cast<const T*>(this); }
		uint32_t getCID() const override { return T::GetCID(); }
	};
	//! IModelインタフェースの子ノードイテレータ
	template <class T>
	struct PtrItr {
		const uint8_t*	ptr;
		const size_t	stride;

		PtrItr(): ptr(nullptr), stride(0) {}
		PtrItr(const void* p, size_t st):
			ptr(reinterpret_cast<const uint8_t*>(p)), stride(st) {}
		PtrItr& operator ++ () {
			ptr += stride;
			return *this;
		}
		bool operator == (const PtrItr& pi) const {
			return ptr == pi.ptr; }
		bool operator != (const PtrItr& pi) const {
			return !(operator == (pi)); }
		const T& operator * () const {
			return *(operator -> ()); }
		const T* operator -> () const {
			return get(); }
		const T* get() const {
			return reinterpret_cast<const T*>(ptr); }
		int operator - (const PtrItr& p) const {
			auto diff = ptr - p.ptr;
			return diff / stride;
		}
	};
	//! ラインの位置を示す
	enum class LINEPOS {
		Begin,		//!< 始点
		End,		//!< 終点
		OnLine,		//!< ライン上
		NotHit
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
		static void set(CACHE& c, const std::tuple<Ts...>& tup) {}
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
				constexpr int idx = CTTAG::template Find<T>::result;
				return _getInfo<T>(std::is_reference<decltype(T::get(std::declval<CORE>()))>());
			}
			const CORE& getCoreRef() const { return _core; }
			CORE& getCoreRef() { return _core; }
	};
	DEF_CHECKMETHOD_OV(mdlhit, hit)
	using ColFunc = bool (*)(const void*, const void*);

	template <class CTG>
	struct Narrow_Init0 {
		template <class T0, class T1>
		static bool Invoke(std::true_type, const T0& t0, const T1& t1) {
			// 専用アルゴリズムでの当たり判定
			return t0.hit(t1); }
		template <class T0, class T1>
		static bool Invoke(std::false_type, const T0& t0, const T1& t1) {
			//TODO: GJKアルゴリズムでの当たり判定
			Assert(Trap, false, "not implemented yet (GJK algorithm)") throw 0; }
		template <class T0, class T1>
		static bool _CFRaw(const T0& t0, const T1& t1, std::true_type) {
			return Invoke(CheckMethod_mdlhit<T0, const T1&>(), t0, t1); }
		template <class T0, class T1>
		static bool _CFRaw(const T0& t0, const T1& t1, std::false_type) {
			return Invoke(CheckMethod_mdlhit<T1, const T0&>(), t1, t0); }
		//! リスト登録用関数ラッパ
		template <int N0, int N1, class B0, class B1,
					class = typename std::enable_if<spn::TType<B0,B1>::t_and::value>::type>
		static bool CFRaw(const void* m0, const void* m1) {
			using T0 = typename CTG::template At<N0>::type;
			using T1 = typename CTG::template At<N1>::type;
			return _CFRaw(*reinterpret_cast<const T0*>(m0),
						  *reinterpret_cast<const T1*>(m1),
						  typename spn::NType<T0::GetCID(), T1::GetCID()>::less_eq());
		}
		template <int N0, int N1, class B0, class B1,
					class = int,
					class = typename std::enable_if<spn::TType<B0,B1>::t_nand::value>::type>
		static bool CFRaw(const void* m0, const void* m1) {
			return false;
		}
	};
	template <class CTG, int WIDE, int N>
	struct Narrow_Init {
		static void Init(ColFunc* dst) {
			constexpr int Left = N >> WIDE,
						Right = N&((1<<WIDE)-1);
			*dst = Narrow_Init0<CTG>::template CFRaw<Left, Right,
				typename spn::NType<Left, CTG::size>::less, typename spn::NType<Right, CTG::size>::less>;
			Narrow_Init<CTG, WIDE, N-1>::Init(dst-1);
		}
	};
	template <class CTG, int WIDE>
	struct Narrow_Init<CTG, WIDE, -1> {
		static void Init(ColFunc* dst) {}
	};

	//! Narrow-phase 判定関数群
	template <class CTG, class IM>
	struct Narrow {
		constexpr static int WideBits = spn::CTBit::MSB_N<CTG::size>::result + 1,
							ArraySize = 1<<(WideBits*2);
		static ColFunc cs_cfunc[ArraySize];

		//! 当たり判定を行う関数をリストにセットする
		static void Initialize() {
			Narrow_Init<CTG, WideBits, ArraySize-1>::Init(cs_cfunc+ArraySize-1);
		}
		//! 当たり判定を行う関数ポインタを取得
		static ColFunc GetCFunc(int id0, int id1) {
			return cs_cfunc[(id0 << WideBits) | id1];
		}
		template <class T0, class T1,
			class = typename std::enable_if<!std::is_pointer<T0>::value>::type,
			class = typename std::enable_if<!std::is_pointer<T1>::value>::type>
		static bool Hit(const T0& t0, const T1& t1) {
			constexpr int id0 = CTG::template Find<T0>::result,
						id1 = CTG::template Find<T1>::result;
			return GetCFunc(id0, id1)(&t0, &t1);
		}
		//! 2つの物体(階層構造可)を当たり判定
		static bool Hit(const IM* mdl0, const IM* mdl1) {
			if(GetCFunc(mdl0->getCID(), mdl1->getCID())(mdl0->getCore(), mdl1->getCore()))
				return HitL(mdl1, mdl0, false);
			return false;
		}
		//! mdl0を展開したものとmdl1を当たり判定
		/*! \param[in] mdl0 展開する方のインタフェース
			\param[in] mdl1 展開されない方のインタフェース
			\param[in] bSwap 左右を入れ替えて判定している場合はtrue */
		static bool HitL(const IM* mdl0, const IM* mdl1, bool bSwap) {
			auto in = mdl0->getInner();
			if(in.first != in.second) {
				// ここで判定関数を取得
				// Innerに含まれる子オブジェクトは全て同じ型 という前提
				ColFunc cf = GetCFunc(in.first->getCID(), mdl1->getCID());
				const void* core1 = mdl1->getCore();
				do {
					if(cf(in.first->getCore(), core1)) {
						if(HitL(mdl1, in.first.get(), false))
							return true;
					}
				} while(++in.first != in.second);
			} else {
				if(bSwap)
					return false;
				return HitL(mdl1, mdl0, true);
			}
			return false;
		}
	};
	template <class CTG, class IM>
	ColFunc Narrow<CTG, IM>::cs_cfunc[ArraySize];
}
