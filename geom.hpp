#pragma once
#include "spinner/type.hpp"
#include "spinner/misc.hpp"
#include "spinner/matrix.hpp"
#include "spinner/plane.hpp"

namespace boom {
	using spn::Vec2;
	using spn::Vec3;
	using spn::Mat33;
	using spn::Mat43;
	using spn::Mat44;
	using spn::AMat44;
	using spn::AMat43;
	using spn::Plane;

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
			mutable uint8_t		_buff[CTRET::template SumN<CTRET::size-1, GetSize_Ref0>::result];
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
				if(!Bit::Check(_cacheFlag, FlagGet0<T>::get())) {
					Bit::Set(_cacheFlag, FlagGet0<T>::get());
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
				spn::Bit::Clear(_cacheFlag, spn::CType<Tags...>::template SumN<sizeof...(Tags)-1, FlagGet0>::result);
			}
		public:
			template <class T>
			const T_OP<T>& _getInfo() const {
				constexpr int idx = CTTAG::template Find<T>::result;
				return _getInfo<T>(std::false_type());
			}
	};
}
