#pragma once
#include "spinner/type.hpp"

namespace boom {
	using spn::SelectType;

	template <class CORE>
	class Cache {
		using CT = typename CORE::CT;
		using CTD = typename CT::template Another<>::result;
		template <class T>
		using Wrapper = typename CORE::template Wrap<T>;

		mutable uint32_t				_rflag = ~0;
		mutable typename CTD::AsTuple	_data;
		protected:
			Wrapper<CORE>	_core;

			void setCacheFlag() {}
			template <class T, class... Ts>
			void setCacheFlag() {
				_rflag |= (1<<CT::template Find<T>::result);
				setCacheFlag<Ts...>();
			}
			void setCacheFlagAll() {
				_rflag = ~0;
			}

		public:
			Cache() {}
			template <class... Ts>
			Cache(Ts&&... ts): _core(std::forward<Ts>(ts)...) {}

			struct FlagTrue {};
			struct FlagFalse {};
			template <class TAG>
			struct TV {
				constexpr static int POS = CT::template Find<TAG>::result,
									FLAG = 1<<POS;
			};
			template <class TAG>
			using flag_type = typename SelectType<CT::template Has<TAG>::result, FlagTrue, FlagFalse>::type;
			template <class TAG>
			using TypeD = decltype(TAG().get(Wrapper<CORE>()));
			template <class TAG>
			using TypeC = const decltype(TAG()())&;
			template <class TAG>
			using Detect = typename spn::SelectTypeT<CT::template Has<TAG>::result, TypeC, TypeD>::template type<TAG>;

			const CORE& getCore() const { return _core; }
			template <class TAG>
			decltype(std::get<TV<TAG>::POS>(_data))& refCache(TAG) const {
				_rflag &= ~TV<TAG>::FLAG;
				return std::get<TV<TAG>::POS>(_data);
			}

			// キャッシュ有効
			template <class TAG>
			TypeC<TAG> _getCache(FlagTrue) const {
				if(_rflag & TV<TAG>::FLAG)
					_core.getInfo(*this, TAG());
				return std::get<TV<TAG>::POS>(_data);
			}
			// キャッシュ無効
			template <class TAG>
			TypeD<TAG> _getCache(FlagFalse) const {
				return TAG().get(_core);
			}

			template <class TAG>
			Detect<TAG> getCache(TAG) const {
				return _getCache<TAG>(flag_type<TAG>());
			}
	};
	//! ベースクラスにキャッシュ対象値の取得メソッドを加える
	template <class BASE>
	struct CoreWrap : BASE {
		CoreWrap() {}
		template <class... Ts>
		CoreWrap(Ts&&... ts): BASE(std::forward<Ts>(ts)...) {}
		using BASE::getInfo;
		template <class INFO, class TAG>
		void getInfo(INFO& info, TAG tag) const {
			info.refCache(tag) = tag.get(*this);
		}
		operator const BASE& () const { return *this; }
	};
	template <class BASE>
	struct CoreRaw : BASE {
		CoreRaw() {}
		template <class... Ts>
		CoreRaw(Ts&&... ts): BASE(std::forward<Ts>(ts)...) {}
		template <class INFO, class TAG>
		void getInfo(INFO& info, TAG tag) const {
			info.refCache(tag) = tag.get(*this);
		}
		operator const BASE& () const { return *this; }
	};
}