#pragma once
/*! 2D, 3Dを問わないコリジョンに関する定義 */
#include "spinner/optional.hpp"

namespace boom {
	//! サポートされていない関数を読んだ時の実行時エラーを送出
	#define INVOKE_ERROR Assert(Trap, false, "not supported function: ", __func__) throw 0;
	//! IModelから指定の型にキャスト出来ればその参照を返す関数のデフォルト実装
	#define DEF_CASTFUNC(typ) virtual spn::Optional<typ&> castAs##typ() { return spn::none; } \
			virtual spn::Optional<const typ&> castAs##typ() const { auto ref = const_cast<IModel*>(this)->castAs##typ(); \
			if(ref) return *ref; return spn::none; }
	template <class T, class CT>
	struct ITagP_base {
		static uint32_t GetCID() { return CT::template Find<T>::result; }
	};
	template <class T, class CT, class MDL>
	struct IModelP_base : ITagP_base<T,CT>, MDL {
		virtual uint32_t getCID() const override { return ITagP_base<T,CT>::GetCID(); }
	};

	//! ラインの位置を示す
	enum class LINEPOS {
		Begin,		//!< 始点
		End,		//!< 終点
		OnLine,		//!< ライン上
		NotHit
	};
}
