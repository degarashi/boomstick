#pragma once
#include "ntree.hpp"

namespace boom {
	namespace ntree {
		//! 配列によるエントリ実装
		/*! \tparam Ent		エントリクラス
			\tparam NDiv	分割度
			\tparam Dim		次元数(2 or 3) */
		template <class Ent, int NDiv, int Dim>
		class CTEnt_Array : public CTEnt_Base<NDiv, Dim> {
			public:
				using base_t = CTEnt_Base<NDiv, Dim>;
				using Entry = Ent;
			private:
				Entry	_ent[base_t::N_Ent];	//!< 線形木
			public:
				bool hasEntry(MortonId n) const {
					AssertP(Trap, n<countof(_ent))
					// 配列なので常にエントリは存在する
					return true;
				}
				const Entry& getEntry(MortonId n) const {
					AssertP(Trap, n<countof(_ent))
					return _ent[n];
				}
				Entry& refEntry(MortonId n) {
					AssertP(Trap, n<countof(_ent))
					return _ent[n];
				}
				void remEntry(MortonId n) {
					// 固定配列なので削除しない
					AssertP(Trap, _ent[n].isEmpty())
				}
				void increment(MortonId num) {
					AssertP(Trap, _ent[num].nLower >= 0)
					++_ent[num].nLower;
				}
				void decrement(MortonId num) {
					// カウンタが0になっても何もしない
					AssertP(Trap, _ent[num].nLower > 0)
					--_ent[num].nLower;
				}
				void clear() {
					for(auto& e : _ent)
						e.clear();
				}
		};
	}
}
