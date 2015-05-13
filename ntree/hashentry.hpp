#pragma once
#include "ntree.hpp"

namespace boom {
	namespace ntree {
		//! ハッシュマップによるエントリ実装
		/*! \tparam Ent		エントリクラス
			\tparam NDiv	分割度
			\tparam NL		次元数(2 or 3) */
		template <class Ent, int NDiv, int Dim>
		class CTEnt_Hash : public CTEnt_Base<NDiv,Dim> {
			public:
				using Entry = Ent;
			private:
				using ObjHash = std::unordered_map<MortonId, Entry>;
				ObjHash		_ent;
			public:
				CTEnt_Hash() {
					// ルートノードだけは作成しておく
					_ent[0];
				}
				bool hasEntry(MortonId n) const {
					return _ent.count(n) == 1;
				}
				const Entry& getEntry(MortonId n) const {
					return _ent.at(n);
				}
				Entry& refEntry(MortonId n) {
					return _ent[n];
				}
				void remEntry(MortonId n) {
					auto itr = _ent.find(n);
					AssertP(Trap, itr->second.isEmpty())
					if(n != 0)
						_ent.erase(itr);
				}
				void increment(MortonId num) {
					AssertP(Trap, _ent[num].getLowerCount() >= 0)
					_ent[num].incrementLowerCount();
				}
				void decrement(MortonId num) {
					// カウンタが0になったらエントリを削除 (ルートは消さない)
					auto itr = _ent.find(num);
					AssertP(Trap, itr!=_ent.end() && itr->second.getLowerCount()>0)
					itr->second.decrementLowerCount();
					if(itr->second.isEmpty() && num!=0) {
						AssertP(Trap, itr->second.getObjList().empty())
						_ent.erase(itr);
					}
				}
				void clear() {
					_ent.clear();
				}
		};
	}
}
