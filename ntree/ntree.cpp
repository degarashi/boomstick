#include "ntree.hpp"

namespace boom {
	namespace ntree {
		// ---------------------- CTreeEntry ----------------------
		CTreeEntry::CTreeEntry():
			nLower(0)
		{}
		void CTreeEntry::clear() {
			olist.clear();
			nLower = 0;
		}
		bool CTreeEntry::remObj(spn::SHandle obj) {
			auto itr = std::find_if(olist.cbegin(), olist.cend(),
						[obj](const VolEntry& e){ return e.hObj == obj; });
			AssertP(Trap, itr != olist.cend())
			olist.erase(itr);
			return isEmpty();
		}
		bool CTreeEntry::isNodeEmpty() const {
			return olist.empty();
		}
		bool CTreeEntry::isEmpty() const {
			return isNodeEmpty() && nLower==0;
		}
		// ---------------------- CTreeEntryM ----------------------
		void CTreeEntryM::clear() {
			mlist.clear();
			CTreeEntry::clear();
		}
		bool CTreeEntryM::isNodeEmpty() const {
			return olist.empty() && mlist.empty();
		}
		bool CTreeEntryM::isEmpty() const {
			return CTreeEntry::isEmpty() && mlist.empty();
		}

		// ---------------------- CTreeObjStack ----------------------
		CTreeObjStack::CTreeObjStack() {
			_nstk.push(Ent(0,0));
		}
		void CTreeObjStack::addBlock(const CTreeEntry& ent, bool bAdd) {
			addBlock(ent.olist, bAdd);
		}
		void CTreeObjStack::popBlock() {
			auto top = _nstk.top();
			_nstk.pop();

			_obj.erase(_obj.end()-top.nPop, _obj.end());
		}
		bool CTreeObjStack::isTopEmpty() const {
			return int(_obj.size()) <= _nstk.top().baseIdx;
		}
		std::tuple<const VolEntry*,int> CTreeObjStack::getObj() const {
			if(_obj.empty())
				return std::make_tuple((const VolEntry*)nullptr, 0);

			int bi = _nstk.top().baseIdx;
			return std::make_tuple(&_obj[0] + bi, int(_obj.size())-bi);
		}

		// ---------------------- CTreeObjStackM ----------------------
		void CTreeObjStackM::addBlock(const CTreeEntryM& ent, bool bAdd) {
			_stk.addBlock(ent.olist, bAdd);
			_stkM.addBlock(ent.mlist, bAdd);
		}
		void CTreeObjStackM::popBlock() {
			_stk.popBlock();
			_stkM.popBlock();
		}
		std::tuple<const VolEntry*, int> CTreeObjStackM::getObj() const {
			return _stk.getObj();
		}
		std::tuple<const VolEntry*, int> CTreeObjStackM::getMap() const {
			return _stkM.getObj();
		}
		bool CTreeObjStackM::isTopEmpty() const {
			return _stk.isTopEmpty() && _stkM.isTopEmpty();
		}
	}
}
