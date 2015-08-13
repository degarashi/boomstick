#include "ntree.hpp"

namespace boom {
	namespace ntree {
		// ---------------------- CTreeEntry ----------------------
		CTreeEntry::CTreeEntry():
			_nLower(0)
		{}
		void CTreeEntry::clear() {
			_olist.clear();
			_nLower = 0;
		}
		void CTreeEntry::addObj(CMask /*mask*/, const VolEntry& ve) {
			_olist.emplace_back(ve);
		}
		void CTreeEntry::_remObj(VolVec& v, CacheId cid) {
			auto itr = std::find_if(v.begin(), v.end(),
						[cid](const VolEntry& e){ return e.cacheId == cid; });
			AssertP(Trap, itr != v.cend())
			v.erase(itr);
		}
		bool CTreeEntry::remObj(CMask /*mask*/, CacheId cid) {
			_remObj(_olist, cid);
			return isEmpty();
		}
		bool CTreeEntry::isNodeEmpty() const {
			return _olist.empty();
		}
		bool CTreeEntry::isEmpty() const {
			return isNodeEmpty() && _nLower==0;
		}
		const VolVec& CTreeEntry::getObjList() const {
			return _olist;
		}
		int CTreeEntry::getLowerCount() const {
			return _nLower;
		}
		void CTreeEntry::incrementLowerCount() {
			++_nLower;
			AssertP(Trap, _nLower > 0)
		}
		void CTreeEntry::decrementLowerCount() {
			--_nLower;
			AssertP(Trap, _nLower >= 0)
		}
		// ---------------------- CTreeEntryM ----------------------
		bool CTreeEntryM::IsTypeB(CMask mask) {
			return mask & 0x80000000;
		}
		void CTreeEntryM::addObj(CMask mask, const VolEntry& v) {
			if(IsTypeB(mask))
				_mlist.emplace_back(v);
			else
				base_t::addObj(mask, v);
		}
		bool CTreeEntryM::remObj(CMask mask, CacheId cid) {
			if(IsTypeB(mask))
				_remObj(_mlist, cid);
			else
				base_t::remObj(mask, cid);
			return isEmpty();
		}
		void CTreeEntryM::clear() {
			_mlist.clear();
			base_t::clear();
		}
		bool CTreeEntryM::isNodeEmpty() const {
			return base_t::isNodeEmpty() && _mlist.empty();
		}
		bool CTreeEntryM::isEmpty() const {
			return base_t::isEmpty() && _mlist.empty();
		}
		const VolVec& CTreeEntryM::getMapList() const {
			return _mlist;
		}
		// ---------------------- CTreeObjStack ----------------------
		CTreeObjStack::CTreeObjStack() {
			_nstk.push(Ent{0,0});
		}
		void CTreeObjStack::addBlock(const CTreeEntry& ent, bool bAdd) {
			addBlock(ent.getObjList(), bAdd);
		}
		void CTreeObjStack::addBlock(const VolVec& ol, bool bAdd) {
			int nAdd = ol.size();
			if(bAdd)
				_nstk.top().nPop += nAdd;
			else
				_nstk.push(Ent{nAdd, int(_obj.size())});

			int cur = _obj.size();
			_obj.resize(cur + nAdd);
			for(auto& obj : ol)
				_obj[cur++] = obj;
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
			_stk.addBlock(ent.getObjList(), bAdd);
			_stkM.addBlock(ent.getMapList(), bAdd);
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
