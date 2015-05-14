#pragma once
#include "ntree.hpp"

namespace boom {
	namespace ntree {
		//! マップ用判定と通常オブジェクトを分けたコリジョンクラス
		/*!	属性フラグが0x80000000の時はBへ、それ以外はAに登録
			A->A, A->Bでは判定が行われるが B->Bはされない
			\tparam CTDim		(モートンID変換などを担当する)次元クラス
			\tparam CTEnt		エントリ管理クラステンプレート
			\tparam NDiv		分割数
		*/
		/*! マップ用判定同士は判定しない */
		template <class CTDim, template <class,int,int> class CTEnt, int NDiv>
		class NTreeMap : public _NTree<CTDim,
									CTEnt,
									NDiv,
									CTreeEntryM,
									NTreeMap<CTDim,CTEnt, NDiv>>
		{
			private:
				template <class, template<class,int,int> class, int, class, class>
				friend class _NTree;
				using base_t = _NTree<CTDim, CTEnt, NDiv, CTreeEntryM, NTreeMap<CTDim,CTEnt,NDiv>>;
				using this_t = NTreeMap;

			protected:
				template <class Notify>
				int _doCollision(CMask mask, const typename base_t::BVolume& bv, const CTreeEntryM& cur, const Notify& ntf) const {
					int count = 0;
					auto fnChk = [&](const VolVec& v){
						for(auto& obj : v) {
							auto& c = base_t::_getCache(obj.cacheId);
							if((mask & c.mask) &&
								bv.hit(c.bvolume))
							{
								ntf(c.hObj);
								++count;
							}
						}
					};
					fnChk(cur.getObjList());
					if(!CTreeEntryM::IsTypeB(mask))
						fnChk(cur.getMapList());
					return count;
				}
				template <class Notify>
				int _doCollision(const typename CTreeEntryM::ItrStack& stk, const CTreeEntryM& cur, const Notify& ntf) const {
					const VolVec &ol = cur.getObjList(),
								&ml = cur.getMapList();
					auto fnChk = [this](const VolEntry& v0, const VolEntry& v1){
						auto	&c0 = base_t::_getCache(v0.cacheId),
								&c1 = base_t::_getCache(v1.cacheId);
						return (c0.mask & c1.mask) && c0.bvolume.hit(c1.bvolume);
					};
					auto fnNtf = [this, &ntf](const VolEntry& v0, const VolEntry& v1){
						auto	&c0 = base_t::_getCache(v0.cacheId),
								&c1 = base_t::_getCache(v1.cacheId);
						ntf(c0.hObj, c1.hObj);
					};
					// Objリストとブロック内Objとの判定
					auto ret = stk.getObj();
					int count = base_t::_HitCheck(std::get<0>(ret), std::get<0>(ret)+std::get<1>(ret), ol, fnChk, fnNtf);
					// ObjリストとマップObjとの判定
					count += base_t::_HitCheck(std::get<0>(ret), std::get<0>(ret)+std::get<1>(ret), ml, fnChk, fnNtf);
					auto ret2 = stk.getMap();
					count += base_t::_HitCheck(std::get<0>(ret2), std::get<0>(ret2)+std::get<1>(ret2), ol, fnChk, fnNtf);
					if(!ol.empty())
						count += base_t::_HitCheck(ol.cbegin(), ol.cend(), ml, fnChk, fnNtf);
					// ブロック内Obj同士の判定
					count += base_t::_HitCheck(ol, fnChk, fnNtf);
					return count;
				}
			public:
				using base_t::base_t;
		};
	}
}
