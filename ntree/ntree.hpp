#pragma once
#include "spinner/bits.hpp"
#include "spinner/common.hpp"
#include "spinner/resmgr.hpp"
#include "../common.hpp"
#include <stack>

namespace boom {
	namespace ntree {
		using MortonId = uint32_t;
		using MortonId_OP = spn::Optional<MortonId>;
		using PositionId = uint32_t;
		struct VolEntry {
			spn::SHandle	hObj;
			//! 各軸の数値をビットフィールドで記録 (Dim依存)
			PositionId		posMin,
							posMax;
		};
		using VolList = spn::noseq_vec<VolEntry>;
		using VolVec = std::vector<VolEntry>;

		class CTreeObjStack;
		//! シングルラインCTreeエントリ
		struct CTreeEntry {
			//! 巡回時にはこのスタックを使う
			using ItrStack = CTreeObjStack;

			VolList		olist;			//!< セルに内包されたオブジェクトリスト
			int			nLower;			//!< このセルより下に幾つオブジェクトが存在するか

			CTreeEntry();
			void clear();
			bool remObj(spn::SHandle obj);
			bool isNodeEmpty() const;	//!< このノード単体が空か？
			bool isEmpty() const;		//!< 下位ノードを含めて空か？
		};

		//! 木を巡回する際に手持ちのオブジェクトを保持
		class CTreeObjStack {
			struct Ent {
				int		nPop;		//!< pop時に幾つpopするか
				int		baseIdx;
			};
			using Stack = std::stack<Ent>;

			Stack		_nstk;
			VolVec		_obj;

			public:
				CTreeObjStack();
				void addBlock(const CTreeEntry& ent, bool bAdd);
				void addBlock(const VolVec& ol, bool bAdd);

				/*! \param[in] bWirte 下層レイヤーのマス目有効フラグ
					\param[in] centerId */
				template <class DIM>
				void classify(const bool (&bWrite)[DIM::N_LayerSize], typename DIM::Id centerId) {
					auto cur = getObj();
					const int stride = std::get<1>(cur);
					// 処理途中で配列が再アロケートされることを防ぐ
					int wcur = _obj.size();		// 実際に書き込まれたオブジェクトを示すカーソル
					_obj.resize(wcur + stride * DIM::N_LayerSize);
					// 振り分けるポインタを一本の配列に一時的に格納
					std::unique_ptr<const VolEntry*>	da(new const VolEntry*[stride * DIM::N_LayerSize]);
					const VolEntry** daP[DIM::N_LayerSize];
					for(int i=0 ; i<DIM::N_LayerSize ; i++)
						daP[i] = &da.get()[i*stride];

					// オブジェクトを振り分け
					cur = getObj();
					DIM::Classify(daP, std::get<0>(cur), stride, centerId);

					// 枝の有効リストを見ながらVolEntry実体をコピー, タグ付け
					for(int i=0 ; i<DIM::N_LayerSize ; i++) {
						if(bWrite[i]) {
							auto base = da.get() + i*stride;

							// daPに幾つオブジェクトが割り振られたか
							int nObj = daP[i] - base;
							_nstk.push(Ent{nObj, wcur});
							for(int j=0 ; j<nObj ; j++) {
								AssertP(Trap, (*base)->hObj.valid())
								_obj[wcur++] = *(*base++);
							}
						}
					}
					// 本来の配列の長さに設定
					AssertP(Trap, wcur <= _obj.size())
					_obj.resize(wcur);
				}
				bool isTopEmpty() const;
				void popBlock();
				//! スタックトップのVolEntryリストを受け取る
				std::tuple<const VolEntry*, int> getObj() const;
		};

		class CTreeObjStackM;
		//! マップ用オブジェクトを含むCTreeエントリ
		struct CTreeEntryM : CTreeEntry {
			using ItrStack = CTreeObjStackM;
			VolList		mlist;

			void clear();
			bool isNodeEmpty() const;
			bool isEmpty() const;
		};
		class CTreeObjStackM {
			CTreeObjStack _stk,
						  _stkM;
			public:
				void addBlock(const CTreeEntryM& ent, bool bAdd);
				void popBlock();
				std::tuple<const VolEntry*, int> getObj() const;
				std::tuple<const VolEntry*, int> getMap() const;
				bool isTopEmpty() const;

				template <class DIM>
				void classify(const bool (&bWrite)[DIM::N_LayerSize], typename DIM::Id centerId) {
					_stk.classify<DIM>(bWrite, centerId),
					_stkM.classify<DIM>(bWrite, centerId);
				}
		};

		//! broad-phase collision manager (N-tree)
		/*!	属性フラグが0x80000000の時はBへ、それ以外はAに登録
			A->A, A->Bでは判定が行われるが B->Bはされない
			\tparam CTDim		(モートンID変換などを担当する)次元クラス
			\tparam CTEnt		エントリ管理クラステンプレート
			\tparam NDiv		分割数
			\tparam Ent			エントリーの中身
		*/
		template <class CTDim, template<class,int,int> class CTEnt, int NDiv, class Ent=CTreeEntry>
		class NTree {
			public:
				constexpr static int NDim = CTDim::N_Dim,
									N_LayerSize = CTDim::N_LayerSize;
				using Entry_t = Ent;
				using CTDim_t = CTDim;
				using CTEnt_t = CTEnt<Entry_t, NDiv, NDim>;
				using BVolume = typename CTDim_t::BVolume;
				using IDType = spn::SHandle;
				using IDTypeV = std::vector<IDType>;

				int _saturateW(float w) const {
					return Saturate(static_cast<int>(w*_unitSizeInv),
									0, CTEnt_t::N_Width-1);
				}

			private:
				struct MCEnt {
					MortonId	mortonId;
					BVolume		bvolume;
					uint32_t	posMin,
								posMax;
				};
				CTDim_t			_dim;
				CTEnt_t			_ent;
				const float		_unitSizeInv;

				using FGetBV = std::function<BVolume (spn::SHandle)>;
				const FGetBV	_fGetBV;
				using MortonCache = std::unordered_map<spn::SHandle, MCEnt>;
				MortonCache		_mmap;

				template <class Notify>
				int iterateChk(CMask cmask, const BVolume& bv, Notify&& ntf) const {
					if(getEntry(0).isEmpty())
						return 0;

					using Id = typename NTree::CTDim_t::Id;
					auto mid = CTDim_t::ToMortonId(bv, CTEnt_t::N_Width, _unitSizeInv);
					VolEntry ve{spn::SHandle(),
								std::get<2>(mid).value,
								std::get<3>(mid).value};

					int count = 0;
					struct Pair {
						int toProc;
						Id	center;
					};
					std::stack<Pair> stkId;
					// ルートノードをプッシュ
					int curWidth = NTree::CTEnt_t::N_Width/2;		// 現在の走査幅 (常に2の乗数)
					stkId.push(Pair{0, Id(curWidth)});
					while(!stkId.empty()) {
						auto p = stkId.top();
						stkId.pop();
						if(p.toProc < 0) {
							if(p.toProc == -2) {
								curWidth *= 2;
								if(curWidth == 0)
									curWidth = 1;
							}
							continue;
						}
						const auto& curEnt = getEntry(p.toProc);
						// 判定を行うかの判断 (assert含む)
						if(!curEnt.isNodeEmpty()) {
							for(auto& obj : curEnt.olist) {
								if(bv.hit(_mmap.at(obj.hObj).bvolume)) {
									std::forward<Notify>(ntf)(obj.hObj);
									++count;
								}
							}
						} else {
							AssertP(Trap, curEnt.nLower > 0)
						}
						// 下に枝があればそれを処理
						if(curEnt.nLower > 0) {
							stkId.push(Pair{-2, 0});
							curWidth /= 2;

							// 手持ちリストをOctave個のリストに分類(重複あり)にしてスタックに積む
							Id			lowerCenterId[NTree::N_LayerSize];
							NTree::CTDim_t::CalcLowerCenterId(lowerCenterId, p.center, curWidth);
							// 先に枝が有効かどうか確認
							int tcount = 0;
							for(int i=0 ; i<NTree::N_LayerSize ; i++) {
								int idx = p.toProc*NTree::N_LayerSize+1+i;		// 子ノードのインデックス
								if(hasEntry(idx)) {
									auto& ent = getEntry(idx);
									if(!ent.isEmpty()) {
										NTree::CTDim_t::Classify(ve, lowerCenterId[i], [i, idx2=idx, &lowerCenterId, &p, &stkId](const VolEntry& ve, int idx){
											stkId.push(Pair{idx2, lowerCenterId[i]});
										});
										++tcount;
										continue;
									}
								}
							}
							AssertP(Trap, tcount>0)
						}
					}
					return count;
				}
				//! コリジョン判定の為に木を巡る
				//! オブジェクト分類してaddBlockする
				/*! @param[out] dst 要素数はエントリに対応した数を用意
					@return コリジョン判定が実施された回数 */
				template <class Chk, class Notify>
				int iterate(Chk&& chk, Notify&& ntf) const {
					if(getEntry(0).isEmpty())
						return 0;

					int count = 0;
					typename NTree::CTEnt_t::Entry::ItrStack stk;		// 現在持っているオブジェクト集合

					using Id = typename NTree::CTDim_t::Id;
					struct Pair {
						int	toProc;			// これから処理する枝ノード. 負数は一段上に上がる
						Id center;			// 中心座標(ノードが負数の時は無効)

						Pair(int toproc, Id cent): toProc(toproc), center(cent) {}
					};
					std::stack<Pair> stkId;

					int curWidth = NTree::CTEnt_t::N_Width/2;		// 現在の走査幅 (常に2の乗数)

					// ルートノードをプッシュ
					stkId.push(Pair(0, Id(curWidth)));

					while(!stkId.empty()) {
						auto p = stkId.top();
						stkId.pop();
						if(p.toProc < 0) {
							if(p.toProc == -2) {
								stk.popBlock();		// オブジェクト集合の状態を1つ戻す
								curWidth *= 2;
								if(curWidth == 0)
									curWidth = 1;
							}
							continue;
						}

						const auto& curEnt = getEntry(p.toProc);
						// 判定を行うかの判断 (assert含む)
						if(!curEnt.isNodeEmpty()) {
							count += _doCollision(stk, curEnt,
									std::forward<Chk>(chk), std::forward<Notify>(ntf));
							// オブジェクト集合を加える
							stk.addBlock(curEnt, true);
						} else {
							AssertP(Trap, curEnt.nLower > 0)
						}

						// 下に枝があればそれを処理
						if(curEnt.nLower > 0) {
							stkId.push(Pair(-2, 0));
							curWidth /= 2;
							// 手持ちリストをOctave個のリストに分類(重複あり)にしてスタックに積む
							Id			lowerCenterId[NTree::N_LayerSize];
							bool		bWrite[NTree::N_LayerSize] = {};
							NTree::CTDim_t::CalcLowerCenterId(lowerCenterId, p.center, curWidth);
							// 先に枝が有効かどうか確認
							int tcount = 0;
							for(int i=0 ; i<NTree::N_LayerSize ; i++) {
								int idx = p.toProc*NTree::N_LayerSize+1+i;		// 子ノードのインデックス
								if(hasEntry(idx)) {
									auto& ent = getEntry(idx);
									if(!ent.isEmpty()) {
										stkId.push(Pair(idx, lowerCenterId[i]));
										++tcount;
										bWrite[i] = true;
										continue;
									}
								}
							}
							stk.template classify<typename NTree::CTDim_t>(bWrite, p.center);
							AssertP(Trap, tcount>0)
						} else
							stk.popBlock();		// オブジェクト集合の状態を1つ戻す
					}
					return count;
				}

				//! 上位セルのカウンタをインクリメント
				void _incrementUpper(MortonId num) {
					while(num != 0) {
						num = (num-1)>>(CTDim_t::N_Dim);
						AssertP(Trap, _ent.refEntry(num).nLower >= 0)
						_ent.increment(num);
					}
				}
				//! 上位セルのカウンタをデクリメント
				void _decrementUpper(MortonId num) {
					while(num != 0) {
						num = (num-1)>>(CTDim::N_Dim);
						AssertP(Trap, _ent.refEntry(num).nLower > 0)
						_ent.decrement(num);
					}
				}
				//! リスト総当たり判定
				template <class Itr, class T1, class Chk, class Notify>
				static int _HitCheck(Itr itr0, Itr itrE0, const T1& oL1, Chk&& chk, Notify&& ntf) {
					int count = 0;
					while(itr0 != itrE0) {
						auto& o0 = *itr0;
						for(auto& o1 : oL1) {
							if(std::forward<Chk>(chk)(o0.hObj, o1.hObj)) {
								// 通知など
								std::forward<Notify>(ntf)(o0.hObj, o1.hObj);
							}
						}
						++itr0;
						count += oL1.size();
					}
					return count;
				}
				//! 1つのリスト内で当たり判定
				template <class T0, class Chk, class Notify>
				static int _HitCheck(const T0& oL0, Chk&& chk, Notify&& ntf) {
					int count = 0;
					int nL = oL0.size();
					for(auto itr0=oL0.cbegin() ; itr0!=oL0.cend() ; ++itr0) {
						auto itr1=itr0;
						++itr1;
						for( ; itr1!=oL0.cend() ; ++itr1) {
							++count;
							if(std::forward<Chk>(chk)((*itr0).hObj, (*itr1).hObj)) {
								// 通知など
								std::forward<Notify>(ntf)((*itr0).hObj, (*itr1).hObj);
							}
						}
					}
					return count;
				}
				//! バウンディングボリュームの更新
				/*! \param[in] idv 非nullptrならリストに指定したオブジェクトのみの更新 */
				void _refreshBV(const IDTypeV* idv) {
					if(idv) {
						for(auto& idt : *idv) {
							auto itr = _mmap.find(idt);
							AssertP(Trap, itr!=_mmap.end())
							auto& cache = itr->second;
							auto bv = _fGetBV(itr->first);
							cache.bvolume = bv;
							auto mid = CTDim_t::ToMortonId(bv, CTEnt_t::N_Width, _unitSizeInv);
							if(std::get<2>(mid).value != cache.posMin ||
								std::get<3>(mid).value != cache.posMax)
							{
								// エントリ間の移動
								_moveObject(itr->first, false);
								cache.posMin = std::get<2>(mid).value;
								cache.posMax = std::get<3>(mid).value;
								cache.mortonId = ToMortonId(std::get<0>(mid), std::get<1>(mid));
							}
						}
					} else {
						for(auto& m : _mmap) {
							auto bv = _fGetBV(m.first);
							auto mid = CTDim_t::ToMortonId(bv, CTEnt_t::N_Width, _unitSizeInv);
							auto& cache = m.second;
							cache.bvolume = bv;
							if(std::get<2>(mid).value != cache.posMin ||
								std::get<3>(mid).value != cache.posMax)
							{
								// エントリ間の移動
								_moveObject(m.first, false);
								cache.posMin = std::get<2>(mid).value;
								cache.posMax = std::get<3>(mid).value;
								cache.mortonId = ToMortonId(std::get<0>(mid), std::get<1>(mid));
							}
						}
					}
				}
				MortonId _moveObject(const IDType& idt, bool bAddRem) {
					_remObject(idt, bAddRem);
					return _addObject(idt, bAddRem);
				}
				MortonId _addObject(spn::SHandle hObj, bool bAddMM) {
					AssertP(Trap, hObj.valid())
					auto bv = _fGetBV(hObj);
					auto mid = CTDim_t::ToMortonId(bv, CTEnt_t::N_Width, _unitSizeInv);
					MortonId num = ToMortonId(std::get<0>(mid), std::get<1>(mid));
					if(bAddMM)
						_mmap[hObj] = MCEnt{num, bv, std::get<2>(mid).value, std::get<3>(mid).value};
					// まだエントリがリストを持っていなければ自動的に作成
					auto& entry = _ent.refEntry(num);
					VolEntry ve{hObj, std::get<2>(mid).value, std::get<3>(mid).value};
					entry.olist.emplace_back(std::move(ve));
					_incrementUpper(num);
// 					LogOutput("Add(%1%) %2%", hObj.getIndex(), num);
					return num;
				}
				MortonId _remObject(const IDType& idt, bool bRemMM) {
					AssertP(Trap, idt.valid())
					auto itr = _mmap.find(idt);
					AssertP(Trap, itr != _mmap.end())
					MortonId num = itr->second.mortonId;
					if(bRemMM)
						_mmap.erase(itr);
					AssertP(Trap, _ent.hasEntry(num))

					auto& e = _ent.refEntry(num);
					// リストが空になったらエントリを削除
					if(e.remObj(idt))
						_ent.remEntry(num);
					_decrementUpper(num);
// 					LogOutput("Rem(%1%) %2%", idt.getIndex(), num);
					return num;
				}

			public:
				/*! \param[in] fieldSize	当たり判定対象の一片サイズ */
				NTree(const FGetBV& f, float fsize):
					_unitSizeInv(CTEnt_t::N_Width / fsize),
					_fGetBV(f)
				{}
				void clear() {
					_ent.clear();
				}
				IDType add(spn::SHandle hObj, CMask mask) {
					AssertP(Trap, hObj.valid())
					_addObject(hObj, true);
					return hObj;
				}
				void rem(const IDType& idt) {
					_remObject(idt, true);
				}
				template <class CB>
				void checkCollision(CMask mask, const BVolume& bv, CB&& cb) {
					_refreshBV(nullptr);
					iterateChk(mask, bv, std::forward<CB>(cb));
				}
				template <class CB>
				int broadCollision(CB&& cb, const IDTypeV* idv=nullptr) {
					_refreshBV(idv);
					return iterate([this](const IDType& id0, const IDType& id1){
									auto& bv0 = _mmap.at(id0).bvolume;
									auto& bv1 = _mmap.at(id1).bvolume;
									return bv0.hit(bv1);
								},
								std::forward<CB>(cb));
				}

				//! モートンIDとレベルから配列インデックスを計算
				static int MortonIdtoIndex(MortonId num, int level) {
					return static_cast<int>((std::pow(int(N_LayerSize), level)-1) / (N_LayerSize-1) + num);
				}
				//! BBox最小、最大地点のモートンIdからエントリを格納する場所を算出
				static MortonId ToMortonId(MortonId num0, MortonId num1) {
					auto diff = num0 ^ num1;
					num0 = MortonIdtoIndex(num0, CTEnt_t::N_Size);
					if(diff != 0) {
						auto upLv = (spn::Bit::MSB_N(diff) / CTDim_t::N_Dim)+1;
						do {
							num0 = (num0-1) >> CTDim_t::N_Dim;
						} while(--upLv != 0);
					}
					AssertP(Trap, num0 < (CTEnt_t::N_Ent))
					return num0;
				}
				const Entry_t& getEntry(MortonId num) const {
					return _ent.getEntry(num);
				}
				bool hasEntry(MortonId num) const {
					return _ent.hasEntry(num);
				}
				template <class Chk, class Notify>
				int _doCollision(const typename Entry_t::ItrStack& stk, const CTreeEntry& cur, Chk&& chk, Notify&& ntf) const {
					const VolList& ol = cur.olist;
					// Objリストとブロック内Objとの判定
					auto ret = stk.getObj();
					auto* bgn = std::get<0>(ret);
					auto nObj = std::get<1>(ret);
					int count = _HitCheck(bgn, bgn+nObj, ol, chk, ntf);
					// ブロック内同士の判定
					count += _HitCheck(ol, chk, ntf);
					return count;
				}
		};

		// -------- <<エントリ実装テンプレートクラス>> --------
		/*! \tparam NDiv	分割度
			\tparam Dim		次元数(2 or 3) */
		template <int NDiv, int Dim>
		class CTEnt_Base {
			public:
				constexpr static int N_Size = NDiv,
									N_Width = spn::NPow<2,N_Size>::result,
									N_LayerSize = spn::NPow<2,Dim>::result,
									N_Ent = (spn::NPow<N_LayerSize, NDiv>::result - 1) / (N_LayerSize-1) + spn::NPow<N_Width, Dim>::result;
			protected:
				CTEnt_Base() {
					// 2D木なら4, 3D木なら8. それ以外はエラー
					static_assert(N_LayerSize==4||N_LayerSize==8, "invalid dimension");
				}
		};
	}
}

