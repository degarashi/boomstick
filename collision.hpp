#pragma once
#include "spinner/resmgr.hpp"

namespace boom {
	//! 当たり判定属性マスク
	using CMask = uint32_t;
	using CTime = uint32_t;
	using Int_OP = spn::Optional<int>;

	/*! コリジョン情報を纏めた構造体
		\tparam CMGR	コリジョンマネージャ
		\tparam BV		BoundingVolume Type
		\tparam HLMDL	モデルハンドル(Locked)
		\tparam BCID	BroadCollisionクラスのオブジェクトを特定できるようなデータ型
		\tparam UD		任意のユーザーデータ型
	*/
	template <class CMGR, class BV, class HLMDL, class BCID, class UD>
	class ColMem {
		private:
			CMask	_mask;	//!< 当たり判定対象フラグ
			struct Index {
				CTime	time;	//!< indexが指している累積時間
				Int_OP	front,	//!< 線形リスト先頭インデックス
						last;	//!< 線形リスト末尾インデックス

				//! リストに何もつないでない状態に初期化
				Index(CTime tm):
					time(tm)
				{}
				//! リスト末尾に新しくインデックスを加える
				/*!	\param[in] idx 新しく加えるインデックス値
					\return 直前の先頭インデックス値 */
				Int_OP setLastIndex(int idx) {
					Int_OP ret = last;
					if(!front)
						front = last = idx;
					else
						last = idx;
					return ret;
				}
			};
			using Index_OP = spn::Optional<Index>;
			Index_OP	_cur,
						_prev;
			HLMDL		_hlMdl;
			BCID		_bcid;
			UD			_udata;

		public:
			template <class UD2>
			ColMem(CMask ms, HLMDL hm, UD2&& ud):
				_mask(ms),
				_hlMdl(hm),
				_udata(std::forward<UD2>(ud))
			{}
			void setBCID(const BCID& ns) {
				_bcid = ns;
			}
			const BCID& getBCID() const {
				return _bcid;
			}
			BV getBVolume(CTime t) const {
				BV bv;
				auto& r = _hlMdl.cref();
				AssertP(Trap, r->imn_refresh(t))
				r->im_getBVolume(bv);
				return bv;
			}
			decltype(_hlMdl.get()) getModelHandle() const {
				return _hlMdl.get();
			}
			auto getModel() const -> decltype(_hlMdl.cref().get()) {
				return _hlMdl.cref().get();
			}
			CMask getMask() const {
				return _mask;
			}
			UD& refUserData() {
				return _udata;
			}
			//! 衝突判定履歴を巡回 (コリジョン開始 / 継続中の判定)
			/*! \param cb[in,out] コールバック(CMgr::Hist) */
			template <class CB>
			void getCollision(CB&& cb) const {
				auto& m = CMGR::_ref();
				CTime tm = m.getAccum();
				if(_cur) {
					auto& cur = *_cur;
					// 現在の時刻と等しい時に巡回
					if(cur.time == tm) {
						m.iterateHistCur(*cur.front, [&cb](const auto& h){
							std::forward<CB>(cb)(h);
						});
					}
				}
			}
			//! 衝突判定結果を参照 (コリジョン終了判定)
			/*! \param cb[in,out] コールバック(CMgr::Hist) */
			template <class CB>
			void getEndCollision(CB&& cb) const {
				auto& m = CMGR::_ref();
				CTime tm = m.getAccum();
				if(_prev && _prev->time == tm-1) {
					if(_cur->time == tm) {
						// 前回何かと衝突していて前回も判定済の場合
						// _prevにあって_curに無い物が対象
						int prev = *(_prev->front),
							cur = *(_cur->front);
						m.iterateHistPrev(prev, [cur, &m, &cb](const auto& hPrev){
							bool bFound = false;
							m.iterateHistCur(cur, [&bFound, &hPrev](const auto& hCur){
								if(hCur.hCol == hPrev.hCol) {
									bFound = true;
								}
							});
							if(!bFound)
								std::forward<CB>(cb)(hPrev);
						});
					}
				}
				if(_cur && _cur->time == tm-1) {
					// 今回何とも衝突しなかった場合
					// _cur のリスト全部対象。ただし履歴はColMgr上では過去のものになってるのでPrevを参照
					m.iterateHistPrev(*_cur->front, std::forward<CB>(cb));
				}
			}
			//! 引数の時刻-1と合致する方のリストインデックス先頭を取得
			/*! \return どちらも無効ならspn::none */
			Int_OP getPrevIndex(CTime tm) const {
				if(_cur && _cur->time == tm-1)
					return _cur->front;
				if(_prev && _prev->time == tm-1)
					return _prev->front;
				return spn::none;
			}
			//! カレントリスト末尾に新しくインデックスを加える
			/*! \return 直前の先頭インデックス値 */
			Int_OP setLastIndex(CTime tm, int idx) {
				// もし新しい時刻の設定だったらカーソルの新旧入れ替え
				if(_cur) {
					if(_cur->time < tm) {
						_prev = _cur;
						_cur = spn::construct(tm);
					}
				} else
					_cur = spn::construct(tm);
				return _cur->setLastIndex(idx);
			}
	};
	/*!	\tparam	BC		BroadCollision class
		\tparam	Types	geom2d::Types or geom3d::Types
		\tparam UD		userdata type
	*/
	template <class BC,
				class Types,
				class UD>
	class ColMgr : public spn::ResMgrA<ColMem<ColMgr<BC,Types,UD>,
											typename BC::BVolume,
											typename Types::MMgr::LHdl,
											typename BC::IDType,
											UD>,
										ColMgr<BC,Types,UD>>,
					public Types::Narrow
	{
		public:
			using IModel = typename Types::IModel;									//!< IModel interface
			using MMgr = typename Types::MMgr;										//!< ModelManager
			using user_t = UD;														//!< Userdata
			using BVolume = typename BC::BVolume;									//!< Bounding-volume
			using this_t = ColMgr<BC, Types, user_t>;
			using BroadC = BC;														//!< BroadCollision class
			using BCId = typename BroadC::IDType;									//!< BroadCollision dependant Id
			using BCIdV = typename BroadC::IDTypeV;
			// ---- model handle type ----
			using HMdl = typename MMgr::SHdl;
			using HLMdl = typename MMgr::LHdl;
			using MdlRef = decltype(std::declval<HMdl>().ref());
			using MdlP = decltype(std::declval<MdlRef>().get());
			using CMem = ColMem<this_t, BVolume, HLMdl, BCId, user_t>;						//!< Collision information structure

			using base = spn::ResMgrA<CMem, this_t>;
			// --- collision handle type ----
			using HCol = typename base::SHdl;
			using HLCol = typename base::LHdl;
			// 単方向リスト
			struct Hist {
				HCol	hCol;
				int		nFrame;				//!< 衝突継続したフレーム数
				int		nextOffset;			//!< 次のエントリへのバイトオフセット

				Hist(HCol hc, int nf):
					hCol(hc),
					nFrame(nf),
					nextOffset(0)
				{}
			};
		private:
			using Narrow = typename Types::Narrow;
			using HistVec = std::vector<Hist>;
			using RemList = std::vector<HLCol>;

			BroadC		_broadC;
			CTime		_accum = 0;
			//! history, removelistの現在のインデックス
			int			_swHist = 0;
			HistVec		_hist[2];
			RemList		_rmList[2];
			bool		_bDelete = false;

			HistVec& _getCurHist() { return _hist[_swHist]; }
			HistVec& _getPrevHist() { return _hist[_swHist^1]; }
			const HistVec& _getCurHist() const { return _hist[_swHist]; }
			const HistVec& _getPrevHist() const { return _hist[_swHist^1]; }
			RemList& _getCurRM() { return _rmList[_swHist]; }
			RemList& _getPrevRM() { return _rmList[_swHist^1]; }
			using Hist_OP = spn::Optional<Hist&>;

			//! hc0のリストにhc1のエントリが存在するかチェック
			/*!	\return hc1が見つかればその参照を返す */
			Hist_OP _hasPrevCollision(HCol hc0, HCol hc1) {
				auto& c0 = hc0.cref();
				Hist_OP ret;
				if(auto pidx = c0.getPrevIndex(_accum)) {
					iterateHistPrev(*pidx, [hc1,&ret](Hist& h){
						if(h.hCol == hc1)
							ret = h;
					});
				}
				return ret;
			}
			/*! \param hv	_hist[0] or _hist[1]
				\param idx	巡回を開始するインデックス */
			template <class HVec, class CB>
			static void _IterateHist(HVec& hv, int idx, CB&& cb) {
				AssertP(Trap, idx>=0)
				auto* hst = &hv[idx];
				for(;;) {
					std::forward<CB>(cb)(*hst);
					int np = hst->nextOffset;
					if(np != 0)
						hst = reinterpret_cast<decltype(hst)>(reinterpret_cast<uintptr_t>(hst) + np);
					else break;
				}
			}
			HistVec& _switchHist() {
				_doRelease(false);
				_swHist ^= 1;
				auto& cur = _getCurHist();
				cur.clear();
				return cur;
			}
			void _doRelease(bool bBoth) {
				struct Flag {
					bool& flag;
					Flag(bool& b): flag(b) { flag=true; }
					~Flag() { flag=false; }
				};
				//! 解放を一時的に遅延させていたハンドルを処理
				Flag f(_bDelete);
				_getPrevRM().clear();
				if(bBoth)
					_getCurRM().clear();
			}
			void _makeHist(HCol hc0, HCol hc1) {
				auto& hist = _getCurHist();
				auto fn = [this, &hist](HCol hc0, HCol hc1) {
					auto &c0 = hc0.ref();
					// 前のフレームでの継続フレームリストを探索
					int nf = 0;
					if(auto hist = _hasPrevCollision(hc0, hc1))
						nf = hist->nFrame + 1;
					int nextID = hist.size();
					Int_OP lastID = c0.setLastIndex(_accum, nextID);
					if(lastID)
						hist[*lastID].nextOffset = (nextID - *lastID) * sizeof(Hist);
					hist.push_back(Hist(hc1, nf));
				};
				fn(hc0, hc1);
				fn(hc1, hc0);
			}
			bool release(spn::SHandle sh) override {
				release(HCol::FromHandle(sh));
				return false;
			}
			auto _makeGetBV() {
				return [this](spn::SHandle sh) {
					auto& r = HCol::FromHandle(sh).cref();
					return r.getBVolume(_accum);
				};
			}
		public:
			ColMgr():
				_broadC(_makeGetBV())
			{
				Narrow::Initialize();
			}
			~ColMgr() {
				_doRelease(true);
			}
			CTime getAccum() const {
				return _accum;
			}
			template <class CB>
			void iterateHistCur(int idx, CB&& cb) const { _IterateHist(_getCurHist(), idx, std::forward<CB>(cb)); }
			template <class CB>
			void iterateHistPrev(int idx, CB&& cb) const { _IterateHist(_getPrevHist(), idx, std::forward<CB>(cb)); }
			template <class CB>
			void iterateHistCur(int idx, CB&& cb) { _IterateHist(_getCurHist(), idx, std::forward<CB>(cb)); }
			template <class CB>
			void iterateHistPrev(int idx, CB&& cb) { _IterateHist(_getPrevHist(), idx, std::forward<CB>(cb)); }

			//! 当たり判定対象を追加
			/*! \param ms[in] マスク */
			template <class UD2>
			HLCol addCol(CMask ms, HMdl hm, UD2&& ud=UD()) {
				HLCol hlC = base::emplace(ms, hm, std::forward<UD2>(ud));
				hlC.ref().setBCID(_broadC.add(hlC, ms));
				AssertP(Trap, hm->get()->imn_refresh(_accum), "empty object detected")
				return std::move(hlC);
			}
			//! ハンドル解放処理を一時的に遅延させる (デストラクタ時以外)
			void release(HCol hCol) {
				if(!_bDelete && hCol.count() == 1) {
					// OnCollisionEndの処理がある為、次のupdateまでハンドル解放を遅延
					_getCurRM().push_back(hCol);
					// しかしBroadCからは即刻外す
					_broadC.rem(hCol.cref().getBCID());
				}
				base::release(hCol);
			}
			/*!	\param[in] mask		コリジョンマスク値
				\param[in] mp		判定対象のモデルポインタ
				\param[in] cb		コールバック関数(HCol) */
			template <class CB>
			void checkCollision(CMask ms, MdlP mp, CB&& cb) {
				if(!mp->imn_refresh(_accum))
					return;
				BVolume bv;
				mp->im_getBVolume(bv);
				_broadC.checkCollision(ms, bv,
					[ac=_accum,
						pMdl=mp,
						cb=std::forward<CB>(cb)](spn::SHandle sh)
					{
						auto hc = HCol::FromHandle(sh);
						// 詳細判定
						if(Narrow::Hit(pMdl, hc->getModel(), ac))
							cb(hc);
					}
				);
			}
			template <class CB>
			void checkCollision(CMask ms, HMdl hMdl, CB&& cb) {
				checkCollision(ms, hMdl->get(), std::forward<CB>(cb));
			}
			//! 全てのオブジェクトを相互に当たり判定
			/*!	A->Bの組み合わせでコールバックを呼び、B->Aは呼ばない
				\param[in] cb		コールバック関数(HCol, HCol)
				\param[in] bAdvance	trueの時に時刻を進め、新旧の履歴を切り替える
				\param[in] idv		非nullptrの時はリストに登録してあるオブジェクトのみ更新 */
			template <class CB>
			void update(CB&& cb, bool bAdvance, const BCIdV* idv=nullptr) {
				if(bAdvance) {
					_switchHist();
					++_accum;
				}
				_broadC.broadCollision([bAdvance, this, cb=std::forward<CB>(cb)](spn::SHandle sh0, spn::SHandle sh1){
					// 同じハンドルという事はあり得ない筈
					AssertP(Trap, sh0!=sh1)

					HCol hc0(HCol::FromHandle(sh0)),
						hc1(HCol::FromHandle(sh1));
					// 詳細判定
					// 毎フレームヒットリストを再構築
					if(Narrow::Hit(hc0->getModel(), hc1->getModel(), _accum)) {
						if(bAdvance)
							_makeHist(hc0, hc1);
						cb(hc0, hc1);
					}
				}, idv);
			}
			//! 時間を1進め、衝突履歴の更新
			void update() {
				update([](auto,auto){}, true);
			}
	};
}
