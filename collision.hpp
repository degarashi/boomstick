#pragma once
#include "spinner/resmgr.hpp"

namespace boom {
	//! 当たり判定属性マスク
	using CMask = uint32_t;
	using CTime = uint32_t;

	template <class CMGR, class HLMDL, class NSID, class UD>
	class ColMem {
		CMask	_mask;	//!< 当たり判定対象フラグ
		struct Index {
			CMask 	time;	//!< indexが指している累積時間
			int		front,	//!< 線形リスト先頭インデックス
					last;	//!< 線形リスト末尾インデックス

			//! リストに何もつないでない状態に初期化
			void init(CTime tm) {
				time = tm;
				front = last = -1;
			}
			//! リスト末尾に新しくインデックスを加える
			/*!	\param[in] idx 新しく加えるインデックス値
				\return 直前の先頭インデックス値 */
			int add(int idx) {
				int ret = last;
				if(front < 0)
					front = last = idx;
				else
					last = idx;
				return ret;
			}
		};
		Index		_cur,
					_prev;
		HLMDL		_hlMdl;
		NSID		_nsid;
		UD			_udata;

		public:
			ColMem(CMask ms, HLMDL hm, const UD& ud): _mask(ms), _hlMdl(hm), _udata(ud) {
				_cur.init(0);
				_prev.init(0);
			}
			void setNSID(const NSID& ns) {
				_nsid = ns;
			}
			const NSID& getNSID() const {
				return _nsid;
			}
			auto getBVolume() const -> decltype(_hlMdl.cref()->im_getBVolume()) {
				return _hlMdl.cref()->im_getBVolume();
			}
			auto getModel() const -> decltype(_hlMdl.cref().get()) {
				return _hlMdl.cref().get();
			}
			//! 衝突判定結果を参照 (コリジョン開始 / 継続中の判定)
			/*! \param cb[in,out] コールバック(CMgr::Hist) */
			template <class CB>
			void getCollision(CB cb) const {
				auto& m = CMGR::_ref();
				CTime tm = m.getAccum();
				if(_cur.time >= tm && _cur.front>=0)
					m.iterateHistCur(_cur.front, cb);
			}
			//! 衝突判定結果を参照 (コリジョン終了判定)
			/*! \param cb[in,out] コールバック(CMgr::Hist) */
			template <class CB>
			void getEndCollision(CB cb) const {
				auto& m = CMGR::_ref();
				CTime tm = m.getAccum();
				if(_cur.time == tm) {
					// 今回何かと衝突した場合
					if(_prev.time == tm-1 && _prev.front >= 0) {
						// _prevにあって_curに無い物が対象
						m.iterateHistPrev(_prev.front, [&cb](const typename CMGR::Hist& h) {
							// _prevで一度参照されたエントリにはフラグが立ててある -> (nFrame==0)
							if(h.nFrame > 0)
								cb(h);
						});
					}
				} else if(_cur.time == tm-1 && _cur.front >= 0) {
					// 今回何とも衝突しなかった場合
					// _cur のリスト全部対象。ただし履歴はColMgr上では過去のものになってるのでPrevを参照
					m.iterateHistPrev(_cur.front, cb);
				}
			}
			CMask getMask() const { return _mask; }
			spn::Optional<int> getPrevIndex(CTime tm) const {
				if(_cur.time == tm-1 &&
					_cur.front >= 0)
					return _cur.front;
				if(_prev.time == tm-1 &&
					_prev.front >= 0)
					return _prev.front;
				return spn::none;
			}
			//! カレントリスト末尾に新しくインデックスを加える
			/*! \return 直前の先頭インデックス値 */
			int addQueue(CTime tm, int idx) {
				if(tm > _cur.time) {
					_prev = _cur;
					_cur.init(tm);
				}
				return _cur.add(idx);
			}
			UD& refUserData() {
				return _udata;
			}
	};
	template <template <class,class> class BC,
				class CTG,
				class MMGR,
				class IM,
				class UD>
	class ColMgr : public spn::ResMgrA<ColMem<ColMgr<BC,CTG,MMGR,IM,UD>,
											typename MMGR::LHdl,
											typename BC<decltype(std::declval<IM>().im_getBVolume()), IM>::IDType,
											UD>,
										ColMgr<BC,CTG,MMGR,IM,UD>>,
					public Narrow<CTG,IM>
	{
		public:
			using BVolume = decltype(std::declval<IM>().im_getBVolume());
			using user_t = UD;
			using this_t = ColMgr<BC, CTG, MMGR, IM, user_t>;
			using BroadC = BC<BVolume, IM>;
			using NSID = typename BroadC::IDType;
			using HMdl = typename MMGR::SHdl;
			using HLMdl = typename MMGR::LHdl;
			using CMem = ColMem<this_t, HLMdl, NSID, user_t>;
			using base = spn::ResMgrA<CMem, this_t>;
			using HCol = typename base::SHdl;
			using HLCol = typename base::LHdl;
			using IModel = IM;
			// 単方向リスト
			struct Hist {
				HCol	hCol;
				int		nFrame;
				int		nextOffset;
				Hist(HCol hc, int nf): hCol(hc), nFrame(nf), nextOffset(0) {}
			};
		private:
			using narrow = Narrow<CTG,IM>;
			using FrameCount = std::vector<std::pair<HCol, int>>;
			using HistVec = std::vector<Hist>;
			using RemList = std::vector<HLCol>;

			BroadC		_broadC;
			CTime		_accum = 0;
			//! 継続フレーム数を一本の配列に保存
			FrameCount	_fcount;
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

			//! hc0のリストにhc1のエントリが存在するかチェック
			/*!	\return hc1が見つかればその参照を返す */
			spn::Optional<Hist&> _hasPrevCollision(HCol hc0, HCol hc1) {
				auto& c0 = hc0.cref();
				auto pidx = c0.getPrevIndex(_accum);
				if(!pidx)
					return spn::none;
				int prevIdx = *pidx;
				Hist* res = nullptr;
				iterateHistPrev(prevIdx, [hc1,&res](Hist& h){
					if(h.hCol == hc1)
						res = &h;
				});
				if(res)
					return res;
				return spn::none;
			}
			template <class HVec, class CB>
			static void _IterateHist(HVec& hv, int idx, CB cb) {
				AssertP(Trap, idx>=0)
				auto* hst = &hv[idx];
				for(;;) {
					cb(*hst);
					int np = hst->nextOffset;
					if(np != 0)
						hst = reinterpret_cast<decltype(hst)>(reinterpret_cast<uintptr_t>(hst) + np);
					else break;
				}
			}
			HistVec& _switchHist() {
				_doRelease(false);
				_swHist ^= 1;
				_hist[_swHist].clear();
				return _hist[_swHist];
			}
			void _doRelease(bool bBoth) {
				//! 解放を一時的に遅延させていたハンドルを処理
				_bDelete = true;
				_getPrevRM().clear();
				if(bBoth)
					_getCurRM().clear();
				_bDelete = false;
			}
		public:
			ColMgr(): _broadC([this](spn::SHandle sh){
								return HCol(sh).cref().getBVolume(); })
			{
				narrow::Initialize();
			}
			~ColMgr() {
				_doRelease(true);
			}
			CTime getAccum() const {
				return _accum;
			}
			template <class CB>
			void iterateHistCur(int idx, CB cb) const { _IterateHist(_getCurHist(), idx, cb); }
			template <class CB>
			void iterateHistPrev(int idx, CB cb) const { _IterateHist(_getPrevHist(), idx, cb); }
			template <class CB>
			void iterateHistCur(int idx, CB cb) { _IterateHist(_getCurHist(), idx, cb); }
			template <class CB>
			void iterateHistPrev(int idx, CB cb) { _IterateHist(_getPrevHist(), idx, cb); }

			//! 当たり判定対象を追加
			/*! \param ms[in] マスク */
			HLCol addCol(CMask ms, HMdl hm, const UD& ud=UD()) {
				HLCol hlC = base::emplace(ms, hm, ud);
				hlC.ref().setNSID(_broadC.add(hlC, ms));
				return std::move(hlC);
			}
			//! ハンドル解放処理を一時的に遅延させる (デストラクタ時以外)
			void release(HCol hCol) {
				if(!_bDelete && hCol.count() == 1) {
					// 次のupdateまでハンドル解放を遅延
					_getCurRM().push_back(hCol);
					// しかしBroadCからは即刻外す
					_broadC.rem(hCol.cref().getNSID());
				}
				base::release(hCol);
			}
			//! 時間を1進め、衝突履歴の更新
			void update() {
				auto& hist = _switchHist();

				++_accum;
				_broadC.broadCollision([this, &hist](spn::SHandle sh0, spn::SHandle sh1){
					// 同じハンドルという事はあり得ない筈
					AssertP(Trap, sh0!=sh1)

					HCol hc0(sh0), hc1(sh1);
					// 詳細判定
					if(narrow::Hit(hc0.ref().getModel(), hc1.ref().getModel())) {
						auto fn = [this, &hist](HCol hc0, HCol hc1){
							auto &c0 = hc0.ref(),
								&c1 = hc1.ref();
							// 前のフレームでの継続フレームリストを探索
							int nf = 0;
							if(auto hist = _hasPrevCollision(hc0, hc1)) {
								nf = hist->nFrame;
								// 衝突が継続中なので目印としてnFrameを0にしておく
								hist->nFrame = 0;
							}
							int nextID = hist.size();
							int lastID = c0.addQueue(_accum, nextID);
							if(lastID >= 0)
								hist[lastID].nextOffset = (nextID - lastID) * sizeof(Hist);
							hist.push_back(Hist(hc1, nf+1));
						};
						fn(hc0, hc1);
						fn(hc1, hc0);
					}
				});
			}
	};
}
