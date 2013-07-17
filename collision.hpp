/*! 2D, 3Dを問わないコリジョンに関する定義 */
#include <vector>
#include <algorithm>
#include <cassert>
#include <boost/optional.hpp>
#include "spinner/layerbit.hpp"
#include "spinner/noseq.hpp"

namespace boom {
	//! ColResultで使用: ダミー情報クラス
	class DummyInfo {
		public:
			using type = void;
			void resize(int /*nMaxObj*/) {}
			void clear() {}
	};

	//! ColResultで使用: 値の平均を計算する
	template <class T>
	class AverageEntry {
		T		_value;
		int		_count;

		public:
			using type = T;
			AverageEntry(): _value(0), _count(0) {}
			void operator()(const T& val, int sign) {
				_value += val * sign;
				++_count;
			}
			T get() const {
				return _value * spn::_sseRcp22Bit(_count);
			}
	};

	//! ColResultで使用: 値を合算する
	template <class T, class ID>
	class SharedEntry {
		public:
			using type = T;
			using id_type = ID;
		private:
			using EntList = std::vector<T>;
			EntList		_ent;
		public:
			void resize(int nMaxObj) {
				static_assert(std::is_integral<ID>::value, "");
				_ent.resize(nMaxObj);
			}
			void clear() {
				_ent.clear();
			}
			void pushInfo(id_type id0, id_type id1, const typename T::type& t) {
				assert(id0 < _ent.size() && id1 < _ent.size());
				_ent[id0](t, 1);
				_ent[id1](t, -1);
			}
			boost::optional<decltype(_ent[0].get())> getInfo(id_type id0, id_type /*id1*/) const {
				return _ent[id0].get();
			}
	};

	//! ColResultで使用: オブジェクトの組み合わせ1組に対して1つの情報を割り当てる
	template <class T>
	class SharedInfo {
		using type = T;
		using id_type = typename T::id_type;
		using InfoVec = std::vector<T>;
		struct Link {
			uint16_t minID,
					maxID;
			Link() = default;
			Link(uint16_t v): minID(v), maxID(v) {}
		};
		using LinkVec = std::vector<Link>;
		using Cursor = std::pair<uint16_t,uint16_t>;
		using CurVec = std::vector<Cursor>;

		InfoVec	_infoVec;
		LinkVec	_linkVec;
		CurVec	_curVec;
		id_type	_curID,
				_nItem;

		public:
			void resize(int nMaxObj) {
				_curVec.resize(nMaxObj);
				_curID = _nItem = 0;
				std::memset(&_curVec[0], 0, sizeof(uint16_t)*nMaxObj*2);
			}
			void clear() {
				_infoVec.clear();
				_linkVec.clear();
				_curVec.clear();
			}
			// 組み合わせは常にid0<id1で、値が前後したりは考慮しない
			template <class P>
			void pushInfo(id_type id0, id_type id1, P&& t) {
				assert(id0 < id1);
				assert(id0 < _curVec.size());

				_infoVec.emplace_back(t);
				if(_curID != id0) {
					// ID0エントリの切り替え
					if(_nItem > 0) {
						_curVec[_curID] = _curVec[id0] = _linkVec.size();
						_linkVec.emplace_back(id1);
					} else if(!_linkVec.empty()) {
						auto& b = _linkVec.back();
						b.minID = b.maxID = id1;
					}
					_curID = id0;
					_nItem = 0;
				} else {
					// 同一ID0の中でID1の間隔が空いた場合、リンクを分割
					auto& link = _linkVec.back();
					if(link.maxID+1 != id1) {
						_linkVec.emplace_back(id1);
					} else
						++link.maxID;
				}
			}
			//! 組み合わせに関連付けられた情報を取り出す
			boost::optional<const T&> getInfo(id_type id0, id_type id1) const {
				if(id0 < id1)
					std::swap(id0,id1);

				const Cursor& cur = _curVec[id0];
				if(cur.first == cur.second)
					return boost::none;

				uint32_t idx = cur.first,
						ofs = 0;
				while(id1 > _linkVec[idx].maxID) {
					ofs += _linkVec[idx].minID;
					if(++idx >= cur.second)
						return boost::none;
				}
				ofs += _linkVec[idx].minID;
				return _linkVec[ofs];
			}
	};

	//! コリジョン判定の結果を参照しやすい形で格納
	/*! テンプレート引数による固定数実装。適当にクラスを複数用意するなり工夫する */
	template <int MAXN, class INFO=DummyInfo>
	class ColResult {
		public:
			using id_type = typename INFO::id_type;
		private:
			// ------------ 1次変数 ------------
			using LBits = spn::LayerBitArray<MAXN>;
			using LBVec = std::vector<LBits>;
			LBVec	_lbVec;

			// ------------ 2次変数 ------------
			struct Item {
				uint16_t	objID,
							infoID;
			};
			using ItemArray = std::vector<Item>;
			using CursorVec = std::vector<std::pair<uint16_t,uint16_t>>;

			ItemArray		_array;
			CursorVec		_cursor;
			INFO			_info;

		public:
			ColResult() = default;
			ColResult(ColResult&& cr):
				_lbVec(std::forward<LBVec>(cr._lbVec)),
				_array(std::forward<ItemArray>(cr._array)),
				_cursor(std::forward<CursorVec>(cr._cursor)),
				_info(std::forward<INFO>(cr._info))
			{}
			ColResult& operator = (ColResult&& cr) {
				std::swap(_lbVec, cr._lbVec);
				std::swap(_array, cr._array);
				std::swap(_cursor, cr._cursor);
				std::swap(_info, cr._info);
				return *this;
			}
			ColResult& operator = (const ColResult& cr) = default;

			//! オブジェクト数を指定して初期化
			void setNumObjects(int n) {
				_lbVec.resize(n);
				_info.resize(n);
			}
			//! 内部情報をクリア(実際にはメモリ解放されないかもしれない)
			void clear() {
				_lbVec.clear();
				_array.clear();
				_cursor.clear();
				_info.clear();
			}
			//! 共有情報を付けて登録
			/*! id0 < id1 で、必ず小さい番号順で呼ばれる前提 */
			template <class INFOT>
			void setFlagWithInfo(id_type id0, id_type id1, const INFOT& info) {
				setFlag(id0,id1);
				_info.pushInfo(id0, id1, info);
			}
			//! ヒットフラグを立てる
			void setFlag(id_type id0, id_type id1) {
				_lbVec[id0].set(id1);
				_lbVec[id1].set(id0);
			}
			//! あるIDとヒットしているIDを列挙
			template <class CB>
			void iterate(id_type id0, CB cb) const {
				_lbVec[id0].iterateBit(cb);
			}
			//! iterateと同じだが、こちらは共有情報も取得する
			/*! Callback: void operator()(ID id1, boost::optional<const INFO::type&> t);
				よって共有情報が無い場合もありうる */
			template <class CB>
			void iterateWithInfo(id_type id0, CB cb) {
				auto acb = [id0,&cb,this](int id1) {
					cb(id1, _info.getInfo(id0,id1));
				};
				_lbVec[id0].iterate(acb);
			}
			const INFO& getInfo() const {
				return _info;
			}
	};

	//! 剛体シミュレーションで使用する係数
	struct RCoeff {
		float	spring,	//!< スプリング係数
				dumper,	//!< ダンパ係数
				fricD,	//!< 動摩擦係数
				fricS,	//!< 静止摩擦係数
				fricMS;	//!< 最大静止摩擦係数
	};

	//! broad-phase collision manager (round robin)
	/*! A->A, A->Bでは判定が行われるが B->Bはされない */
	template <class NODE>
	class BroadC_RoundRobin {
		public:
			enum TYPE {
				TYPE_A,
				TYPE_B,
				NUM_TYPE
			};
		private:
			using Nodes = spn::noseq_list<NODE>;
			Nodes	_node[NUM_TYPE];

		public:
			using const_iterator = typename Nodes::const_iterator;
			using id_type = typename Nodes::id_type;
			int getNum(TYPE typ) const { return _node[typ].size(); }
			int getNumObj() const { return getNum(TYPE_A) + getNum(TYPE_B); }

			id_type add(TYPE typ, NODE obj) {
				return _node[typ].add(obj);
			}
			void rem(TYPE typ, id_type id) {
				_node[typ].rem(id);
			}
			const_iterator cbegin(TYPE typ) const {
				return _node[typ].cbegin();
			}
			const_iterator cend(TYPE typ) const {
				return _node[typ].cend();
			}
			template <class CRes>
			int collision(CRes& cr) const {
				int nA = _node[TYPE_A].size(),
					nB = _node[TYPE_B].size(),
					count = 0;
				typename CRes::narrow_type nchk;
				cr.clear();
				cr.setNumObjects(nA + nB);
				// A -> A
				auto itrB_a = _node[TYPE_A].begin(),
					itrE_a = _node[TYPE_A].end(),
					itrE_a1 = itrE_a;
				--itrE_a1;

				int iCur=0;
				for(auto itr=itrB_a ; itr!=itrE_a1 ; ++itr,++iCur) {
					const auto& nodeA = *(*itr).cref().get();

					auto itr2 = itr;
					++itr2;
					for(int jCur=iCur+1; itr2!=itrE_a ; ++itr2,++jCur) {
						const auto& nodeA2 = *(*itr2).cref().get();
						if(nchk(nodeA, nodeA2)) {
							cr(iCur, jCur, nodeA, nodeA2, nchk.getInfo());
							++count;
						}
					}
				}
				// A -> B
				if(!_node[TYPE_B].empty()) {
					auto itrB_b = _node[TYPE_B].begin(),
							itrE_b = _node[TYPE_B].end(),
							itrE_b1 = itrE_b;
					--itrE_b1;

					iCur=0;
					for(auto itr=itrB_a ; itr!=itrE_a ; ++itr,++iCur) {
						const auto& nodeA = *(*itr).cref().get();
						int jCur=0;
						for(auto itr2=itrB_b ; itr2!=itrE_b ; ++itr2,++jCur) {
							const auto& nodeB = *(*itr2).cref().get();
							if(nchk(nodeA,nodeB)) {
								cr(iCur, jCur, nodeA, nodeB, nchk.getInfo());
								++count;
							}
						}
					}
				}
				return count;
			}
	};
}
