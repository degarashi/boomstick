/*! 2D, 3Dを問わないコリジョンに関する定義 */
#include <vector>
#include <algorithm>
#include <cassert>
#include <boost/optional.hpp>
#include "spinner/layerbit.hpp"
#include "spinner/noseq.hpp"

namespace boom {
	namespace c_ent {
		//! c_info::Pairsで使用: 値の平均を計算する
		template <class T>
		class Average {
			T		_value;
			int		_count;

			public:
				using type = T;
				Average(): _value(0), _count(0) {}
				void operator ()(const T& val) {
					_value += val;
					++_count;
				}
				T get() const {
					return _value * spn::_sseRcp22Bit(_count);
				}
		};
		//! c_info::Pairsで使用: 値の合計を出す
		template <class T>
		class Sum {
			T		_value;
			public:
				using type = T;
				Sum(): _value(0) {}
				void operator ()(const T& val) { _value += val; }
				const T& get() const { return _value; }
		};
	}
	//! マイナス演算子で負数ではなく別の値を返す
	template <class T>
	struct DualValue {
		T	_value[2];

		T& operator() (int n) { return _value[n]; }
		operator const T& () const { return _value[0]; }
		const T& operator - () const { return _value[1]; }
	};
	namespace c_info {
		//! ColResultで使用: ダミー情報クラス
		/*! 変数をセットされても何もしない */
		class Dummy {
			public:
				using type = void;
				void resize(int /*nMaxObj*/) {}
				void clear() {}
				template <class ID, class DAT>
				void pushInfo(ID /*id0*/, ID /*id1*/, const DAT& /*dat*/) {}
				template <class ID>
				int getInfo(ID /*id0*/, ID /*id1*/) const { throw std::runtime_error("not supported"); }
		};
		//! ColResultで使用: オブジェクト毎に値を合算する
		template <class T, class ID>
		class Pairs {
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
				template <class T2>
				void pushInfo(id_type id0, id_type id1, const T2& t) {
					assert(id0 < _ent.size() && id1 < _ent.size());
					_ent[id0](t);
					_ent[id1](-t);
				}
				boost::optional<decltype(_ent[0].get())> getInfo(id_type id0, id_type /*id1*/) const {
					return _ent[id0].get();
				}
		};
		//! ColResultで使用: オブジェクトの組み合わせ1組に対して1つの情報を割り当てる
		template <class T>
		class Shared {
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

					_infoVec.emplace_back(std::forward<P>(t));
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
	}

	//! コリジョン判定の結果を参照しやすい形で格納
	/*! テンプレート引数による固定数実装。適当にクラスを複数用意するなり工夫する */
	template <int MAXN, class INFO=c_info::Dummy>
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
	template <class NODE, class CAST=NODE>
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
			using iterator = typename Nodes::iterator;
			using id_type = typename Nodes::id_type;
			int getNum(TYPE typ) const { return _node[typ].size(); }
			int getNumObj() const { return getNum(TYPE_A) + getNum(TYPE_B); }

			id_type add(TYPE typ, NODE obj) {
				return _node[typ].add(obj);
			}
			void rem(TYPE typ, id_type id) {
				_node[typ].rem(id);
			}
			const_iterator cbegin(TYPE typ) const { return _node[typ].cbegin(); }
			const_iterator cend(TYPE typ) const { return _node[typ].cend(); }
			iterator begin(TYPE typ) { return _node[typ].begin(); }
			iterator end(TYPE typ) { return _node[typ].end(); }

			template <class CB>
			void iterate(CB cb) {
				for(int i=0 ; i<NUM_TYPE ; i++) {
					auto typ = static_cast<TYPE>(i);
					for(auto itr=begin(typ) ; itr!=end(typ) ; itr++)
						cb(*itr);
				}
			}

			//! リストに溜め込まずに直接コールバックを呼ぶ
			/*! \param[in] cr		コールバック関数
				\param[in] nchk		詳細判定関数 */
			template <class CRes, class NChk>
			int broadCollision(CRes& cr, NChk nchk) {
				int count = 0;
				auto itrB_a = _node[TYPE_A].begin(),
					itrE_a = _node[TYPE_A].end(),
					itrE_a1 = itrE_a;
				--itrE_a1;

				if(!_node[TYPE_A].empty()) {
					// A -> A
					int iCur=0;
					for(auto itr=itrB_a ; itr!=itrE_a1 ; ++itr,++iCur) {
						const auto& nodeA = static_cast<CAST>(*itr);

						auto itr2 = itr;
						++itr2;
						for(int jCur=iCur+1; itr2!=itrE_a ; ++itr2,++jCur) {
							const auto& nodeA2 = static_cast<CAST>(*itr2);
							if(nchk(nodeA, nodeA2)) {
								cr(iCur, jCur, *itr, *itr2, nchk.getInfo());
								++count;
							}
						}
					}
				}
				// A -> B
				if(!_node[TYPE_B].empty()) {
					auto itrB_b = _node[TYPE_B].begin(),
							itrE_b = _node[TYPE_B].end(),
							itrE_b1 = itrE_b;
					--itrE_b1;

					int nA = _node[TYPE_A].size();
					int iCur=0;
					for(auto itr=itrB_a ; itr!=itrE_a ; ++itr,++iCur) {
						const auto& nodeA = static_cast<CAST>(*itr);
						int jCur=0;
						for(auto itr2=itrB_b ; itr2!=itrE_b ; ++itr2,++jCur) {
							const auto& nodeB = static_cast<CAST>(*itr2);
							if(nchk(nodeA, nodeB)) {
								cr(iCur, jCur, *itr, *itr2, nchk.getInfo());
								++count;
							}
						}
					}
				}
				return count;
			}

			//! 判定結果を一時的に貯めておく用
			template <class Aux>
			struct AuxInfo {
				int			id0, id1;
				const CAST *pR0, *pR1;
				Aux			aux;

				template <class A>
				AuxInfo(int i0, int i1, const CAST* r0, const CAST* r1, A&& a):
					id0(i0), id1(i1), pR0(r0), pR1(r1), aux(std::forward<A>(a)) {}
			};
			template <class Aux>
			using CPList = std::vector<AuxInfo<Aux>>;
			template <class Aux>
			struct AuxCB {
				CPList<Aux>	buff;
				template <class A>
				void operator() (int id0, int id1, CAST& c0, CAST& c1, A&& a) {
					buff.emplace_back(id0, id1, &c0, &c1, std::forward<A>(a));
				}
			};

			//! 判定結果をリスト出力する
			template <class NChk>
			CPList<typename NChk::aux_type> broadCollision(NChk nchk) {
				AuxCB<typename NChk::aux_type> cb;
				broadCollision(cb, nchk);
				return std::move(cb.buff);
			}
	};
}
