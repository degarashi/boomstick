#pragma once
#include "collision.hpp"
#include "spinner/resmgr.hpp"

namespace boom {
	//! broad-phase collision manager (round robin)
	/*!	属性フラグが0x80000000の時はBへ、それ以外はAに登録
		A->A, A->Bでは判定が行われるが B->Bはされない */
	template <class BV, class IM>
	class BroadC_RoundRobin {
		public:
			enum Type {
				TypeA,
				TypeB,
				NumType
			};
			struct IDType {
				uint32_t	id;
				Type		typ;
			};
		private:
			using IModel = IM;
			using BVolume = BV;
			using FGetBV = std::function<BVolume (spn::SHandle)>;
			struct Node {
				uint32_t		mask;
				spn::SHandle	hCol;
				BVolume			volume;
			};
			using Nodes = spn::noseq_list<Node, uint32_t>;
			Nodes			_node[NumType];
			const FGetBV 	_fGetBV;

			void _refreshBV() {
				for(int i=0 ; i<NumType ; i++) {
					for(auto& nd : _node[i])
						nd.volume = _fGetBV(nd.hCol);
				}
			}
			template <class CB>
			static int _Proc(const Node& nd0, const Node& nd1, CB cb) {
				// 属性マスクによる判定
				if(nd0.mask & nd1.mask) {
					// 境界ボリュームチェック
					if(nd0.volume.hit(nd1.volume)) {
						cb(nd0.hCol, nd1.hCol);
						return 1;
					}
				}
				return 0;
			}
		public:
			BroadC_RoundRobin(FGetBV cb): _fGetBV(cb) {}
			IDType add(spn::SHandle sh, uint32_t mask) {
				Type typ = (mask & 0x80000000) ? TypeB : TypeA;
				return IDType{_node[typ].add(Node{mask, sh}), typ};
			}
			void rem(const IDType& idt) {
				_node[idt.typ].rem(idt.id);
			}
			//! リストに溜め込まずに直接コールバックを呼ぶ
			/*!	\param[in] ac_t		累積時間
				\param[in] cb		コールバック関数(SHandle,SHandle) */
			template <class CB>
			int broadCollision(CB cb) {
				_refreshBV();

				int count = 0;
				auto itrB_a = _node[TypeA].begin(),
					itrE_a = _node[TypeA].end(),
					itrE_a1 = itrE_a;
				--itrE_a1;

				int nA = _node[TypeA].size();
				if(nA > 0) {
					// A -> A
					for(auto itr=itrB_a ; itr!=itrE_a1 ; ++itr) {
						const auto& nodeA = *itr;
						for(auto itr2=itr+1 ; itr2!=itrE_a ; ++itr2)
							count += _Proc(nodeA, *itr2, cb);
					}
					// A -> B
					if(!_node[TypeB].empty()) {
						auto itrB_b = _node[TypeB].begin(),
						itrE_b = _node[TypeB].end();

						for(auto itr=itrB_a ; itr!=itrE_a ; ++itr) {
							const auto& nodeA = *itr;
							for(auto itr2=itrB_b ; itr2!=itrE_b ; ++itr2)
								count += _Proc(nodeA, *itr2, cb);
						}
					}
				}
				return count;
			}
	};
}
