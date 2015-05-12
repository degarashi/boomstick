#pragma once
#include "collision.hpp"

namespace boom {
	//! broad-phase collision manager (round robin)
	/*!	属性フラグが0x80000000の時はBへ、それ以外はAに登録
		A->A, A->Bでは判定が行われるが B->Bはされない
		\tparam BV	バウンディングボリューム
	*/
	template <class BV>
	class BroadC_RoundRobin {
		public:
			enum Type {
				TypeA,
				TypeB,
				NumType
			};
			using NS_id = uint32_t;
			//! 外部(CollisionMgr)から内部Nodeを区別するId
			struct IDType {
				NS_id	id;
				Type	typ;
			};
			using IDTypeV = std::vector<IDType>;
			using BVolume = BV;
		private:
			struct Node {
				CMask			mask;	//!< 登録時に渡されたマスク値
				spn::SHandle	hCol;	//!< コリジョンマネージャで使われるリソースハンドル
				BVolume			volume;
			};
			using Nodes = spn::noseq_list<Node, std::allocator, NS_id>;
			Nodes			_node[NumType];
			using FGetBV = std::function<BVolume (spn::SHandle)>;
			const FGetBV 	_fGetBV;

			//! バウンディングボリュームの更新
			/*! \param[in] idv 非nullptrならリストに指定したオブジェクトのみの更新 */
			void _refreshBV(const IDTypeV* idv) {
				if(idv) {
					// 指定したオブジェクトのみを対象
					for(auto& idt : *idv) {
						Node& node = _node[idt.typ].get(idt.id);
						node.volume = _fGetBV(node.hCol);
					}
				} else {
					// 全てのオブジェクトが対象
					for(int i=0 ; i<NumType ; i++) {
						for(auto& nd : _node[i])
							nd.volume = _fGetBV(nd.hCol);
					}
				}
			}
			template <class CB>
			static int _Proc(const Node& nd0, const Node& nd1, CB&& cb) {
				// 属性マスクによる判定
				if(nd0.mask & nd1.mask) {
					// 境界ボリュームチェック
					if(nd0.volume.hit(nd1.volume)) {
						std::forward<CB>(cb)(nd0.hCol, nd1.hCol);
						return 1;
					}
				}
				return 0;
			}
			static Type _DetectType(CMask m) {
				return (m & 0x80000000) ? TypeB : TypeA;
			}
		public:
			BroadC_RoundRobin(FGetBV cb, float /*fieldSize*/, float /*fieldOfs*/): _fGetBV(cb) {}
			IDType add(spn::SHandle sh, CMask mask) {
				Type typ = _DetectType(mask);
				return IDType{_node[typ].add(Node{mask, sh}), typ};
			}
			void rem(const IDType& idt) {
				_node[idt.typ].rem(idt.id);
			}
			//! リストに登録してある全ての物体に対して衝突判定
			/*!	\param[in] mask		コリジョンマスク値
				\param[in] bv		判定対象のバウンディングボリューム
				\param[in] cb		コールバック関数(spn::SHandle) */
			template <class CB>
			void checkCollision(CMask mask, const BVolume& bv, CB&& cb) {
				_refreshBV(nullptr);
				auto fnChk = [mask, &bv, cb=std::forward<CB>(cb)](auto& nd) {
					for(auto& obj : nd) {
						if(mask & obj.mask) {
							if(bv.hit(obj.volume))
								cb(obj.hCol);
						}
					}
				};
				if(_DetectType(mask) == TypeA) {
					// TypeBと判定
					fnChk(_node[TypeB]);
				}
				// TypeAと判定
				fnChk(_node[TypeA]);
			}
			//! リストに溜め込まずに直接コールバックを呼ぶ
			/*!	\param[in] ac_t		累積時間
				\param[in] cb		コールバック関数(SHandle,SHandle)
				\param[in] idv		非nullptr時はリストに登録してあるオブジェクトのみBVの更新を行う
				\return 衝突したペア数 */
			template <class CB>
			int broadCollision(CB&& cb, const IDTypeV* idv=nullptr) {
				_refreshBV(idv);

				int count = 0;
				auto itrB_a = _node[TypeA].begin(),
					itrE_a = _node[TypeA].end(),
					itrE_a1 = std::prev(itrE_a, 1);

				int nA = _node[TypeA].size();
				if(nA > 0) {
					// A -> A
					for(auto itr=itrB_a ; itr!=itrE_a1 ; ++itr) {
						const auto& nodeA = *itr;
						for(auto itr2=itr+1 ; itr2!=itrE_a ; ++itr2)
							count += _Proc(nodeA, *itr2, std::forward<CB>(cb));
					}
					// A -> B
					if(!_node[TypeB].empty()) {
						auto itrB_b = _node[TypeB].begin(),
						itrE_b = _node[TypeB].end();

						for(auto itr=itrB_a ; itr!=itrE_a ; ++itr) {
							const auto& nodeA = *itr;
							for(auto itr2=itrB_b ; itr2!=itrE_b ; ++itr2)
								count += _Proc(nodeA, *itr2, std::forward<CB>(cb));
						}
					}
				}
				return count;
			}
	};
}
