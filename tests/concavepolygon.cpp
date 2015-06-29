#include "test2D.hpp"

namespace boom {
	namespace test2d {
		geo2d::PointL ConcavePolygon::getPoints() const {
			// 始点を探す -> どれが一辺でも隣接ポリゴンを持ってない物
			int polyId = -1,
				edgeId;
			int nPoly = _poly.size();
			for(int i=0 ; i<nPoly ; i++) {
				auto* p = _poly[i];
				for(int j=0 ; j<3 ; j++) {
					if(!p->data[j].pNeighbor) {
						polyId = i;
						edgeId = j;
						break;
					}
				}
				if(polyId >= 0)
					break;
			}
			Assert(Trap, polyId >= 0)
			int count = 0;
			for(int i=0 ; i<nPoly ; i++) {
				for(auto& pe : _poly[i]->data) {
					if(!pe.pNeighbor)
						++count;
				}
			}

			// 時計回りに外周をたどってリスト出力
			std::vector<int> index;
			index.reserve(_vtx.size());
			auto* pFirst = _poly[polyId];
			index.push_back(pFirst->index[edgeId]);
			// 次のポリゴンへ
			auto* pe = pFirst->data.data() + ((edgeId+1)%3);
			if(!pe->pNeighbor) {
				index.push_back(pFirst->index[(edgeId+1)%3]);
				pe = pFirst->data.data() + ((edgeId+2)%3);
			}
			auto* p = pe->pNeighbor;
			edgeId = pe->edgeId;
			while(p != pFirst) {
				// 隣にポリゴンを持っていたらそのまま進む
				auto* pe = p->data.data() + ((edgeId+1)%3);
				if(pe->pNeighbor) {
					p = pe->pNeighbor;
					edgeId = pe->edgeId;
					continue;
				}
				// なければ最初の頂点を追加した後に次へ
				index.push_back(p->index[(edgeId+1)%3]);
				pe = p->data.data() + ((edgeId+2)%3);
				if(!pe->pNeighbor) {
					index.push_back(p->index[(edgeId+2)%3]);
					pe = p->data.data() + edgeId;
				}
				p = pe->pNeighbor;
				edgeId = pe->edgeId;
			}
			Assert(Trap, index.size() == _vtx.size())

			int sz = _vtx.size();
			geo2d::PointL pl(sz);
			for(int i=0 ; i<sz ; i++)
				pl[i] = _vtx[index[i]];
			return std::move(pl);
		}
		bool ConcavePolygon::hit(const spn::Vec2& p, float threshold) const {
			// 中身は単なる三角形ポリゴンのリストなので個別に判定
			for(auto& pl : _poly) {
				if(geo2d::Poly(_vtx[pl->index[0]],
								_vtx[pl->index[1]],
								_vtx[pl->index[2]]).hit(p, threshold))
					return true;
			}
			return false;
		}
		ConcavePolygon ConcavePolygon::Random(const FRandF& rff, const FRandI& rfi, const spn::RangeF& rV, int nPoly) {
			// ランダムに3角ポリゴンを繋げていく
			nPoly = std::max(1, nPoly);
			auto fnRV = [&](){ return Vec2::Random(rff, rV); };

			ConcavePolygon res;
			auto& vtx = res._vtx;
			auto& poly = res._poly;
			// 最初のポリゴンを生成
			for(;;) {
				vtx = {fnRV(), fnRV(), fnRV()};
				float dist = std::min(vtx[0].distance(vtx[1]), vtx[0].distance(vtx[2]));
				dist = std::min(dist, vtx[1].distance(vtx[2]));
				if(dist > 1e1f)
					break;
			}
			int idx[3] = {0,1,2};
			if(!geo2d::Poly(vtx[0], vtx[1], vtx[2]).isCW())
				std::swap(idx[1], idx[2]);
			poly.emplace_back(new IdxT({{idx[0], idx[1], idx[2]}}, {}));
			--nPoly;

			while(nPoly > 0) {
				// ランダムでポリゴンを選び、隣に何も繋がっていない辺を選ぶ
				int polyId,
					edgeId = -1;
				for(;;) {
					polyId = rfi({0, int(poly.size()-1)});
					auto* p = poly[polyId];
					for(int i=0 ; i<3 ; i++) {
						if(!p->data[i].pNeighbor) {
							edgeId = i;
							break;
						}
					}
					if(edgeId >= 0)
						break;
				}
				auto* p = poly[polyId];
				geo2d::Poly tpoly(vtx[p->index[0]], vtx[p->index[1]], vtx[p->index[2]]);
				auto &v0 = tpoly.point[edgeId],
					&v1 = tpoly.point[(edgeId+1)%3];
				auto dir = (v1 - v0).normalization();
				dir = Vec2(-dir.y, dir.x);
				auto edgeCenter = (v0 + v1)/2;

				// 外側に変異
				// 変異量を特に制限しない
				float moveDist = rff({1e-2f, 1e1f});
				auto newVtx = edgeCenter + dir * moveDist;
				// 一旦ポリゴンを作り…
				geo2d::Poly tpoly0(vtx[p->index[edgeId]], newVtx, vtx[p->index[(edgeId+1)%3]]);
				Assert(Trap, tpoly0.isCW())
				// 他のポリゴンと重なっていたらやり直し
				bool bResult = false;
				for(auto* p : poly) {
					geo2d::Poly tpoly(vtx[p->index[0]], vtx[p->index[1]], vtx[p->index[2]]);
					Assert(Trap, tpoly.isCW())
					if((bResult = tpoly.hit(tpoly0)))
						break;
				}
				if(bResult)
					continue;

				// 新しい頂点を登録
				int newVid = vtx.size();
				vtx.push_back(newVtx);
				// ポリゴンを形成
				poly.push_back(new IdxT({{p->index[edgeId],
										newVid,
										p->index[(edgeId+1)%3]}},
										{{Neighbor<IdxT>{nullptr, 0}, {nullptr,0}, {p, edgeId}}}));
				auto& pe = p->data[edgeId];
				pe.pNeighbor = poly.back();
				pe.edgeId = 2;
				--nPoly;
			}
			return std::move(res);
		}
	}
}
