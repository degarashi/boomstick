#include "geom3D.hpp"

namespace boom {
	namespace geo3d {
		// ------------------------- Idx3 -------------------------
		Idx3::Idx3() {}
		Idx3::Idx3(int i0, int i1, int i2, int pid): _pid(pid) {
			_id[0] = i0;
			_id[1] = i1;
			_id[2] = i2;
		}
		Idx3::Idx3(int i0, int i1, int i2, const Vec3List& vtx, int pid): Idx3(i0,i1,i2,pid) {
			_plane = Plane::FromPts(vtx[i0], vtx[i1], vtx[i2]);
		}
		void Idx3::iterateEdge(std::function<void (int,int,Idx3*)> cb) {
			for(int i=0 ; i<3 ; i++)
				cb(_id[i], _id[(i+1)%3], _neighbor[i]);
		}
		void Idx3::flip() {
			std::swap(_id[0], _id[1]);
		}
		const Plane& Idx3::getPlane() const { return _plane; }
		int Idx3::getID(int n) const { return _id[n]; }
		std::pair<int,int> Idx3::getEdge(int n) const {
			return std::make_pair(int(_id[n]), int(_id[(n+1)%3]));
		}
		void Idx3::rewriteNeighbor(int idx, Idx3* p, bool bOther) {
			if(idx == _id[0]) {
				std::swap(_neighbor[0], p);
				idx = _id[1];
			} else if(idx == _id[1]) {
				std::swap(_neighbor[1], p);
				idx = _id[2];
			} else {
				assert(idx == _id[2]);
				std::swap(_neighbor[2], p);
				idx = _id[0];
			}

			if(bOther)
				p->rewriteNeighbor(idx, this, false);
		}
		void Idx3::initNeighbor(Idx3* p0, Idx3* p1, Idx3* p2) {
			_neighbor[0] = p0;
			_neighbor[1] = p1;
			_neighbor[2] = p2;
		}
		Idx3* Idx3::getNeighbor(int n) { return _neighbor[n]; }
		int Idx3::findEdge(int id) const {
			for(int i=0 ; i<3 ; i++) {
				if(id == _id[i])
					return i;
			}
			return -1;
		}
		// ------------------------- CnvP -------------------------
		CnvP::CnvP(const int (&initial)[4], const Vec3List& vtx): _vtx(vtx), _poly(0x100), _pcCur(0), _pfCur(0) {
			int nV = vtx.size();
			_vcCur = nV-4;
			for(int i=0 ; i<nV ; i++)
				_vcand[i] = i;
			for(int k=0 ; k<4 ; k++) {
				for(int i=0 ; i<nV ; i++) {
					if(_vcand[i] == initial[k]) {
						if(--nV > i)
							_vcand[i] = _vcand[nV];
						break;
					}
				}
			}

			// 最初に使った4頂点は候補から外す
			std::memset(&_useVID[0], 0, sizeof(_useVID));
			const int c_idx[4][3] = {{0,1,2}, {3,1,0}, {2,1,3}, {0,2,3}};
			const int c_neighbor[4][3] = {{1,2,3}, {2,0,3}, {0,1,3}, {0,2,1}};
			for(auto& idx : c_idx)
				_addPoly(initial[idx[0]], initial[idx[1]], initial[idx[2]]);
			for(int i=0 ; i<4 ; i++)
				_pcand[i]->initNeighbor(_pcand[c_neighbor[i][0]].get(), _pcand[c_neighbor[i][1]].get(), _pcand[c_neighbor[i][2]].get());
		}
		void CnvP::_addRefVID(int vid) {
			AssertP(Trap, vid < int(_vtx.size()));
			if(++_useVID[vid] == 1) {
				AssertP(Trap, !_bfVID.check(vid));
				_bfVID.set(vid);
			}
		}
		void CnvP::_unRefVID(int vid) {
			AssertP(Trap, vid < int(_vtx.size()));
			if(--_useVID[vid] == 0) {
				AssertP(Trap, _bfVID.check(vid));
				_bfVID.reset(vid);
			}
		}
		CnvP::PolyF::Ptr& CnvP::_addPoly(int i0, int i1, int i2) {
			auto ptr = _poly.get(i0,i1,i2,_vtx);
			Idx3* p = ptr.get();
			// 頂点番号登録
			_addRefVID(i0);
			_addRefVID(i1);
			_addRefVID(i2);

			// 候補に加える
			return _pcand[_pcCur++] = ptr;
		}
		void CnvP::_remPoly(Idx3* p) {
			for(int i=0 ; i<3 ; i++) {
				// 頂点番号解除
				_unRefVID(p->getID(i));
				// 隣接ポインタの抹消
				p->rewriteNeighbor(p->getID(i), nullptr, true);
			}

			// 候補リストから消す
			int nC = _pcCur;
			for(int i=0 ; i<nC ; i++) {
				if(_pcand[i].get() == p) {
					_poly.put(_pcand[i]);
					// 末尾の要素を突っ込んでおく
					if(--_pcCur > 0)
						_pcand[i] = _pcand[_pcCur];
					return;
				}
			}
			// pfaceからも探す
			int nP = _pfCur;
			for(int i=0 ; i<nP ; i++) {
				if(_pface[i].get() == p) {
					_poly.put(_pface[i]);
					// 末尾の要素を突っ込んでおく
					if(--_pfCur > i)
						_pface[i] = _pface[_pfCur];
					return;
				}
			}
			AssertP(Trap, false);
		}
		Idx3* CnvP::_pickPoly() {
			if(_pcCur == 0)
				return nullptr;
			return _pcand[0].get();
		}
		void CnvP::_remTopPoly() {
			_pface[_pfCur++] = _pcand[0];
			if(--_pcCur > 0)
				_pcand[0] = _pcand[_pcCur];
		}
		void CnvP::_enumEdge(PEVec& dst, IdxP3Set& dstI, const Vec3& vert, Idx3* poly) {
			Idx3* toEnum[3];
			int nE = 0;
			poly->iterateEdge([poly, &vert, &toEnum,&nE,&dst](int i0, int i1, Idx3* neighbor) {
				// 凸物体なので必ず隣接面が存在する
				auto& plane = neighbor->getPlane();
				if(plane.dot(vert) > 0) {
					toEnum[nE++] = neighbor;
				} else {
					// 削除境界を記録
					dst.push_back(PolyEdge(neighbor, i0, i1));
				}
			});
			dstI.insert(poly);
			for(int i=0 ; i<nE ; i++) {
				if(dstI.count(toEnum[i]) == 0)
					_enumEdge(dst, dstI, vert, toEnum[i]);
			}
		}
		bool CnvP::quickHull() {
			// これ以上分割する必要のない場合はポリゴン候補から外す
			Idx3* poly = _pickPoly();
			// 探索ループ
			do {
				// ポリゴンを1つ選ぶ
				// ポリゴン平面より前にある頂点から距離が最大の物を探す
				const auto& plane = poly->getPlane();
				int maxI = 0;
				float maxD = plane.dot(_vtx[_vcand[0]]);
				for(int i=1 ; i<_vcCur ; i++) {
					float d = plane.dot(_vtx[_vcand[i]]);
					if(d > maxD) {
						maxD = d;
						maxI = i;
					}
				}
				if(maxD > 1e-3f) {	// ポリゴン上の頂点が選ばれてしまうので誤差を設ける
					_check();
					// これからポリゴンを張る辺
					// 表を向いているポリゴンを列挙して削除
					PEVec pe;
					IdxP3Set pset;
					_enumEdge(pe, pset, _vtx[_vcand[maxI]], poly);
					for(auto& p : pe) {
						AssertP(Trap, p._srcIdx[0] != _vcand[maxI] && p._srcIdx[1] != _vcand[maxI]);
					}

					// エッジリストをソート
					int nPE = pe.size();
					for(int i=0 ; i<nPE-1 ; i++) {
						int prevI = pe[i]._srcIdx[1];
						for(int j=i+1 ; j<nPE ; j++) {
							if(pe[j]._srcIdx[0] == prevI) {
								if(j > i+1)
									std::swap(pe[i+1], pe[j]);
								prevI = -1;
								break;
							}
						}
						AssertP(Trap, prevI == -1);
					}

					// 古い面を削除
					for(auto& ip : pset)
						_remPoly(ip);
					Idx3 *prevP = nullptr,
						*firstP;
					int prevID;
					AssertP(Trap, pe.size() >= 3 && pset.size() > 0);
					// 新しい面を張る
					for(auto& tp : pe) {
						auto& fr = _addPoly(tp._srcIdx[0], tp._srcIdx[1], _vcand[maxI]);	// 面を追加
						auto* p = fr.get();
						p->rewriteNeighbor(tp._srcIdx[0], tp._dstPoly, false);
						tp._dstPoly->rewriteNeighbor(tp._srcIdx[1], p, false);
						// エッジを繋ぐ (peは順番通り並んでるハズ)
						if(prevP) {
							prevP->rewriteNeighbor(tp._srcIdx[0], p, false);
							p->rewriteNeighbor(_vcand[maxI], prevP, false);
						} else
							firstP = p;
						prevP = p;
						prevID = tp._srcIdx[1];
					}
					firstP->rewriteNeighbor(_vcand[maxI], prevP, false);
					prevP->rewriteNeighbor(prevID, firstP, false);
					_check();

					// 一度選ばれた頂点は絶対に凸包を構成するので候補から外す
					if(--_vcCur > maxI)
						_vcand[maxI] = _vcand[_vcCur];
				} else {
					// これ以上この面について走査する必要は無い
					_remTopPoly();
				}
			} while(_vcCur>0 && (poly=_pickPoly()));
			while(_pcCur > 0)
				_remTopPoly();
			return _bfVID.getNBit() == int(_vtx.size());
		}
		void CnvP::_check() {
			for(int i=0 ; i<_pcCur ; i++) {
				for(int j=0 ; j<3 ; j++) {
					auto* ths = _pcand[i].get();
					auto* p = ths->getNeighbor(j);
					// どのエッジと接続しているか探す
					AssertP(Trap, p->findEdge(ths->getID(j)) >= 0);
				}
			}
		}
		Idx3List CnvP::getResult(Vec3List& dstV) const {
			int nP = _pfCur;
			Idx3List ret(nP);
			Idx3* pDst = &ret[0];

			// qhull前後で頂点数が変わったら有効な頂点だけ残して後は消す
			int nV = _vtx.size();
			int nNV = _bfVID.getNBit();
			if(nNV != nV) {
				AssertP(Trap, nNV < nV);
				// ポリゴンインデックス組換え
				auto tgt = _bfVID.getIndexList<int>();
				std::vector<int> idxcnv(nV);
				// 単純にコピー
				Vec3List vl(nNV);
				for(int i=0 ; i<nNV ; i++) {
					idxcnv[tgt[i]] = i;
					vl[i] = _vtx[tgt[i]];
				}

				int nP = _pfCur;
				for(int i=0 ; i<nP ; i++) {
					auto& p = _pface[i];
					*pDst++ = Idx3(idxcnv[p->getID(0)], idxcnv[p->getID(1)], idxcnv[p->getID(2)], i);
				}
				std::swap(dstV, vl);
			} else {
				// 頂点番号が変わっていないのでベタコピー
				//std::memcpy(&ret[0], &_pface[0], sizeof(Idx3)*nP);
				for(int i=0 ; i<nP ; i++)
					*pDst++ = *_pface[i].get();
			}
			return std::move(ret);
		}
		// ------------------------- ConvexP -------------------------
		ConvexP::ConvexP(): _rflg(~0) {}
		ConvexP::ConvexP(const Vec3* src, int n): ConvexP() {
			_vtx.resize(n);
			for(int i=0 ; i<n ; i++)
				_vtx[i] = src[i];
		}
		ConvexP::ConvexP(Vec3List&& src): ConvexP() {
			_vtx.swap(src);
		}
		ConvexP::ConvexP(ConvexP&& c) {
			swap(c);
		}
		void ConvexP::swap(ConvexP& c) noexcept {
			std::swap(_rflg, c._rflg);
			std::swap(_vGCenter, c._vGCenter);
			std::swap(_vtx, c._vtx);
			std::swap(_pface, c._pface);
		}
		int ConvexP::getNVtx() const { return _vtx.size(); }
		void ConvexP::popVtx() {
			_vtx.pop_back();
			_rflg |= RFL_GCENTER|RFL_POLYFACE;
		}
		Vec3 ConvexP::support(const Vec3& dir) const {
			int nV = getNVtx();
			int idx = 0;
			float maxD = _vtx[0].dot(dir);
			for(int i=1 ; i<nV ; i++) {
				float d = _vtx[i].dot(dir);
				if(d > maxD) {
					maxD = d;
					idx = i;
				}
			}
			return _vtx[idx];
		}
		const Vec3& ConvexP::bs_getGCenter() const {
			if(spn::Bit::ChClear(_rflg, RFL_GCENTER)) {
				int nV = getNVtx();
				Vec3 c(_vtx[0]);
				for(int i=1 ; i<nV ; i++)
					c += _vtx[i];
				_vGCenter = c/nV;
			}
			return _vGCenter;
		}
		const Vec3& ConvexP::bs_getCenter() const { return bs_getGCenter(); }
		Sphere ConvexP::bs_getBVolume() const { INVOKE_ERROR }
		float ConvexP::bs_getArea() const { INVOKE_ERROR }
		Mat33 ConvexP::bs_getInertia() const { INVOKE_ERROR }
		AABB ConvexP::bs_getAABB() const {
			constexpr float fmax = std::numeric_limits<float>::max(),
							fmin = std::numeric_limits<float>::lowest();
			Vec3 vmin(fmax),
				vmax(fmin);
			// 単純に各頂点の最大最小をとる
			for(auto& v : _vtx) {
				vmin.selectMin(v);
				vmax.selectMax(v);
			}
			return AABB(vmin, vmax);
		}

		void ConvexP::addVtx(const Vec3& v) {
			_rflg |= RFL_GCENTER|RFL_POLYFACE;
			_vtx.push_back(v);
		}
		ConvexP ConvexP::operator - (const Vec3& ofs) const {
			return (*this) + (-ofs);
		}
		ConvexP& ConvexP::operator -= (const Vec3& ofs) {
			return (*this) += -ofs;
		}
		ConvexP& ConvexP::operator += (const Vec3& ofs) {
			auto tmp = *this + ofs;
			_vtx.swap(tmp._vtx);
			return *this;
		}
		ConvexP ConvexP::operator + (const Vec3& ofs) const {
			int nV = _vtx.size();
			Vec3List vt(nV);
			for(int i=0 ; i<nV ; i++)
				vt[i] = _vtx[i] + ofs;
			return ConvexP(std::move(vt));
		}
		Vec3 ConvexP::DualTransform(const Plane& plane) {
			float d = plane.d;
			Assert(Trap, plane.d <= 0);
			if(std::fabs(d) < 1e-3f)
				d = (d<0) ? -1e-3f : 1e-3f;

			return plane.getNormal() / -d;
		}
		Plane ConvexP::DualTransform(const Vec3& p) {
			float len = p.length();
			len = 1.0f/len;
			return Plane(p*len, -len);
		}
		Vec3List ConvexP::DualTransform(const PlaneList& plL) {
			Vec3List vlL(plL.size());
			auto* pDst = &vlL[0];
			for(auto& p : plL)
				*pDst++ = DualTransform(p);
			return std::move(vlL);
		}
		Vec3List ConvexP::exportDualTransform() {
			auto& plL = getPolyFace();
			Vec3List vlL(plL.size());
			auto* pDst = &vlL[0];
			for(auto& itr : plL)
				*pDst++ = DualTransform(Plane::FromPts(_vtx[itr.getID(0)], _vtx[itr.getID(1)], _vtx[itr.getID(2)]));
			return std::move(vlL);
		}
		const Vec3& ConvexP::getVtx(int n) const { return _vtx[n]; }
		const Vec3List& ConvexP::getVtxArray() const { return _vtx; }

		bool ConvexP::quickHull() {
// 			mgr_profiler.start("quickHull");
			spn::Bit::Clear(_rflg, RFL_POLYFACE);

			int nV = getNVtx();
			// 頂点が4点なら必然的に凸型である
			if(nV <= 4) {
				_pface.clear();
				_pface.push_back(Idx3(0,1,2, _vtx, 0));
				if(nV == 4) {
					int iFlip = 0;
					// 4番目の頂点で法線を決める
					if(_pface.back().getPlane().dot(_vtx[3]) > 0) {
						_pface.back().flip();
						iFlip = 1;
					}
					const int c_idx[3][2] = {{1,0}, {0,2}, {2,1}};
					for(auto& idx : c_idx)
						_pface.push_back(Idx3(3, idx[iFlip], idx[iFlip^1], _vtx, 0));
				}
// 				mgr_profiler.end("quickHull");
				return false;
			}

			// 頂点候補リストを初期化
			std::vector<int> cand(nV);
			for(int i=0 ; i<nV ; i++)
				cand[i] = i;

			// 最初の4面体を形成
			int fsmp[4];
			// X軸について最大最小の頂点を通過する線分と、それについて距離が最大の頂点を探す => 3角形の定義
			auto itrC = cand.begin();
			float minX = _vtx[*itrC].x,
				maxX = minX;
			decltype(itrC) minI = itrC,
							maxI = minI;
			while(++itrC != cand.end()) {
				float x = _vtx[*itrC].x;
				if(minX > x) {
					minX = x;
					minI = itrC;
				} else if(maxX < x) {
					maxX = x;
					maxI = itrC;
				}
			}
			int minID = *minI,
				maxID = *maxI;
			fsmp[0] = minID;
			fsmp[1] = maxID;
			// 構成するのに使った頂点は候補から外す
			if(maxID > minID)
				--maxID;

			cand.erase(cand.begin()+minID);
			cand.erase(cand.begin()+maxID);

			Line lineX(_vtx[fsmp[0]], Vec3(_vtx[fsmp[1]]-_vtx[fsmp[0]]).normalization());
			itrC = cand.begin();
			float maxD = lineX.dist_sq(_vtx[*itrC]);
			decltype(itrC) maxLI = itrC;
			while(++itrC != cand.end()) {
				float d = lineX.dist_sq(_vtx[*itrC]);
				if(d > maxD) {
					maxD = d;
					maxLI = itrC;
				}
			}
			fsmp[2] = *maxLI;
			cand.erase(maxLI);

			// 3角形の片方の面法線で一番遠いものを選ぶ => 4面体を形成
			// 全て裏面にある事も考えられるので最大と最小を記録
			Plane plane = Plane::FromPts(_vtx[fsmp[0]], _vtx[fsmp[1]], _vtx[fsmp[2]]);
			itrC = cand.begin();
			maxD = plane.dot(_vtx[*itrC]);
			maxLI = itrC;
			decltype(itrC) minLI = itrC;
			float minD = maxD;
			while(++itrC != cand.end()) {
				float d = plane.dot(_vtx[*itrC]);
				if(d < minD) {
					minD = d;
					minLI = itrC;
				} else if(d > maxD) {
					maxD = d;
					maxLI = itrC;
				}
			}
			// ポリゴンの裏側に何か点があればそれを採用
			if(minD < 0) {
				fsmp[3] = *minLI;
				cand.erase(minLI);
			} else {
				// でなければポリゴンをフリップして表側のを採用
				fsmp[3] = *maxLI;
				std::swap(fsmp[0], fsmp[1]);
				cand.erase(maxLI);
			}

// 			mgr_profiler.start("qhIter");
			// 4面体をCnvPに入力
			CnvP cnvp(fsmp, _vtx);
			cnvp.quickHull();
			auto res = cnvp.getResult(_vtx);
			_pface.swap(res);
// 			mgr_profiler.end("qhIter");
// 			mgr_profiler.end("quickHull");
			return res.size() != _pface.size();
		}
		bool ConvexP::hasPoint(const Vec3& p) {
			const Idx3List& idxL = getPolyFace();
			for(auto& idx : idxL) {
				if(Plane::FromPts(_vtx[idx.getID(0)], _vtx[idx.getID(1)], _vtx[idx.getID(2)]).dot(p) > 1e-5f)
					return false;
			}
			return true;
		}
		void ConvexP::dualTransform(Vec3List& dst, const Vec3& dir, const Vec3& pos) {
			const Idx3List& idxL = getPolyFace();
			for(auto& p : idxL) {
				Plane plane = p.getPlane();
				if(plane.getNormal().dot(dir) > 0)
					continue;

				plane.d = std::min(plane.dot(pos), -5e-5f);
				dst.push_back(DualTransform(plane));
			}
		}
		const Idx3List& ConvexP::getPolyFace() {
			if(spn::Bit::ChClear(_rflg, RFL_POLYFACE))
				quickHull();
			return _pface;
		}
	}
}
