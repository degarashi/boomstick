#include "geom2D.hpp"

namespace boom {
	namespace geo2d {
		Convex::Convex(const PointL& pl): point(pl) {}
		Convex::Convex(PointL&& pl): point(pl) {}
		Convex::Convex(std::initializer_list<Vec2> v): point(v.size()) {
			auto itrD = point.begin();
			auto itr = v.begin();
			while(itr != v.end())
				*itrD++ = *itr++;
		}
		Convex::Convex(Convex&& c): point(std::forward<PointL>(c.point)) {}
		Convex& Convex::operator = (Convex&& c) {
			point = std::move(c.point);
			return *this;
		}
		Vec2 Convex::bs_getCenter() const {
			return std::get<2>(area_inertia_center());
		}
		float Convex::bs_getInertia() const {
			return std::get<1>(area_inertia_center());
		}

		namespace {
			bool Back(const Vec2& v0, const Vec2& v1, const Vec2& v2) {
				Vec2 to0(v0 - v1),
					to2(v2 - v1);
				float c = to0.cw(to2),
					d = to0.dot(to2);
				return c >= 0 || (c > -1e-3f && d > 1e-3f);
			}
		}
		Convex Convex::FromConcave(const PointL& src) {
			int nV = src.size();
			Assert(Trap, nV >= 3)

			// X軸についてソート
			PointL tsrc(src);
			std::sort(tsrc.begin(), tsrc.end(), [](const Vec2& v0, const Vec2& v1){ return v0.x < v1.x; });

			PointL pts(nV+1);
			Vec2* pDst = &pts[0];
			*pDst++ = tsrc[0];
			*pDst++ = tsrc[1];
			for(int rc=2 ; rc<nV ; rc++) {
				while(Back(pDst[-2], pDst[-1], tsrc[rc])) {
					if(--pDst == &pts[1])
						break;
				}
				*pDst++ = tsrc[rc];
			}
			*pDst++ = tsrc[nV-2];
			Vec2* pTgt = pDst-1;
			for(int rc=nV-3 ; rc>=0 ; rc--) {
				while(Back(pDst[-2], pDst[-1], tsrc[rc])) {
					if(--pDst == pTgt)
						break;
				}
				*pDst++ = tsrc[rc];
			}
			// 末尾がダブる為、削る
			--pDst;

			Assert(Trap, pDst <= &pts[0]+nV)
			pts.resize(pDst - &pts[0]);
			return Convex(std::move(pts));
		}

		float Convex::bs_getArea() const {
			AreaSum as;
			iterate(as);
			return as.result;
		}
		void Convex::addOffset(const Vec2& ofs) {
			for(auto& p : point)
				p += ofs;
		}
		Vec2 Convex::support(const Vec2& dir) const {
			Vec2 result = point[0];
			float dMax = point[0].dot(dir);
			int nV = point.size();
			for(int i=1 ; i<nV ; i++) {
				float d = point[i].dot(dir);
				if(d > dMax) {
					result = point[i];
					dMax = d;
				}
			}
			return result;
		}
		Segment Convex::getOuterSegment(int n) const {
			return Segment(point[(n+1)%point.size()], point[n]);
		}
		Line Convex::getOuterLine(int n) const {
			return Line(point[n], (point[spn::CndSub(n+1, point.size())] - point[n]).normalization());
		}

		Convex Convex::operator * (const AMat32& m) const {
			int nP = point.size();
			PointL pl(nP);
			for(int i=0 ; i<nP ; i++)
				pl[i] = point[i].asVec3(1) * m;
			return std::move(pl);
		}
		Convex& Convex::operator *= (const AMat32& m) {
			for(auto& p : point)
				p = p.asVec3(1) * m;
			return *this;
		}
		Convex& Convex::operator += (const Vec2& ofs) {
			for(auto& p : point)
				p += ofs;
			return *this;
		}
		void Convex::distend(float width, float mindist) {
			int nP = point.size();
			PointL tmp(nP);
			// 単純に重心からの方向ベクトルで算出
			auto center = bs_getCenter();
			for(int i=0 ; i<nP ; i++) {
				auto dir = point[i] - center;
				float dist = dir.normalize();
				dist = std::max(mindist, dist+width);
				tmp[i] = center + dir*dist;
			}
			*this = FromConcave(tmp);
		}

		bool Convex::checkCW() const {
			int nV = point.size();
			if(nV < 3)
				return false;
			Vec2 c = (point[0] + point[1] + point[2]) * (1.f/3);
			for(int i=0 ; i<nV-1 ; i++) {
				if((point[i+1] - point[i]).cw(c-point[i]) < 0)
					return false;
			}
			return (point[0] - point[nV-1]).cw(c-point[0]) >= 0;
		}
		void Convex::adjustLoop() {
			// point全部が凸包に使用される頂点と仮定
			Vec2 c = (point[0] + point[1] + point[2]) * (1.f/3);
			int nV = point.size();
			using AngPair = std::pair<spn::RadF, Vec2>;
			std::vector<AngPair> tmp(nV);
			for(int i=0 ; i<nV ; i++)
				tmp[i] = AngPair{spn::AngleValue(point[i]-c), point[i]};
			std::sort(tmp.begin(), tmp.end(), [](const AngPair& a0, const AngPair& a1) { return a0.first < a1.first; });

			for(int i=0 ; i<nV ; i++)
				point[i] = tmp[i].second;
		}

		Convex Convex::GetOverlappingConvex(const Convex& m0, const Convex& m1, const Vec2& inner) {
			Assert(Trap, m0.hit(inner) && m1.hit(inner), "invalid inner point")

			// DualTransformで直線から点に変換
			int nV0 = m0.getNPoints(),
				nV1 = m1.getNPoints();
			std::vector<Vec2> vTmp(std::max(nV0,nV1));
			std::vector<Vec2> v2(nV0+nV1);
			for(int i=0 ; i<nV0 ; i++)
				vTmp[i] = m0.getPoint(i) - inner;
			for(int i=0 ; i<nV0 ; i++)
				v2[i] = Dual2(vTmp[i], vTmp[spn::CndSub(i+1, nV0)]);

			for(int i=0 ; i<nV1 ; i++)
				vTmp[i] = m1.getPoint(i) - inner;
			for(int i=0 ; i<nV1 ; i++)
				v2[i+nV0] = Dual2(vTmp[i], vTmp[spn::CndSub(i+1, nV1)]);

			// 凸包を求める
			Convex cc = Convex::FromConcave(v2);
			// DualTransformで直線から点に変換
			int nVD = cc.point.size();
			v2.resize(nVD);
			for(int i=0 ; i<nVD ; i++)
				v2[i] = -Dual2(cc.point[i], cc.point[spn::CndSub(i+1, nVD)]) + inner;

			return std::move(v2);
		}
		int Convex::getNPoints() const { return point.size(); }
		Vec2 Convex::getPoint(int n) const { return point[n]; }

		namespace {
			int BinSearch(const PointL& point, int nV, const Convex& mdl, int a, int b, bool flip) {
				for(;;) {
					if(a+1 >= b)
						break;
					int c = (a+b)/2;
					if(mdl.hit(point[spn::CndRange(c, nV)]) ^ flip)
						b = c;
					else
						a = c;
				}
				return spn::CndRange(a, nV);
			}
		}
		std::pair<bool,PointL> Convex::getOverlappingPoints(const Convex& mdl, const Vec2& inner) const {
			auto res = checkPosition(inner);
			if(res.first != ConvexPos::Outer) {
				int nV = point.size();
				int a,b, begI;
				// 2分探索で衝突が始まる地点を探す
				if(mdl.hit(point[res.second])) {
					a = res.second - nV;
					b = res.second;
					// 衝突開始インデックス(これ自体は衝突していない)
					begI = BinSearch(point, nV, mdl, a,b,false);
				} else {
					begI = res.second;
					if(!mdl.hit(point[spn::CndSub(res.second+1, nV)]))
						return std::make_pair(false, PointL());
				}

				// 衝突が終わる地点を探す
				a = begI + 1;
				b = a + nV;
				int endI = BinSearch(point,nV,mdl,a,b,true);		// 衝突終了インデックス(これ自体衝突している)

				if(begI == endI) {
					// 全部出力
					return std::make_pair(true, point);
				}
				endI = spn::CndSub(endI+1, nV);
				endI = spn::CndAdd(endI, begI, nV)+1;
				PointL pts(endI - begI);
				auto* pDst = &pts[0];
				// 開始地点から終了地点までを出力
				while(begI != endI) {
					*pDst++ = point[spn::CndSub(begI,nV)];
					++begI;
				}
				return std::make_pair(false, std::move(pts));
			}
			return std::make_pair(false, PointL());
		}

		std::pair<ConvexPos,int> Convex::checkPosition(const Vec2& pos, float threshold) const {
			// 適当に内部点を算出
			Vec2 inner = (point[0] + point[1] + point[2]) * (1.f/3);
			Vec2 toP(pos - inner);
			if(toP.len_sq() < NEAR_THRESHOLD_SQ) {
				// 重心がちょうどposと重なってしまったら少しずらす
				inner.lerp(point[0], 0.5f);
			}

			// 内部のどの三角形に該当するか2分探索
			int nV = point.size();
			struct Tmp {
				GAP_VECTOR(vec, 2,
						   (int index)
						   (float cw))

				Tmp& operator = (const Tmp& t) {
					reg_store_ps(vec.m, reg_load_ps(t.vec.m));
					return *this;
				}
			};
			// tmpA = point[0]へのベクトル
			Tmp tmpA, tmpB, tmpC;
			tmpA.vec = point[0]-inner;
			tmpA.index = 0;
			tmpA.cw = tmpA.vec.cw(toP);
			// tmpB = point[nV]へのベクトル
			tmpB.index = nV;
			tmpB.vec = tmpA.vec;
			tmpB.cw = tmpA.cw;

			while(tmpA.index+1 < tmpB.index) {
				// 捜査対象範囲の半分で区切る => tmpC
				tmpC.index = (tmpA.index+tmpB.index) / 2;
				tmpC.vec = point[tmpC.index] - inner;
				tmpC.cw = tmpC.vec.cw(toP);
				float crAC = tmpA.vec.cw(tmpC.vec);

				if(tmpA.cw >= 0) {
					// posはtmpAの右側
					if(crAC >= 0) {
						// tmpCはtmpAの右側
						if(tmpC.cw <= 0) {
							// posはtmpCの左側
							tmpB = tmpC;
						} else
							tmpA = tmpC;
					} else
						tmpB = tmpC;
				} else {
					// posはtmpAの左側
					if(crAC >= 0) {
						// tmpCはtmpAの右側
						tmpA = tmpC;
					} else {
						if(tmpC.cw >= 0) {
							// posはtmpCの右側
							tmpA = tmpC;
						} else
							tmpB = tmpC;
					}
				}
			}
			tmpB.index = spn::CndSub(tmpB.index, nV);
			float d = (point[tmpB.index] - point[tmpA.index]).cw(pos - point[tmpA.index]);
			ConvexPos posf;
			if(d > threshold)
				posf = ConvexPos::Inner;
			else if(d < -threshold)
				posf = ConvexPos::Outer;
			else
				posf = ConvexPos::OnLine;
			return std::make_pair(posf, tmpA.index);
		}
		bool Convex::hit(const Vec2& p, float t) const {
			return checkPosition(p, t).first != ConvexPos::Outer;
		}
		std::tuple<float,float,Vec2> Convex::area_inertia_center() const {
			int nL = point.size();
			AreaList al(nL);
			iterate(std::ref(al));
			float invarea = spn::Rcp22Bit(al.sum);
			float areaInv3 = invarea * (1.0f/3),
					areaInv6 = invarea * (1.0f/6);
			Vec2 center(0,0);
			float inertia = 0;
			iterate([&, areaInv3, areaInv6](int n, const Vec2& p0, const Vec2& p1) {
				center += (p0 + p1) * al.areaL[n] * areaInv3;
				inertia += al.areaL[n] * areaInv6 * (p0.dot(p0) + p0.dot(p1) + p1.dot(p1));
			});
			inertia -= center.len_sq();
			return std::make_tuple(al.sum, inertia, center);
		}
		Convex2 Convex::splitTwo(const Line& ls) const {
			int nV = point.size();
			PointL pt0(nV*2), pt1(nV*2);		// 最初に最大容量確保しておき、後で縮める
			Vec2* pDst[2] = {&pt0[0], &pt1[0]};

			constexpr float DOT_THRESHOLD = 1e-5f;
			// 右側は0, ライン上は前回値, 左側は1を返す
			auto fcheck = [&ls](int prev, const Vec2& pt) -> int {
				float c = ls.dir.ccw(pt - ls.pos);
				if(c >= DOT_THRESHOLD)
					return 0;
				if(c <= -DOT_THRESHOLD)
					return 1;
				return prev;
			};
			// 走査の始点はライン上でない頂点とする
			int cur = 0;
			int flag;
			while(cur != nV) {
				const auto& p = point[cur++];
				flag = fcheck(-1, p);
				if(flag != -1)
					break;
				// とりあえず片方のリストに貯めておく
				*pDst[0]++ = p;
			}
			// 2点以上がライン上なら分割なし
			if(cur > 2) {
				if(flag == 0x01) {
					// すべて右側
					return Convex2({}, point);
				} else {
					// すべて左側
					return Convex2(point, {});
				}
			}
			if(flag == 0x01) {
				std::swap(pt0, pt1);
				std::swap(pDst[0], pDst[1]);
			}
			int prev = flag;

			auto fadd = [&pDst, &ls](const Vec2& pPrev, const Vec2& pCur, int flg) {
				switch(flg) {
					case 0x03:		// Left -> Left
						*pDst[1]++ = pPrev;
						break;
					case 0x02: {	// Left -> Right
						auto res = Segment(pPrev, pCur).crossPoint(ls);
						*pDst[1]++ = pPrev;
						if(res.second != LinePos::OnLine)
							res.first = pCur;
						if(pPrev.dist_sq(res.first) >= 1e-6f)
							*pDst[1]++ = res.first;
						*pDst[0]++ = res.first;
						break; }
					case 0x01: {	// Right -> Left
						auto res = Segment(pPrev, pCur).crossPoint(ls);
						if(res.second != LinePos::OnLine)
							res.first = pCur;
						*pDst[0]++ = pPrev;
						if(pPrev.dist_sq(res.first) >= 1e-6f)
							*pDst[0]++ = res.first;
						*pDst[1]++ = res.first;
						break; }
					case 0x00:		// Right -> Right
						*pDst[0]++ = pPrev;
						break;
				}
			};
			for(int i=cur ; i<nV ; i++) {
				// 正数が左側、負数は右側
				const auto& p = point[i];
				prev = fcheck(prev, p);
				flag = ((flag<<1) | prev) & 0x03;
				fadd(point[i-1], p, flag);
			}
			prev = fcheck(prev, point[0]);
			flag = ((flag<<1) | prev) & 0x03;
			fadd(point[nV-1], point[0], flag);

			if(pDst[0] - &pt0[0] < 2)
				pDst[0] = &pt0[0];
			else if(pDst[1] - &pt1[0] < 2)
				pDst[1] = &pt1[0];
			else {
				if(cur == 2) {
					// 最初の1点がライン上だったので頂点を複製
					if((flag & 1) == 1)
						*pDst[1]++ = pt0[0];
					else
						*pDst[0]++ = pt1[0];
				}
			}

			pt0.resize(pDst[0] - &pt0[0]);
			pt1.resize(pDst[1] - &pt1[0]);
			return Convex2(std::move(pt0), std::move(pt1));
		}
		Convex Convex::split(const Line& ls) {
			auto res = splitTwo(ls);
			std::swap(res.first, *this);
			return std::move(res.second);
		}
		void Convex::splitThis(const Line& l) {
			split(l);
		}
		Circle Convex::bs_getBVolume() const {
			// 多分遅いアルゴリズムだが、今はこれで我慢
			// 全ての3点の組み合わせを調べる
			int nV = point.size();
			AssertP(Trap, nV >= 3)
			Circle c;
			c.vCenter = (point[0] +  point[1])/2;
			c.fRadius = point[0].distance(point[1])/2;
			for(int i=0 ; i<nV-2 ; i++) {
				for(int j=i+1 ; j<nV-1 ; j++) {
					for(int k=j+1 ; k<nV ; k++) {
						PolyM p(point[i], point[j], point[k]);
						c.appendBoundary(&p);
					}
				}
			}
			return c;
		}
		std::tuple<bool,Vec2,Vec2> Convex::checkCrossingLine(const Line& ls) const {
			//TODO: 効率の良い算出方法を探す
			Vec2 pt[2];
			Vec2* ppt = pt;
			int nV = point.size();
			for(int i=0 ; i<nV ; i++) {
				auto res = Segment(point[i], point[spn::CndSub(i+1, nV)]).crossPoint(ls);
				if(res.second != LinePos::NotHit) {
					*ppt++ = res.first;
					if(ppt == pt+2)
						break;
				}
			}
			return std::make_tuple(ppt == pt+2, pt[0], pt[1]);
		}
		Vec3 Dual(const Plane& plane) {
			float d = plane.d;
			AssertP(Trap, plane.d <= 0)
			if(std::fabs(d) < 1e-3f)
				d = (d<0) ? -1e-3f : 1e-3f;
			return plane.getNormal() / -d;
		}
		Plane Dual(const Vec3& p) {
			float invlen = spn::Rcp22Bit(p.length());
			return Plane(p*invlen, -invlen);
		}
		namespace {
			constexpr float MIN_DUAL2DIST = 1e-5f;
		}
		Vec2 Dual2(const Vec2& v0, const Vec2& v1) {
			Vec2 dir = (v1-v0).normalization();
			float dist = dir.cw(-v0);
			dist = std::max(MIN_DUAL2DIST, dist);
			return dir * spn::Rcp22Bit(dist);
		}
		Line Dual(const Vec2& v) {
			Vec2 t0(0,-v.y),
				t1(1, v.x-v.y);
			t1 = (t1-t0).normalization();
			return Line(t0, t1);
		}
		Vec2 Dual(const Line& ls) {
			const Vec2& t0 = ls.pos;
			Vec2 t1 = ls.dir + t0;
			return Vec2((t0.y-t1.y) / (t0.x-t1.x),
						-t0.ccw(t1) / (t1.x-t0.x));
		}

		std::ostream& operator << (std::ostream& os, const Convex& c) {
			os << "Convex(2d) [";
			int idx=0;
			for(auto& p : c.point)
				os << idx++ << ": " << p << std::endl;
			return os << ']';
		}
	}
}
