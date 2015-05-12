#include "dim2.hpp"
#include "../geom2D.hpp"
#include "../bits.hpp"

namespace boom {
	namespace ntree {
		// ---------------------- CTDim_2D ----------------------
		CTDim_2D::Id::Id(uint32_t val):
			v0(val),
			v1(val)
		{}
		CTDim_2D::Id::Id(uint32_t val0, uint32_t val1) {
			v0 = val0;
			v1 = val1;
		}
		CTDim_2D::Id CTDim_2D::Id::operator + (Id id) const {
			return Id(v0 + id.v0,
						v1 + id.v1);
		}
		CTDim_2D::Id CTDim_2D::Id::operator - (Id id) const {
			return Id(v0 - id.v0,
						v1 - id.v1);
		}
		CTDim_2D::Id CTDim_2D::Id::operator / (int s) const {
			return Id(v0 / s,
						v1 / s);
		}
		// ---------------------- CTDim_2D ----------------------
		std::tuple<MortonId,MortonId,CTDim_2D::Id,CTDim_2D::Id> CTDim_2D::ToMortonId(const BVolume& objBox, int nwidth, float unit, float ofs) {
			auto sat = [nwidth,unit,ofs](float w) -> uint32_t {
				return spn::Saturate(static_cast<int>((w+ofs)*unit), 0, nwidth-1); };

			auto& ab = objBox;
			uint32_t tmin[2] = {sat(ab.minV.x), sat(ab.minV.y)},
					tmax[2] =  {sat(ab.maxV.x), sat(ab.maxV.y)};
			return std::make_tuple(ToMortonId(tmin[0], tmin[1]),
									ToMortonId(tmax[0], tmax[1]),
									Id(tmin[0], tmin[1]),
									Id(tmax[0], tmax[1]));
		}
		MortonId CTDim_2D::ToMortonId(uint32_t x, uint32_t y) {
			return SeparateBits1(x) | (SeparateBits1(y) << 1);
		}
		CTDim_2D::Id CTDim_2D::ToXYPos(MortonId id) {
			return Id(ComposeBits1(id & 0x55555555),
						ComposeBits1((id & 0xaaaaaaaa) >> 1));
		}
		int CTDim_2D::Classify(const VolEntry& ve, Id centerId, const CBClassify& cb) {
			const uint32_t cx = centerId.v0,
							cy = centerId.v1;
			Id pmin = *reinterpret_cast<const Id*>(&ve.posMin),
				pmax = *reinterpret_cast<const Id*>(&ve.posMax);
			AssertP(Trap, pmin.v0 <= pmax.v0 && pmin.v1 <= pmax.v1)
			int count = 0;
			if(pmin.v0 < cx) {
				// 左リストに登録
				if(pmin.v1 < cy) {
					// 左上
					cb(ve, 0);
					++count;
				}
				if(pmax.v1 >= cy) {
					// 左下
					cb(ve, 2);
					++count;
				}
			}
			if(pmax.v0 >= cx) {
				// 右リストに登録
				if(pmin.v1 < cy) {
					// 右上
					cb(ve, 1);
					++count;
				}
				if(pmax.v1 >= cy) {
					// 右下
					cb(ve, 3);
					++count;
				}
			}
			return count;
		}
		void CTDim_2D::Classify(const VolEntry** (&dst)[N_LayerSize], const VolEntry* src, int nSrc, Id centerId) {
			auto cb =[&dst](const VolEntry& ve, int idx){
							*dst[idx]++ = &ve;
						};
			// オブジェクト振り分け
			for(int i=0 ; i<nSrc ; i++) {
				int count = Classify(src[i], centerId, cb);
				AssertP(Trap, count > 0)
			}
		}
		void CTDim_2D::CalcLowerCenterId(Id (&idDst)[N_LayerSize], Id id, int width) {
			uint32_t cx = id.v0,
					cy = id.v1;
			// 中心座標計算
			idDst[0] = Id(cx - width, cy - width);
			idDst[1] = Id(cx + width, cy - width);
			idDst[2] = Id(cx - width, cy + width);
			idDst[3] = Id(cx + width, cy + width);
		}
	}
}
