#include "test2D.hpp"

namespace boom {
	namespace test2d {
		MonoPolygon::MonoPolygon(DistV&& dL, DistV&& dR, const Vec2& ori, const Vec2& dir, float wofs):
			_vOrigin(ori),
			_vDir(dir),
			_distL(std::move(dL)),
			_distR(std::move(dR)),
			_widthOffset(wofs)
		{}
		MonoPolygon MonoPolygon::Random(const FRandF& rff, const FRandI& rfi, const spn::RangeF& rV, const spn::RangeF& rLen, int nV) {
			nV = std::max(3, nV);
			// ランダムにスイープ方向を決める
			auto dir = Vec2::RandomDir(rff);
			// 原点座標
			auto origin = Vec2::Random(rff, rV);
			// 長さ
			float length = rff(rLen);

			int nLeft = rfi({1, nV-2}),
				nRight = rfi({0, nV-2-nLeft});
			auto fnMakeLengthList = [&rff](const int n, const float len) {
				std::vector<Vec2>	distL(n+1);
				float sumL = 0;
				for(auto& d : distL) {
					float num = rff({1e-2f, 1e1f});
					sumL += num;
					d.x = num;
					d.y = rff({0, 1e1f});
				}
				distL.back().y = 0;
				for(auto& d : distL) {
					d.x /= sumL;
					d.x *= len;
				}
				return distL;
			};
			float width_offset = rff(rLen * 1e-1f);
			auto distL = fnMakeLengthList(nLeft, length),
				distR = fnMakeLengthList(nRight, length);
			return MonoPolygon(std::move(distL),
								std::move(distR),
								origin,
								dir,
								width_offset);
		}
		geo2d::PointL MonoPolygon::getPoints() const {
			int nLeft = _distL.size()-1,
				nRight = _distR.size()-1;
			geo2d::PointL pts(2 + nLeft + nRight);
			auto* ptr = pts.data();
			spn::Vec2 dir90{-_vDir.y, _vDir.x};
			// 始点を追加
			*ptr++ = _vOrigin + dir90*_widthOffset;
			float cur = 0;
			// 左側にランダムな数の頂点(最低1)
			for(int i=0 ; i<nLeft ; i++) {
				auto& dist = _distL[i];
				cur += dist.x;
				*ptr++ = _vOrigin + _vDir*cur + dir90*(dist.y + _widthOffset);
			}
			cur += _distL.back().x;
			// 終点に頂点を配置
			*ptr++ = _vOrigin + _vDir*cur + dir90*_widthOffset;
			// 右側にランダムな数の頂点
			cur -= _distR.back().x;
			for(int i=nRight-1 ; i>=0 ; i--) {
				auto& dist = _distR[i];
				*ptr++ = _vOrigin + _vDir*cur - dir90*(dist.y - _widthOffset);
				cur -= dist.x;
			}
			Assert(Trap, pts.data()+pts.size() == ptr)
			return pts;
		}
		const spn::Vec2& MonoPolygon::getDir() const {
			return _vDir;
		}
		bool MonoPolygon::hit(const spn::Vec2& p, float threshold) const {
			Vec2 dir90(-_vDir.y, _vDir.x);
			auto toP = p - (_vOrigin + dir90*_widthOffset);
			float d_vert = _vDir.dot(toP),
				d_horz = std::sqrt(toP.len_sq() - spn::Square(d_vert));
			if(d_vert < 0)
				return false;
			auto fnGetHDist = [](const auto& distV, float d_vert, float invalid){
				int nL = distV.size();
				float cur_y = 0;
				for(int i=0 ; i<nL ; i++) {
					auto &dist = distV[i];
					if(d_vert <= dist.x)
						return spn::Lerp(cur_y, dist.y, d_vert / dist.x);
					d_vert -= dist.x;
					cur_y = dist.y;
				}
				return invalid;
			};
			if(dir90.dot(toP) < 0)
				d_horz *= -1;
			float dL = fnGetHDist(_distL, d_vert, -1.f),
					dR = fnGetHDist(_distR, d_vert, -1.f);
			return spn::IsInRange(d_horz, -dR-threshold, dL+threshold);
		}
	}
}
