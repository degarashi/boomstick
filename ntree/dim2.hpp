#pragma once
#include "ntree.hpp"

namespace boom {
	namespace geo2d {
		struct AABB;
	}
	namespace ntree {
		// -------- <<コリジョンツリー次元ポリシークラス>> --------
		//! 4分木用
		class CTDim_2D {
			public:
				struct Id {
					union {
						struct {
							unsigned int v0 : 16;
							unsigned int v1 : 16;
						};
						uint32_t value;
					};
					Id() = default;
					Id(uint32_t val);
					Id(uint32_t val0, uint32_t val1);
					Id operator + (Id id) const;
					Id operator - (Id id) const;
					Id operator / (int s) const;
				};
				using Id_OP = spn::Optional<Id>;

				constexpr static int N_Dim = 2,
									N_LayerSize = spn::NPow<2,N_Dim>::result;
				using BVolume = ::boom::geo2d::AABB;

				//! geo2d::IModelからモートンIdを算出
				/*! \return [0]=Min(XY)のモートンId
							[1]=Max(XY)のモートンId
							[2]=Min(XY)軸インデックス
							[3]=Max(XY)軸インデックス */
				static std::tuple<MortonId,MortonId,Id,Id> ToMortonId(const BVolume& objBox, int nwidth, float unit, float ofs);
				//! 軸別インデックス値からモートンIDに変換
				static MortonId ToMortonId(uint32_t x, uint32_t y);
				//! モートンIDから軸別インデックス値へ変換
				static Id ToXYPos(MortonId id);
				//! オブジェクトの軸分類
				/*! \param[out] dst			オブジェクトの振り分け先 (4マス) 重複あり
					\param[in] src			振り分けるオブジェクト配列
					\param[in] nSrc			srcの要素数
					\param[in] centerId		振り分け中心点(軸別のマス目インデックス) */
				static void Classify(const VolEntry** (&dst)[N_LayerSize], const VolEntry* src, int nSrc, Id centerId);
				using CBClassify = std::function<void (const VolEntry&, int)>;
				static int Classify(const VolEntry& ve, Id centerId, const CBClassify& cb);
				//! 1つ下の階層の中心インデックス値を算出
				/*! \param[out] idDst	1つ下の階層の中心座標インデックス
					\param[in] id		現在の階層の中心座標インデックス
					\param[in] width	現在の階層のマス幅(width*2で全体) */
				static void CalcLowerCenterId(Id (&idDst)[N_LayerSize], Id id, int width);
		};
	}
}
