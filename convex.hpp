#pragma once
#include "spinner/misc.hpp"
#include "geom3D.hpp"

namespace boom {
	namespace geo3d {
		//! 三角形ポリゴン
		class Polygon {
			//! 更新フラグ
			enum RFLG {
				NORMALFLAG = 0x01,
				PLANEFLAG = 0x02,
				CENTERFLAG = 0x04,
				ALLFLAG = 0xff
			};
			private:
				Vec3		_vtx[3];				//! 頂点3つ
				uint32_t	_ud;					//! ユーザー任意データ
				// キャッシュ変数
				mutable Vec3		_vNormal, _vCenter;
				mutable Plane		_plane;
				mutable uint32_t	_rflg;

				Vec3 _calcLineHit(const Vec3& pos, const Vec3& dir) const;
				static bool _IsValidRange(const Vec3& res);

			public:
				Polygon();
				Polygon(const Vec3& v0, const Vec3& v1, const Vec3& v2);
				Polygon(const Vec3* vsrc);

				void init(const Vec3& v0, const Vec3& v1, const Vec3& v2);
				void init(const Vec3* vsrc);
				void init(const Plane& plane);

				const Vec3& getNormal() const;
				const Vec3& getVtx(int n) const;
				const Plane& getPlane() const;
				Plane getEdgePlane(int n) const;
				Segment getEdge(int n) const;
				void setVtx(int n, const Vec3& src);
				float calcRadius() const;
				float calcAcreage() const;

				// ユーザーデータ関係
				void setUserData(uint32_t dat);
				uint32_t getUserData() const;

				bool hit(const Ray& ray) const;
				bool hit(const Line& ls) const;
				std::pair<Vec3,bool> nearest(const Vec3& p) const;
				bool isOnTriangleSpace(const Vec3& p) const;
				//! 平面でポリゴンを分割し、結果を引数に返す
				/*! 簡単の為にConvex<Vec3>で処理してから結果をPolygonに変換
					自身は表側に位置するが他は保障されない
					@param[out] dst 出力先配列
					@param[in] plane 分割する平面
					@return 分割された数(表側, 裏側) */
				std::pair<int,int> split(Polygon (&dst)[3], const Plane& plane) const;

				Polygon operator * (const AMat43& m) const;
				const Vec3& getGCenter3D() const;
				Vec3 support(const Vec3& dir) const;
		};

		//! ポリゴンと平面の関係 (for Point)
		enum class PsType {
			Front,
			Back,
			OnPlane,
			Bridge,
			Invalid
		};
		//! ポリゴンと平面の関係 (for Line)
		enum class LsType {
			Front,
			Back,
			Bridge,
			Invalid
		};
		struct ConvexUD_Col {
			uint32_t		id;
		};
		//! 平面上の凸包
		/*! 様々な頂点形式を持たせるためにテンプレート実装
			\param VType 	頂点型
			\param UD		ユーザーデータ型 */
		template <class VT, class UD=spn::none_t>
		class Convex {
			using UDType = UD;
			using VType = VT;
			using VList = std::vector<VType>;
			using FList = std::vector<float>;
			using Type = Convex<VType, UDType>;
			enum RFLG {
				RFLG_CENTER = 0x01,			//!< 重心点
				RFLG_NORMAL = 0x02,			//!< 法線
				RFLG_DEGENERATE = 0x04,		//!< 頂点重複チェック
				RFLG_CPOINT = 0x08,			//!< AABBの中心
				RFLG_PLANE = 0x10,
				RFLG_ALL = 0x1f
			};
			private:
				VList	_vtx;
				UDType	_ud;

				// キャッシュ
				mutable Vec3		_vCenter, _vNormal, _vCPoint;
				mutable Plane		_plane;
				mutable uint32_t	_rflg;
			public:
				Convex();
				Convex(int n);
				Convex(const Convex& c);
				Convex(Convex&& c);

				void init();
				//! 頂点数を変更(メモリ領域を予約, またはシュリンク)
				void setNVtx(int nv);
				//! ダブった頂点を省略
				void degeneration();
				bool checkValid();				//!< degene & linear
				bool checkConvex() const;		//!< 凸ポリゴンになっているか判定
				bool checkLinear() const;		//!< 全ての頂点が同一平面上にあるか判定
				bool isOnPlane(const Plane& p, float th) const;
				//! 全ての頂点に対しての内積を取得
				FList checkPlaneFlags(const Plane& plane) const;
				//! 指定の面に頂点を寄せる
				void adjustToPlane(const Plane& plane);

				UDType& refUserData();
				const UDType& getUserData() const;

				void setVertexArray(const VType* src, int nv);
				void splitThis(const Plane& plane, Convex& bDst);	//!< 前面のポリゴンがこのポリゴンに適用される
				void splitThis(const Plane& plane);					//!< 背面は捨てる
				void split(const Plane& plane, Convex& fDst, Convex& bDst) const;

				//! ポリゴンが引数の平面のどちら側にあるか
				PsType detectPlaneSide(const Plane& plane) const;
				constexpr static float CREATE_PLANE_DIST = 1;
				//! ポリゴン生成関数
				static Convex FromPlane(const Plane& plane, float dist=CREATE_PLANE_DIST);
				static Convex FromVtx(const VType& v0, const VType& v1, const VType& v2, const VType& v3);
				static Convex FromVtx(const VType& v0, const VType& v1, const VType& v2);

				void addVtx(const VType& v);		//!< 末尾に頂点を足す
				void popVtx();						//!< 末尾の頂点を取り除く

				// ---- Getter functions ----
				int getNVtx() const;
				const Vec3& getCPoint() const;
				const Plane& getPlane() const;
				const Vec3& getCenter() const;
				const Vec3& getNormal() const;
				Vec3 getSumNormal() const;			//!< (頂点が平面上でなくても作用)

				const VType& getVtx(int n) const;
				const Vec3& getPos(int n) const;
				VType& refVtx(int n);
				Vec3& refPos(int n);

				int getNPoly() const;				//!< 3角形ポリゴンで何枚分か
				Polygon getPolygon(int n) const;	//!< 3角形ポリゴンを取り出す
				void extractPolygons(VList& vl, IndexList& il) const;

				Convex& operator = (const Convex& c);
				Vec3x2 support(const Vec3& vc) const;	//!< サポート射像
				Plane getEdgePlane(int eid) const;
				Vec3 getEdgeNormal(int eid) const;
				bool selectSuitableOrigin();		//!< TriangleFanにする時に一番無理のない頂点を始点に選択
				void selectOrigin(int idx);			//!< 始点を選択
				void invert();						//!< 頂点を逆順にする
				Convex loopExtract(int vtxI) const;	//!< 指定の頂点を始点とした凸ポリゴン

				std::pair<VList, IndexList> subdivide(float dot_area) const;

				//! ポリゴンボリュームを計算(中心座標から一番遠い頂点の距離が半径)
				float calcRadius() const;
				Sphere calcSphere() const;

				//! 行列変換
				Convex& operator *= (const AMat43& m);
				Convex operator * (const AMat43& m) const;
				//! 最小座標，最大座標計算
				Vec3x2 calcMinMax() const;
				//! 行列変換した後の最小座標，最大座標計算
				Vec3x2 calcMinMax(const AMat43& m) const {
					Convex c = *this * m;
					return c.calcMinMax();
				}
				//! 歪みのない座標系に変換した後の最小座標，最大座標
				Vec3x2 calcMinMax(const Vec3& xAxis, const Vec3& yAxis, const Vec3& zAxis) const;

				void addOfsVec(const Vec3& ofs);		//!< オフセットを加える
				Convex cloneReverse() const;			//!< ポリゴンを裏返して出力
				spn::Optional<Vec3> lineCollision(const Vec3& vBeg, const Vec3& vEnd, float offset=0, bool planeThres=true) const;	//!< 線分との当たり判定
				bool hit(const Sphere& s) const;		//!< バウンディングボリュームとの当たり判定

				void swap(Convex& c) throw();
		};
		using ColCv = Convex<Vec3, ConvexUD_Col>;
	}
}
