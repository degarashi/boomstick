//! 3D形状自体の実装 (hit関数以外)
#include "geom3D.hpp"

namespace boom {
	namespace geo3d {
		AMat43 IModel::s_identityMat;

		float Sphere::bs_getArea() const {
			return 0;
		}
	}
}
