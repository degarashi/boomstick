//! テストコードの記述
#include "geom3D.hpp"

int main() {
	boom::geo3d::Narrow na;
	boom::geo3d::Narrow::Initialize();
	boom::geo3d::SphereM sp;
	boom::geo3d::CapsuleM cp;
	boom::geo3d::Narrow::Hit(&sp, &cp);
	return 0;
}
