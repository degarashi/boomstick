//! テストコードの記述
#include "geom3D.hpp"
#include "geom2D.hpp"

int main() {
	boom::geo3d::Narrow na;
	boom::geo3d::Narrow::Initialize();
	boom::geo3d::SphereM sp;
	boom::geo3d::CapsuleM cp;
	boom::geo3d::Narrow::Hit(&sp, &cp);
	boom::geo2d::Narrow na2;
	boom::geo2d::Narrow::Initialize();

	boom::geo2d::CircleM m;
	return 0;
}
