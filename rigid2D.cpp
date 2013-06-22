#include "geom2D.hpp"

namespace boom {
	namespace geo2d {
		// -------------------------- RPose --------------------------
		RPose::Value::Value(RPose& rp):
			Pose2D::Value(rp), vel(rp._vel), acc(rp._acc), rotVel(rp._rotVel), rotAcc(rp._rotAcc)
		{}
		RPose::Value::~Value() {
			reinterpret_cast<RPose*>(&_pose)->_setAsChanged();
		}

		RPose::RPose(): _finalInv() {
			_velAccum = _invAccum = getAccum()-1;
		}
		RPose::RPose(const RPose& rp): _finalInv() {
			// trivialなctorなのでベタコピー
			std::memcpy(this, &rp, sizeof(rp));
		}
		void RPose::identity() {
			Pose2D::identity();
			_vel = _acc = Vec2(0);
			_rotVel = _rotAcc = 0;
			_setAsChanged();
		}
		void RPose::_setAsChanged() {
			_velAccum = _invAccum = getAccum()-1;
		}

		void RPose::setVelocity(const Vec2& v) {
			_vel = v;
			++_velAccum;
		}
		void RPose::setRotVel(float m) {
			_rotVel = m;
			++_velAccum;
		}
		void RPose::setRotAccel(float a) {
			_rotAcc = a;
			++_velAccum;
		}

		const Vec2& RPose::getVelocity() const { return _vel; }
		const Vec2& RPose::getAccel() const { return _acc; }
		float RPose::getRotVel() const { return _rotVel; }
		float RPose::getRotAccel() const { return _rotAcc; }

		Vec2 RPose::getVelocAt(const Vec2& at) const {
			Vec2 ret(at - getOffset());
			return _vel + ret * cs_mRot90[0];
		}
		Vec2x2 RPose::getVelocities(const Vec2& at) const {
			return {_vel, getVelocAt(at)};
		}
		Vec2 RPose::toLocal(const Vec2& wpos) const {
			const auto& m = getFinalInv();
			auto v = wpos.asVec3(1) * m;
			return Vec2(v.x, v.y);
		}
		Vec2 RPose::toLocalDir(const Vec2& wdir) const {
			return wdir * spn::AMat22::Rotation(-getAngle());
		}
		Vec2 RPose::toWorld(const Vec2& lpos) const {
			const auto& m = getFinal();
			return lpos.asVec3(1) * m;
		}
		Vec2 RPose::toWorldDir(const Vec2& ldir) const {
			return ldir * spn::AMat22::Rotation(getAngle());
		}
		uint32_t RPose::getVelocityAccum() const {
			return _velAccum;
		}

		const AMat33& RPose::getFinalInv() const {
			auto ac = getAccum();
			if(_invAccum != ac) {
				_invAccum = ac;
				getFinal().convertA33().inversion(_finalInv);
			}
			return _finalInv;
		}
		RPose RPose::lerp(const RPose& p1, float t) const {
			RPose retR;
			new(&retR) Pose2D(((Pose2D&)*this).lerp(p1, t));

			retR._vel = _vel.l_intp(p1._vel, t);
			retR._acc = _acc.l_intp(p1._acc, t);
			retR._rotVel = spn::Lerp(_rotVel, p1._rotVel, t);
			retR._rotAcc = spn::Lerp(_rotAcc, p1._rotAcc, t);
			return retR;
		}
		RPose::Value RPose::refValue() {
			return Value(*this);
		}

		// -------------------------- Rigid --------------------------
		Rigid::Rigid(IModel::csptr sp): _spModel(sp) {}
		Rigid::Rigid(IModel::csptr sp, const RPose& pose): _spModel(sp), _pose(pose) {}

		Vec2 Rigid::support(const Vec2& dir) const {
			// dirを座標変換
			Vec2 tdir = _pose.toLocalDir(dir);
			return _spModel->support(tdir);
		}
		Vec2 Rigid::center() const {
			return _spModel->center();
		}
		uint32_t Rigid::getCID() const {
			return _spModel->getCID();
		}
		bool Rigid::isInner(const Vec2& pos) const {
			return false;
		}
		RPose& Rigid::refPose() {
			return _pose;
		}
		const RPose& Rigid::getPose() const {
			return _pose;
		}
		const IModel& Rigid::getModel() const {
			return *_spModel.get();
		}

		void Rigid::addR(const SPResist& sp) {
			for(int i=0 ; i<4 ; i++) {
				if(!_resist[i]) {
					_resist[i] = sp;
					return;
				}
			}
			assert(false);
		}
		RForce::F Rigid::resist(ColResult::CItrP itr) const {
			RForce::F res = {};
			for(auto& sp : _resist) {
				if(!sp)
					break;
				sp->resist(res, itr, *this);
			}
			return res;
		}

		// -------------------------- IResist --------------------------
		namespace resist {
			Gravity::Gravity(const Vec2& v): _grav(v) {}
			void Gravity::resist(RForce::F& acc, ColResult::CItrP itr, const Rigid& r) const {
				acc.linear += _grav;
			}

			Impact::Impact(const RCoeff& rc): _coeff(rc) {}
			void Impact::resist(RForce::F& acc, ColResult::CItrP itr, const Rigid& r) const {
				while(itr.first != itr.second) {
					//TODO: 衝突平面の計算
					auto f = CalcForce(r, *(itr.first->rigid),
									   itr.second->inner, _coeff, StLineCore(Vec2(0,0), Vec2(0,1)));
					acc += f.sdump;
					acc += f.fricD;
					++itr.first;
				}
			}
		}

		namespace itg {
			// -------------------------- Eular --------------------------
			int Eular::numOfIteration() const { return 1; }
			void Eular::advance(int pass, const RList& rlist, const ColResult& cr, float dt) {
				int nR = rlist.size();
				for(int i=0 ; i<nR ; i++) {
					auto* ptr = rlist[i].get();
					auto st = ptr->refPose().refValue();
					// 現フレームの加速度
					auto acc = ptr->resist(cr.getItem(i));
					st.acc += acc.linear;
					st.rotAcc += acc.torque;
					// 次のフレームの位置
					st.ofs += st.vel * dt;
					st.ang += st.rotVel * dt;
					// 次のフレームの速度
					st.vel += st.acc * dt;
					st.rotVel += st.rotAcc * dt;
				}
			}
			// -------------------------- ImpEular --------------------------
			int ImpEular::numOfIteration() const { return 2; }
			void ImpEular::beginIteration(int n) {
				// 剛体数分のメモリ*2を確保
				_tvalue.resize(n * 2);
			}
			void ImpEular::endIteration() {
				// メモリ解放
				decltype(_tvalue) tmp;
				std::swap(_tvalue, tmp);
			}
			void ImpEular::advance(int pass, const RList& rlist, const ColResult& cr, float dt) {
				auto *tv0 = &_tvalue[0];
				float dth2 = dt/2;
				int nR = rlist.size();
				if(pass == 0) {
					// value = 1つ前の(計算上の)状態
					for(int i=0 ; i<nR ; i++) {
						auto* ptr = rlist[i].get();
						auto dat = ptr->refPose().refValue();
						auto& ps = tv0[i];

						// 内部のメモリに書き込むと同時に出力
						ps = dat;
						auto acc = ptr->resist(cr.getItem(i));
						dat.ofs += dat.vel * dt;
						dat.vel += dat.acc * dt;
						dat.ang += dat.rotVel * dt;
						dat.rotVel += dat.rotAcc * dt;
						// (衝突判定結果は1フレーム遅れて出る為)
						ps.acc = acc.linear;
						ps.rotAcc = acc.torque;
					}
				} else {
					for(int i=0 ; i<nR ; i++) {
						auto* ptr = rlist[i].get();
						auto dat = ptr->refPose().refValue();
						auto& ps0 = tv0[i];

						dat.ofs = ps0.ofs + (ps0.vel + dat.vel) * dth2;
						dat.vel = ps0.vel + (ps0.acc + dat.acc) * dth2;
						dat.ang = ps0.ang + (ps0.rotVel + dat.rotVel) * dth2;
						dat.rotVel = ps0.rotVel + (ps0.rotAcc + dat.rotAcc) * dth2;
						auto acc = ptr->resist(cr.getItem(i));
						dat.acc = acc.linear;
						dat.rotAcc = acc.torque;
					}
				}
			}
			// -------------------------- RK4 --------------------------
			int RK4::numOfIteration() const { return 4; }
			void RK4::beginIteration(int n) {
				_tvalue.resize(n * 4);
			}
			void RK4::endIteration() {
				decltype(_tvalue) tmp;
				std::swap(_tvalue, tmp);
			}
			// TODO: 4回分も当たり判定していたら遅そうなので2回に留める案
			void RK4::advance(int pass, const RList& rlist, const ColResult& cr, float dt) {
				float dt2 = dt/2,
						dt6 = dt/6;
				int nR = rlist.size();
				auto *tv0 = &_tvalue[0],
					*tv1 = &_tvalue[nR],
					*tv2 = &_tvalue[2*nR],
					*tv3 = &_tvalue[3*nR];

				switch(pass) {
					case 0:
						for(int i=0 ; i<nR ; i++) {
							auto* ptr = rlist[i].get();
							auto dat = ptr->refPose().refValue();
							auto& ps = tv0[i];

							ps = dat;
							auto acc = ptr->resist(cr.getItem(i));		// 処理前の加速度
							ps.acc = acc.linear;
							ps.rotAcc = acc.torque;
							dat.ofs += dat.vel * dt2;
							dat.vel += ps.acc * dt2;
							dat.ang += dat.rotVel * dt2;
							dat.rotVel += ps.rotAcc * dt2;
						}
						break;
					case 1:
						for(int i=0 ; i<nR ; i++) {
							auto* ptr = rlist[i].get();
							auto dat = ptr->refPose().refValue();
							auto &ps0 = tv0[i],
								&ps1 = tv1[i];

							ps1 = dat;
							auto acc = ptr->resist(cr.getItem(i));		// ps1の加速度
							ps1.acc = acc.linear;
							ps1.rotAcc = acc.torque;
							dat.ofs = ps0.ofs + ps1.vel * dt2;
							dat.vel = ps0.vel + ps1.acc * dt2;
							dat.ang = ps0.ang + ps1.rotVel * dt2;
							dat.rotVel = ps0.rotVel + ps1.rotAcc * dt2;
						}
						break;
					case 2:
						for(int i=0 ; i<nR ; i++) {
							auto* ptr = rlist[i].get();
							auto dat = ptr->refPose().refValue();
							auto &ps0 = tv0[i],
								&ps2 = tv2[i];

							ps2 = dat;
							auto acc = ptr->resist(cr.getItem(i));		// ps2の加速度
							ps2.acc = acc.linear;
							ps2.rotAcc = acc.torque;
							dat.ofs = ps0.ofs + ps2.vel * dt;
							dat.vel = ps0.vel + ps2.acc * dt;
							dat.ang = ps0.ang + ps2.rotVel * dt;
							dat.rotVel = ps0.rotVel + ps2.rotAcc * dt;
						}
						break;
					case 3:
						for(int i=0 ; i<nR ; i++) {
							auto* ptr = rlist[i].get();
							auto dat = ptr->refPose().refValue();
							auto &ps0 = tv0[i],
								&ps1 = tv1[i],
								&ps2 = tv2[i],
								&ps3 = tv3[i];

							ps3 = dat;
							auto acc = ptr->resist(cr.getItem(i));
							ps3.acc = acc.linear;
							ps3.rotAcc = acc.torque;

							dat.ofs = ps0.ofs + (ps0.vel + ps1.vel*2 + ps2.vel*2 + ps3.vel) * dt6;
							dat.vel = ps0.vel + (ps0.acc + ps1.acc*2 + ps2.acc*2 + ps3.acc) * dt6;
							dat.ang = ps0.ang + (ps0.rotVel + ps1.rotVel*2 + ps2.rotVel*2 + ps3.rotVel) * dt6;
							dat.rotVel = ps0.rotVel + (ps0.rotAcc + ps1.rotAcc*2 + ps2.rotAcc*2 + ps3.rotAcc) * dt6;
							dat.acc = ps0.acc;
							dat.rotAcc = ps0.rotAcc;
						}
						break;
				}
			}
		}
		// -------------------------- ColResult --------------------------
		ColResult::ColResult() {
			clear();
		}
		ColResult::ColResult(ColResult&& cr):
			_array(std::forward<ItemArray>(cr._array)),
			_cursor(std::forward<CursorMap>(cr._cursor)),
			_curID(cr._curID),
			_from(cr._from)
		{}
		ColResult& ColResult::operator = (ColResult&& cr) {
			std::swap(_array, cr._array);
			std::swap(_cursor, cr._cursor);
			_curID = cr._curID;
			_from = cr._from;
			return *this;
		}
		void ColResult::setCurrent(int id) {
			if(id != _curID) {
				int nA = _array.size();
				if(_from < nA)
					_cursor.emplace(_curID, std::make_pair(uint16_t(_from), uint16_t(nA)));
				_curID = id;
				_from = nA;
			}
		}
		void ColResult::pushItem(const Rigid* r, const Vec2& p) {
			assert(_curID >= 0);
			_array.push_back(Item{r, p});
		}
		ColResult::CItrP ColResult::getItem(int id) const {
			auto itr = _cursor.find(id);
			if(itr == _cursor.end())
				return CItrP();
			auto pitr = _array.begin()+itr->second.first;
			return CItrP(pitr, pitr+itr->second.second);
		}
		void ColResult::clear() {
			_array.clear();
			_cursor.clear();
			_curID = -1;
			_from = 0;
		}

		// -------------------------- RigidMgr --------------------------
		RigidMgr::RigidMgr(IItg::csptr itg): _itg(itg) {}
		void RigidMgr::add(Rigid::csptr sp) {
			_rlist.push_back(sp);
		}
		void RigidMgr::_checkCollision() {
			_cresult.clear();
			// TODO: 4分木を使ってマシな動作速度にする
			int nR = _rlist.size();
			for(int i=0 ; i<nR-1 ; i++) {
				const auto* pr0 = _rlist[i].get();
				const auto& mdl0 = pr0->getModel();
				auto c0 = mdl0.getBCircle()
							* pr0->getPose().getFinal();
				_cresult.setCurrent(i);

				for(int j=i+1 ; j<nR ; j++) {
					const auto* pr1 = _rlist[j].get();
					const auto& mdl1 = pr1->getModel();
					auto c1 = mdl1.getBCircle()
								* pr1->getPose().getFinal();

					// 境界球による判定
					if(c0.hit(c1)) {
						// GJKによる判定
						GSimplex gs(mdl0, mdl1);
						if(gs.getResult()) {
							// ついでに内部点も格納
							_cresult.pushItem(pr1, gs.getInner());
						}
					}
				}
			}
			_cresult.setCurrent(-1);
		}
		void RigidMgr::simulate(float dt) {
			IItg* itg = _itg.get();
			int nR = _rlist.size(),
				nItr = itg->numOfIteration();

			itg->beginIteration(nR);
			for(int i=0 ; i<nItr ; i++) {
				// 当たり判定結果は一回分だけとっておけば良い
				_checkCollision();
				itg->advance(i, _rlist, _cresult, dt);
			}
			itg->endIteration();
		}
		int RigidMgr::getNRigid() const {
			return _rlist.size();
		}
		Rigid::csptr RigidMgr::getRigid(int n) const {
			return _rlist[n];
		}
	}
}