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
			_identitySingle();
		}
		RPose::RPose(const RPose& rp): _finalInv() {
			// trivialなctorなのでベタコピー
			std::memcpy(this, &rp, sizeof(rp));
		}

		void RPose::identity() {
			Pose2D::identity();
			_identitySingle();
		}
		void RPose::_identitySingle() {
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
			return _vel + ret * cs_mRot90[0] * _rotVel;
		}
		Vec2x2 RPose::getVelocities(const Vec2& at) const {
			return {_vel, getVelocAt(at)};
		}
		Vec2 RPose::toLocal(const Vec2& wpos) const {
			const auto& m = getToLocal();
			auto v = wpos.asVec3(1) * m;
			return Vec2(v.x, v.y);
		}
		Vec2 RPose::toLocalDir(const Vec2& wdir) const {
			return wdir * spn::AMat22::Rotation(-getAngle());
		}
		Vec2 RPose::toWorld(const Vec2& lpos) const {
			const auto& m = getToWorld();
			return lpos.asVec3(1) * m;
		}
		Vec2 RPose::toWorldDir(const Vec2& ldir) const {
			return ldir * spn::AMat22::Rotation(getAngle());
		}
		uint32_t RPose::getVelocityAccum() const {
			return _velAccum;
		}

		const AMat32& RPose::getToLocal() const {
			auto ac = getAccum();
			if(_invAccum != ac) {
				_invAccum = ac;
				AMat33 tmp(getToWorld().convertA33());
				tmp.invert();
				_finalInv = tmp.convertA32();
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

		// -------------------------- TR_Mat --------------------------
		TR_Mat::TR_Mat(const RPose& rp): _mToLocal(rp.getToLocal()), _mToWorld(rp.getToWorld()) {}
		TR_Mat::TR_Mat(const TR_Mat& t, tagInverse): _mToLocal(t._mToWorld), _mToWorld(t._mToLocal) {}
		TR_Mat::TR_Mat(const RPose& rp, tagInverse): _mToLocal(rp.getToWorld()), _mToWorld(rp.getToLocal()) {}

		Vec2 TR_Mat::toLocal(const Vec2& v) const { return v.asVec3(1) * _mToLocal; }
		Vec2 TR_Mat::toLocalDir(const Vec2& v) const { return v.asVec3(0) * _mToLocal; }
		Vec2 TR_Mat::toWorld(const Vec2& v) const { return v.asVec3(1) * _mToWorld; }
		Vec2 TR_Mat::toWorldDir(const Vec2& v) const { return v.asVec3(0) * _mToWorld; }
		const AMat32& TR_Mat::getToLocal() const { return _mToLocal; }
		const AMat32& TR_Mat::getToWorld() const { return _mToWorld; }

		// -------------------------- TModel --------------------------
		template <class MDL, class BASE>
		TModel<MDL,BASE>::TModel(const MDL& mdl, const BASE& base): BASE(base), _model(mdl) {}
		template <class MDL, class BASE>
		TModel<MDL,BASE>::TModel(const MDL &mdl): _model(mdl) {}

		#define DEF_TMODEL(typ)	template <class MDL,class BASE> typ TModel<MDL,BASE>
		DEF_TMODEL(const IModel&)::_getModel(const IModel& mdl) const { return mdl; }
		DEF_TMODEL(const IModel&)::_getModel(IModel::csptr sp) const { return *sp.get(); }
		DEF_TMODEL(Vec2)::support(const Vec2& dir) const {
			// dirをローカルに座標変換してサポート写像した後、またワールド座標系に戻す
			Vec2 pos = _getModel(_model).support(BASE::toLocalDir(dir));
			return BASE::toWorld(pos);
		}
		DEF_TMODEL(Vec2)::center() const {
			Vec2 res = _getModel(_model).center();
			return BASE::toWorld(res);
		}
		DEF_TMODEL(bool)::isInner(const Vec2& pos) const {
			Vec2 tpos = BASE::toLocal(pos);
			return _getModel(_model).isInner(tpos);
		}
		DEF_TMODEL(const MDL&)::getModel() const {
			return _model;
		}
		DEF_TMODEL(uint32_t)::getCID() const {
			return _getModel(_model).getCID();
		}
		DEF_TMODEL(CircleCore)::getBCircle() const {
			CircleCore c = _getModel(_model).getBCircle();
			return c * BASE::getToWorld();
		}
		DEF_TMODEL(PointL)::getOverlappingPoints(const IModel& mdl, const Vec2& inner) const {
			// innerとmdlをこちらのローカル座標系に変換
			TModelR<TR_Mat> t_mdl(mdl, (const BASE&)*this);
			Vec2 t_inner(BASE::toLocal(inner));
			return _getModel(_model).getOverlappingPoints(t_mdl, t_inner);
		}
		DEF_TMODEL(int)::getNPoints() const { return _getModel(_model).getNPoints(); }
		DEF_TMODEL(Vec2)::getPoint(int n) const {
			Vec2 p = _getModel(_model).getPoint(n);
			return BASE::toWorld(p);
		}
		DEF_TMODEL(IModel::CPos)::checkPosition(const Vec2& pos) const {
			// posをローカル座標系に変換
			Vec2 tpos = BASE::toLocal(pos);
			return _getModel(_model).checkPosition(tpos);
		}
		DEF_TMODEL(Convex2)::splitTwo(const StLineCore& line) const {
			StLineCore tline(BASE::toLocal(line.pos),
							BASE::toLocalDir(line.dir));
			auto res = _getModel(_model).splitTwo(tline);
			res.first *= BASE::getToWorld();
			res.second *= BASE::getToWorld();
			return std::move(res);
		}
		#undef DEF_TMODEL
		template class TModel<IModel::sptr, TR_Mat>;
		template class TModel<const IModel&, TR_Mat>;
		template class TModel<IModel::sptr, RPose>;
		template class TModel<const IModel&, RPose>;

		// -------------------------- Rigid --------------------------
		RPose& Rigid::refPose() { return *this; }
		const RPose& Rigid::getPose() const { return *this; }

		void Rigid::addR(const SPResist& sp) {
			for(int i=0 ; i<4 ; i++) {
				if(!_resist[i]) {
					_resist[i] = sp;
					return;
				}
			}
			assert(false);
		}
		RForce::F Rigid::resist(int index, const RigidMgr::CRes& cr) const {
			RForce::F res = {};
			for(auto& sp : _resist) {
				if(!sp)
					break;
				sp->resist(res, *this, index, cr);
			}
			return res;
		}

		// -------------------------- IResist --------------------------
		namespace resist {
			Gravity::Gravity(const Vec2& v): _grav(v) {}
			void Gravity::resist(RForce::F& acc, const Rigid& /*r*/, int /*index*/, const RigidMgr::CRes& /*cr*/) const {
				acc.linear += _grav;
			}

			void Impact::resist(RForce::F& acc, const Rigid& r, int index, const RigidMgr::CRes& cr) const {
				auto fc = cr.getInfo().getInfo(index, 0);
				auto& ff = fc->sdump;
				auto& ff2 = fc->fricD;
				acc += ff;
				acc += ff2;
			}
		}
		namespace itg {
			// -------------------------- Eular --------------------------
			int Eular::numOfIteration() const { return 1; }
			void Eular::advance(int pass, const RList& rlist, const RigidMgr::CRes& cr, float dt) {
				int nR = rlist.size();
				for(int i=0 ; i<nR ; i++) {
					auto* ptr = rlist[i].get();
					auto st = ptr->refPose().refValue();
					// 現フレームの加速度
					auto acc = ptr->resist(i, cr);
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
			void ImpEular::advance(int pass, const RList& rlist, const RigidMgr::CRes& cr, float dt) {
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
						auto acc = ptr->resist(i, cr);
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
						auto acc = ptr->resist(i, cr);
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
			void RK4::advance(int pass, const RList& rlist, const RigidMgr::CRes& cr, float dt) {
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
							auto acc = ptr->resist(i, cr);		// 処理前の加速度
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
							auto acc = ptr->resist(i, cr);		// ps1の加速度
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
							auto acc = ptr->resist(i, cr);		// ps2の加速度
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
							auto acc = ptr->resist(i, cr);
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

		// -------------------------- RigidMgr --------------------------
		RigidMgr::RigidMgr(IItg::csptr itg): _itg(itg) {}
		void RigidMgr::add(Rigid::csptr sp) {
			_rlist.push_back(sp);
		}
		void RigidMgr::_checkCollision() {
			_cresult.clear();
			// TODO: 4分木を使ってマシな動作速度にする
			int nR = _rlist.size();
			_cresult.setNumObjects(nR);
			for(int i=0 ; i<nR-1 ; i++) {
				const auto* pr0 = _rlist[i].get();
				auto c0 = pr0->getBCircle();

				for(int j=i+1 ; j<nR ; j++) {
					const auto* pr1 = _rlist[j].get();
					auto c1 = pr1->getBCircle();

					// 境界球による判定
					if(c0.hit(c1)) {
						// GJKによる判定
						GSimplex gs(*pr0, *pr1);
						if(gs.getResult()) {
							// ついでに内部点も格納
							//TODO: 衝突平面の計算
							auto f = CalcForce(*pr0, *pr1,
											   gs.getInner(), _coeff, StLineCore(Vec2(0,0), Vec2(0,1)));
							_cresult.setFlagWithInfo(i,j, f);
						}
					}
				}
			}
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
		void RigidMgr::setCoeff(const RCoeff& c) {
			_coeff = c;
		}
		int RigidMgr::getNRigid() const {
			return _rlist.size();
		}
		Rigid::csptr RigidMgr::getRigid(int n) const {
			return _rlist[n];
		}
	}
}