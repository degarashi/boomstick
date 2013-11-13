#include "geom2D.hpp"
#include "rigid2D.hpp"

namespace boom {
	namespace geo2d {
		RForce::F& RForce::F::operator += (const F& f) {
			linear += f.linear;
			torque += f.torque;
			return *this;
		}
		RForce::F& RForce::F::operator *= (float s) {
			linear *= s;
			torque *= s;
			return *this;
		}
		RForce::F RForce::F::operator * (float s) const {
			F ret(*this);
			return ret *= s;
		}

		RForce& RForce::operator +=(const RForce& rf) {
			sdump += rf.sdump;
			fricD += rf.fricD;
			return *this;
		}
		RForce& RForce::operator *= (float s) {
			sdump *= s;
			fricD *= s;
			return *this;
		}
		RForce RForce::operator * (float s) const {
			RForce ret(*this);
			return ret *= s;
		}

		std::ostream& operator << (std::ostream& os, const RForce::F& f) {
			return os << "linear: " << f.linear << std::endl
			<< "torque: " << f.torque;
		}
		std::ostream& operator << (std::ostream& os, const RForce& f) {
			return os << "sdump: " << f.sdump << std::endl
			<< "fricD: " << f.fricD;
		}

		// -------------------------- RPose --------------------------
		RPose::Value::Value(RPose& rp):
			Pose2D::Value(rp), vel(rp._vel), acc(rp._acc), rotVel(rp._rotVel), rotAcc(rp._rotAcc)
		{}
		RPose::Value::~Value() {
			reinterpret_cast<RPose*>(&_pose)->_setAsChanged();
		}
		RPose::Value::Value(const Value& v): spn::Pose2D::Value(v),
			vel(v.vel), acc(v.acc), rotVel(v.rotVel), rotAcc(v.rotAcc) {}
		RPose::Value& RPose::Value::operator = (const TValue& tv) {
			static_cast<spn::Pose2D::Value&>(*this) = static_cast<const spn::Pose2D::TValue&>(tv);
			vel = tv.vel;
			acc = tv.acc;
			rotVel = tv.rotVel;
			rotAcc = tv.rotAcc;
			return * this;
		}
		RPose::TValue& RPose::TValue::operator = (const Value& v) {
			static_cast<spn::Pose2D::TValue&>(*this) = static_cast<const spn::Pose2D::Value&>(v);
			vel = v.vel;
			acc = v.acc;
			rotVel = v.rotVel;
			rotAcc = v.rotAcc;
			return *this;
		}
		RPose::TValue::TValue(const RPose& rp) {
			static_cast<spn::Pose2D::TValue&>(*this) = static_cast<const spn::Pose2D::TValue&>(rp);
			vel = rp.getVelocity();
			acc = rp.getAccel();
			rotVel = rp.getRotVel();
			rotAcc = rp.getRotAccel();
		}

		RPose::RPose(): _finalInv() {
			_identitySingle();
		}
		RPose::RPose(const RPose& rp): _finalInv() {
			// trivialなctorなのでベタコピー
			std::memcpy(this, &rp, sizeof(rp));
		}
		RPose::RPose(const spn::Pose2D& ps): _finalInv() {
			std::memcpy(static_cast<spn::Pose2D*>(this), &ps, sizeof(ps));
			_identitySingle();
		}
		RPose& RPose::operator =(const spn::Pose2D& ps) {
			static_cast<spn::Pose2D&>(*this) = ps;
			return *this;
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
		void RPose::setAccel(const Vec2& a) {
			_acc = a;
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
			retR._rotVel = spn::Lerp(*_rotVel, *p1._rotVel, t);
			retR._rotAcc = spn::Lerp(*_rotAcc, *p1._rotAcc, t);
			return retR;
		}
		RPose::Value RPose::refValue() {
			return Value(*this);
		}

		// -------------------------- TR_Mat --------------------------
		TR_Mat::tagInverse TR_Mat::TagInverse;
		TR_Mat::TR_Mat(const AMat32& m): _mToWorld(m) {
			auto mt = m.convert33();
			mt.invert();
			_mToLocal = mt.convertA32();
		}

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
		TModel<MDL,BASE>::TModel(const MDL &mdl): _model(mdl) {}

		#define DEF_TMODEL(typ)	template <class MDL,class BASE> typ TModel<MDL,BASE>
		DEF_TMODEL(const IModel&)::_getModel(const IModel& mdl) const { return mdl; }
		DEF_TMODEL(const IModel&)::_getModel(const HLMdl& mdl) const { return *mdl.cref().get(); }
		DEF_TMODEL(Vec2)::support(const Vec2& dir) const {
			// dirをローカルに座標変換してサポート写像した後、またワールド座標系に戻す
			Vec2 pos = _getModel(_model).support(BASE::toLocalDir(dir));
			return BASE::toWorld(pos);
		}
		DEF_TMODEL(Vec2)::getCenter() const {
			Vec2 res = _getModel(_model).getCenter();
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
		DEF_TMODEL(float)::getArea(bool bInv) const {
			return _getModel(_model).getArea(bInv);
		}
		DEF_TMODEL(float)::getInertia(bool bInv) const {
			// 返すのは常に重心周りの慣性モーメントなので座標変換は不要
			return _getModel(_model).getInertia(bInv);
		}
		DEF_TMODEL(IModel::PosL)::getOverlappingPoints(const IModel& mdl, const Vec2& inner) const {
			// innerとmdlをこちらのローカル座標系に変換
			TModelR<TR_Mat> t_mdl(mdl, (const BASE&)*this, TR_Mat::TagInverse);
			Vec2 t_inner(BASE::toLocal(inner));
			// 結果をグローバル座標系へ変換
			auto ret = _getModel(_model).getOverlappingPoints(t_mdl, t_inner);
			const auto& mat = BASE::getToWorld();
			for(auto& p : ret.second)
				p = p.asVec3(1) * mat;
			return std::move(ret);
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
			const auto& mat = BASE::getToWorld();
			res.first *= mat;
			res.second *= mat;
			return std::move(res);
		}
		DEF_TMODEL(std::ostream&)::dbgPrint(std::ostream& os) const {
			int nP = getNPoints();
			if(nP > 0) {
				for(int i=0 ; i<nP-1 ; i++)
					os << getPoint(i) << std::endl;
				os << getPoint(nP-1);
			}
			return os;
		}
		#undef DEF_TMODEL
		template class TModel<HLMdl, TR_Mat>;
		template class TModel<const IModel&, TR_Mat>;
		template class TModel<HLMdl, RPose>;
		template class TModel<const IModel&, RPose>;

		// -------------------------- Rigid --------------------------
		RPose& Rigid::refPose() { return *this; }
		const RPose& Rigid::getPose() const { return *this; }

		void Rigid::addR(const SPResist& sp, uint32_t id) {
			for(int i=0 ; i<NUM_RESIST ; i++) {
				if(!_resist[i]) {
					_resist[i] = sp;
					_rID[i] = id;
					return;
				}
			}
			assert(false);
		}
		namespace {
			Rigid::SPResist g_invalidSP;
		}
		const Rigid::SPResist& Rigid::getR(uint32_t id) {
			for(int i=0 ; i<NUM_RESIST ; i++) {
				if(_resist[i] && _rID[i] == id)
					return _resist[i];
			}
			return g_invalidSP;
		}
		void Rigid::remR(uint32_t id) {
			for(int i=0 ; i<NUM_RESIST ; i++) {
				if(!_resist[i])
					break;
				if(_rID[i] == id) {
					_resist[i] = nullptr;
					for(int j=i ; j<NUM_RESIST-1 ; j++) {
						_resist[j] = std::move(_resist[j+1]);
						_rID[j] = _rID[j+1];
					}
					--i;
				}
			}
		}
		void Rigid::remRAll() {
			for(auto& sp : _resist)
				sp = nullptr;
		}

		RForce::F Rigid::resist(int index, const CResult& cr) const {
			RForce::F res = {};
			for(auto& sp : _resist) {
				if(!sp)
					break;
				sp->resist(res, *this, index, cr);
			}
			return res;
		}
		void Rigid::applyForce(const RForce::F& f) {
			RPose& rp = refPose();
			rp.setAccel(rp.getAccel() + f.linear * Rcp22Bit(getArea()));
			rp.setRotAccel(rp.getRotAccel() + f.torque * Rcp22Bit(getInertia()));
		}
		void Rigid::addForce(const Vec2& wpos, const Vec2& f) {
			applyForce(CalcForce(getPose(), wpos, f));
		}
		void Rigid::addLinearForce(const Vec2& f) {
			RPose& rp = refPose();
			rp.setAccel(rp.getAccel() + f * Rcp22Bit(getArea()));
		}
		RForce::F Rigid::CalcForce(const RPose& rp, const Vec2& wpos, const Vec2& f) {
			RForce::F ret;
			// linear
			ret.linear = f;

			Vec2 tpos = rp.getOffset() - wpos;
			float l = tpos.len_sq();
			if(l > 1e-8f) {
				// torque
				ret.torque = f.cw(tpos);
			} else
				ret.torque = 0;
			return ret;
		}

		// -------------------------- IResist --------------------------
		namespace resist {
			Gravity::Gravity(const Vec2& v): _grav(v) {}
			void Gravity::setGravity(const Vec2& v) {
				_grav = v;
			}
			void Gravity::resist(RForce::F& acc, const Rigid& r, int /*index*/, const CResult& /*cr*/) const {
				float inv_area = Rcp22Bit(r.getModel().cref()->getArea(false));
				acc.linear += _grav * Rcp22Bit(inv_area);
			}

			void Impact::resist(RForce::F& acc, const Rigid& r, int index, const CResult& cr) const {
				auto& mdl = *r.getModel().cref().get();
				float inv_area = Rcp22Bit(mdl.getArea(false));
				float inv_inertia = Rcp22Bit(mdl.getInertia(false));
				auto fc = cr.getInfo().getInfo(index, 0);
				auto& ff = fc->sdump;
				auto& ff2 = fc->fricD;
				acc.linear += (ff.linear + ff2.linear) * inv_area;
				acc.torque += (ff.torque + ff2.torque) * inv_inertia;
			}

			Air::Air(float cLinear, float cRot): _cLinear(cLinear), _cRot(cRot) {}
			void Air::setAir(float cLinear, float cRot) {
				_cLinear = cLinear;
				_cRot = cRot;
			}
			void Air::resist(RForce::F& acc, const Rigid& r, int index, const CResult& cr) const {
				const RPose& ps = r.getPose();
				acc.linear -= ps.getVelocity() * _cLinear;
				acc.torque -= ps.getRotVel() * _cRot;
			}
		}
		namespace itg {
			// -------------------------- Eular --------------------------
			int Eular::numOfIteration() const { return 1; }
			void Eular::advance(int pass, int offset, RItr itr, RItr itrE, const CResult& cr, float dt) {
				for(int i=0 ; itr!=itrE ; ++i,++itr) {
					Rigid& r = static_cast<Rigid&>(*itr);
					auto st = r.refPose().refValue();
					// 現フレームの加速度
					auto acc = r.resist(i, cr);
					st.acc += acc.linear * dt / r.getArea(false);
					st.rotAcc += acc.torque * dt / r.getInertia(false);
					// 次のフレームの速度
					st.vel += st.acc * dt;
					st.rotVel += st.rotAcc * dt;
					// 次のフレームの位置
					st.ofs += st.vel * dt;
					st.ang += st.rotVel * dt;
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
			void ImpEular::advance(int pass, int offset, RItr itr, RItr itrE, const CResult& cr, float dt) {
				auto *tv0 = &_tvalue[0];
				float dth2 = dt/2;
				int cur = offset;
				if(pass == 0) {
					// value = 1つ前の(計算上の)状態
					while(itr != itrE) {
						Rigid& r = static_cast<Rigid&>(*itr);
						auto dat = r.refPose().refValue();
						auto& ps = tv0[cur];

						// 内部のメモリに書き込むと同時に出力
						ps = dat;
						auto acc = r.resist(cur, cr);
						dat.ofs += dat.vel * dt;
						dat.vel += dat.acc * dt;
						dat.ang += dat.rotVel * dt;
						dat.rotVel += dat.rotAcc * dt;
						// (衝突判定結果は1フレーム遅れて出る為)
						ps.acc = acc.linear / r.getArea(false);
						ps.rotAcc = acc.torque / r.getInertia(false);

						++itr;
						++cur;
					}
				} else {
					while(itr != itrE) {
						Rigid& r = static_cast<Rigid&>(*itr);
						auto dat = r.refPose().refValue();
						auto& ps0 = tv0[cur];

						dat.ofs = ps0.ofs + (ps0.vel + dat.vel) * dth2;
						dat.vel = ps0.vel + (ps0.acc + dat.acc) * dth2;
						dat.ang = ps0.ang + (ps0.rotVel + dat.rotVel) * dth2;
						dat.rotVel = ps0.rotVel + (ps0.rotAcc + dat.rotAcc) * dth2;
						auto acc = r.resist(cur, cr);
						dat.acc = acc.linear / r.getArea(false);
						dat.rotAcc = acc.torque / r.getInertia(false);

						++itr;
						++cur;
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
			void RK4::advance(int pass, int offset, RItr itr, RItr itrE, const CResult& cr, float dt) {
				float dt2 = dt/2,
						dt6 = dt/6;
				int cur = offset,
					nR = _tvalue.size()/4;
				auto *tv0 = &_tvalue[0],
					*tv1 = &_tvalue[nR],
					*tv2 = &_tvalue[2*nR],
					*tv3 = &_tvalue[3*nR];

				switch(pass) {
					case 0:
						while(itr != itrE) {
							Rigid& r = static_cast<Rigid&>(*itr);
							auto dat = r.refPose().refValue();
							auto& ps = tv0[cur];

							ps = dat;
							dat.ofs += dat.vel * dt2;
							dat.vel += dat.acc * dt2;
							dat.ang += dat.rotVel * dt2;
							dat.rotVel += dat.rotAcc * dt2;
							auto acc = r.resist(cur, cr);		// 処理前の加速度
							ps.acc = acc.linear / r.getArea(false);
							ps.rotAcc = acc.torque / r.getInertia(false);

							++itr;
							++cur;
						}
						break;
					case 1:
						while(itr != itrE) {
							Rigid& r = static_cast<Rigid&>(*itr);
							auto dat = r.refPose().refValue();
							auto &ps0 = tv0[cur],
								&ps1 = tv1[cur];

							ps1 = dat;
							auto acc = r.resist(cur, cr);		// ps1の加速度
							ps1.acc = acc.linear / r.getArea(false);
							ps1.rotAcc = acc.torque / r.getInertia(false);
							dat.ofs = ps0.ofs + ps1.vel * dt2;
							dat.vel = ps0.vel + ps1.acc * dt2;
							dat.ang = ps0.ang + ps1.rotVel * dt2;
							dat.rotVel = ps0.rotVel + ps1.rotAcc * dt2;

							++itr;
							++cur;
						}
						break;
					case 2:
						while(itr != itrE) {
							Rigid& r = static_cast<Rigid&>(*itr);
							auto dat = r.refPose().refValue();
							auto &ps0 = tv0[cur],
								&ps2 = tv2[cur];

							ps2 = dat;
							auto acc = r.resist(cur, cr);		// ps2の加速度
							ps2.acc = acc.linear / r.getArea(false);
							ps2.rotAcc = acc.torque / r.getInertia(false);
							dat.ofs = ps0.ofs + ps2.vel * dt;
							dat.vel = ps0.vel + ps2.acc * dt;
							dat.ang = ps0.ang + ps2.rotVel * dt;
							dat.rotVel = ps0.rotVel + ps2.rotAcc * dt;

							++itr;
							++cur;
						}
						break;
					case 3:
						while(itr != itrE) {
							Rigid& r = static_cast<Rigid&>(*itr);
							auto dat = r.refPose().refValue();
							auto &ps0 = tv0[cur],
								&ps1 = tv1[cur],
								&ps2 = tv2[cur],
								&ps3 = tv3[cur];

							ps3 = dat;
							auto acc = r.resist(cur, cr);
							ps3.acc = acc.linear / r.getArea(false);
							ps3.rotAcc = acc.torque / r.getInertia(false);

							dat.ofs = ps0.ofs + (ps0.vel + ps1.vel*2 + ps2.vel*2 + ps3.vel) * dt6;
							dat.vel = ps0.vel + (ps0.acc + ps1.acc*2 + ps2.acc*2 + ps3.acc) * dt6;
							dat.ang = ps0.ang + (ps0.rotVel + ps1.rotVel*2 + ps2.rotVel*2 + ps3.rotVel) * dt6;
							dat.rotVel = ps0.rotVel + (ps0.rotAcc + ps1.rotAcc*2 + ps2.rotAcc*2 + ps3.rotAcc) * dt6;
							dat.acc = ps0.acc;
							dat.rotAcc = ps0.rotAcc;

							++itr;
							++cur;
						}
						break;
				}
			}
		}

		// -------------------------- RigidMgr --------------------------
		RigidMgr::RigidMgr(IItg::csptr itg): _itg(itg) {}
		void RigidMgr::setCoeff(const RCoeff& c) {
			_coeff = c;
		}
		void RigidMgr::_checkCollision() {
			_cresult.clear();
			_cresult.setNumObjects(_broadC.getNumObj());
			_broadC.broadCollision(*this, NarrowC());
		}
		RigidMgr::id_type RigidMgr::addA(HRig hRig) {
			return _broadC.add(BroadC::TYPE_A, hRig);
		}
		RigidMgr::id_type RigidMgr::addB(HRig hRig) {
			return _broadC.add(BroadC::TYPE_B, hRig);
		}
		void RigidMgr::remA(id_type id) {
			_broadC.rem(BroadC::TYPE_A, id);
		}
		void RigidMgr::remB(id_type id) {
			_broadC.rem(BroadC::TYPE_B, id);
		}
		typename RigidMgr::BroadC::const_iterator RigidMgr::cbeginA() const {
			return _broadC.cbegin(BroadC::TYPE_A);
		}
		typename RigidMgr::BroadC::const_iterator RigidMgr::cendA() const {
			return _broadC.cend(BroadC::TYPE_A);
		}
		void RigidMgr::simulate(float dt) {
			IItg* itg = _itg.get();
			int nR = _broadC.getNumObj(),
				nItr = itg->numOfIteration();
			itg->beginIteration(nR);
			for(int i=0 ; i<nItr ; i++) {
				_broadC.iterate([](ERig& er) {
					// ヒットカウントを初期化
					er.resetHitCount();
				});
				// 当たり判定結果は一回分だけとっておけば良い
				_checkCollision();
				int offset = 0;
				for(int j=0 ; j<BroadC::NUM_TYPE ; j++) {
					auto typ = static_cast<BroadC::TYPE>(j);
					itg->advance(i, offset, _broadC.begin(typ), _broadC.end(typ), _cresult, dt);
					offset += _broadC.getNumObj(typ);
				}
				_broadC.iterate([this](ERig& er) {
					// 不要なエントリを削除
					er.removeOld(_npPool);
				});
			}
			itg->endIteration();
			_broadC.iterate([this](ERig& er) {
				// 次回の為に1フレーム前の姿勢を保存
				er.savePose();
			});
		}
		void RigidMgr::operator ()(int id0, int id1, ERig& er0_, ERig& er1_, const Vec2& inner) {
			using spn::Pose2D;

			HRig hR0 = er0_.hRig,
				hR1 = er1_.hRig;
			ERig *er0 = &er0_,
				*er1 = &er1_;
			// ハンドルのインデックス値が小さいほうがエントリを格納する
			if(hR0.getIndex() > hR1.getIndex()) {
				std::swap(hR0, hR1);
				std::swap(er0, er1);
				std::swap(id0, id1);
			}

			Rigid &r0 = *hR0.ref(),
					&r1 = *hR1.ref();
			auto b0 = er0->addHit(hR1, _npPool);
			boost::optional<GEpa> epa;
			if(b0) {
				// 初回の接触: 一手前の姿勢での最近傍対から接触法線を求める
				const IModel &mdl0 = *r0.getModel().cref(),
							&mdl1 = *r1.getModel().cref();
				TModelR<RPose>	tm0(mdl0, er0->prePose),
								tm1(mdl1, er1->prePose);

				epa = boost::in_place(tm0, tm1);
				// 前回の姿勢で衝突していたら配置当初から衝突していた可能性が高い
				if(epa->getResult() || epa->getNearestPair().first.distance(epa->getNearestPair().second) < 1e-6f) {
					// オブジェクトの中心座標を元に再度トライ
					CircleCore cc0 = r0.getBCircle(),
								cc1 = r1.getBCircle();

					Vec2 dir(cc0.vCenter - cc1.vCenter);
					float dist = dir.length();
					if(dist < 1e-5f)
						dir = Vec2(1,0);
					else
						dir *= Rcp22Bit(dist);

					constexpr float MARGIN = 1.05f;
					tm0.setOfs(cc1.vCenter + dir * (cc0.fRadius + cc1.fRadius) * MARGIN);
					epa = boost::in_place(tm0, tm1);
					assert(!epa->getResult());
				}
			} else {
				// 2回目以降の接触: 前回の法線の向きに物体を適当に離して最近傍対を求める
				Pose2D ps0 = static_cast<const Pose2D&>(r0.getPose()),
						ps1 = static_cast<const Pose2D&>(r1.getPose());
				float dist = (r0.getBCircle().fRadius + r1.getBCircle().fRadius) * 0.005f;
				Vec2 ofs = b0->dir * dist;
				for(;;) {
					ps0.setOfs(ps0.getOffset() + ofs);
					ps1.setOfs(ps1.getOffset() - ofs);

					TModelR<RPose>	tm0(*r0.getModel().cref(), ps0),
									tm1(*r1.getModel().cref(), ps1);
					epa = boost::in_place(tm0, tm1);
					if(!epa->getResult() && epa->getNearestPair().first.distance(epa->getNearestPair().second) > 1e-5f)
						break;
				}
			}
			const Vec2x2& npair = epa->getNearestPair();
			// 中点からhR0方向に伸びる直線が接触平面
			b0->pos = inner;
			b0->dir = (npair.first - npair.second).normalization();

			try {
				StLineCore st(b0->pos, b0->dir);
				auto f = CalcForce(r0, r1, inner, _coeff, st);
				_cresult.setFlagWithInfo(id0, id1, f);
			} catch(const std::exception& e) {
				std::cout << "exception occurred: " << e.what() << std::endl;
			}
		}

		RigidMgr::ERig::ERig(HRig hR): hRig(hR), prePose(hR.cref()->getPose()), nHitCur(0), hitID(~0) {}
		RigidMgr::ERig::operator const Rigid& () const {
			return *hRig.cref().get();
		}
		RigidMgr::ERig::operator Rigid& () {
			return *hRig.ref().get();
		}
		RigidMgr::ERig::operator HRig() const {
			return hRig;
		}
		void RigidMgr::ERig::resetHitCount() {
			nHitPrev = nHitCur;
			nHitCur = nHitNew = 0;
		}
		void RigidMgr::ERig::savePose() {
			prePose = static_cast<const spn::Pose2D&>(hRig.cref()->getPose());
		}
		void RigidMgr::ERig::removeOld(NP_Pool& pool) {
			int nRemove = std::max(nHitPrev - nHitCur - nHitNew, 0);
			if(nRemove > 0) {
				int n = nHitCur+nHitNew;
				NPID id = hitID;
				if(n==0)
					hitID = invalid;
				else {
					NPair* np;
					while(--n >= 0) {
						np = &pool.get(id);
						id = np->nextID;
					}
					np->nextID = invalid;
				}
				// 末尾まで全部削除
				NPID nid;
				do {
					NPair* np = &pool.get(id);
					nid = np->nextID;
					pool.rem(id);
				} while((id=nid) != invalid);
			}
		}
		spn::WithBool<RigidMgr::NPair&> RigidMgr::ERig::addHit(HRig hR, NP_Pool& pool) {
			NPID id = hitID;
			while(id != invalid) {
				NPair* np = &pool.get(id);
				if(np->hAnother == hR) {
					++np->nFrame;
					++nHitCur;
					if(hitID != id) {
						// 一旦リストから解除して・・・
						if(np->prevID != invalid)
							pool.get(np->prevID).nextID = np->nextID;
						if(np->nextID != invalid)
							pool.get(np->nextID).prevID = np->prevID;
						// 先頭に加える
						NPair* npTop = &pool.get(hitID);
						npTop->prevID = id;
						np->nextID = hitID;
						np->prevID = invalid;
						hitID = id;
					}
					return spn::make_wb(false, *np);
				}
				id = np->nextID;
			}
			++nHitNew;
			// 該当のエントリが存在しなかったので新しく作る
			auto ne = pool.alloc();
			id = ne.first;
			NPair& nnp = ne.second;
			if(hitID != invalid) {
				NPair* np = &pool.get(hitID);
				np->prevID = id;
				nnp.nextID = hitID;
			} else
				nnp.nextID = invalid;
			nnp.hAnother = hR;
			nnp.prevID = invalid;
			nnp.nFrame = 0;
			hitID = id;
			return spn::make_wb(true, nnp);
		}

		// ----- 領域の積分計算関数 -----
		namespace {
			class Tmp {
				public:
					struct TmpIn {
						float height;		// 深度
						float velN, velT;	// 衝突断面の相対速度の垂直、水平成分
						float pos[2];		// 直線に射影した2D頂点は1次元の数値になる
						Vec2	posv[2];
						float forceN;		// 点に働く抗力
					};

				private:
					TmpIn	_tmp[2];
					int		_swI = 0;
					const RPose			&_rp0,
										&_rp1;
					const RCoeff		&_coeff;
					StLineCore			_nml;
					Vec2				_div;
					DualRForce			_force = {};

					void _advance(const Vec2& p) {
						_doSwitch();
						auto& cur = _current();
						cur.height = std::fabs(_nml.dir.dot(p-_nml.pos));

						// 物体Bから見た物体Aの相対速度
						Vec2 vel = _rp0.getVelocAt(p) - _rp1.getVelocAt(p);
						float dt = _nml.dir.dot(vel);
						cur.velN = -std::min(dt, 0.f);
						cur.velT = _div.dot(vel.normalization());
						// 物体Aの重心からの相対座標 (直線方向に対して)
						cur.pos[0] = _div.dot(p - _rp0.getOffset());
						cur.pos[1] = _div.dot(p - _rp1.getOffset());
						cur.posv[0] = p;
						cur.posv[1] = p;
						cur.forceN = 0;
					}
					void _doSwitch() { _swI ^= 1; }
					TmpIn& _current() { return _tmp[_swI]; }
					TmpIn& _prev() { return _tmp[_swI^1]; }

					void _calcForce(const ConvexCore& c, float sign0) {
						const auto& pts = c.point;
						int nV = pts.size();
						if(nV < 3)
							return;

						// [0]=Aに働く力、[1]=Bに働く力
						float p_lin = 0,
								p_tor[2] = {},
								p_fdLin = 0,
								p_fdTor[2] = {};

						_advance(pts[0]);
						for(int i=1 ; i<=nV ; i++) {
							int idx = spn::CndSub(i,nV);
							_advance(pts[idx]);
							auto& cur = _current();
							auto& pre = _prev();
							float area = (cur.pos[0] - pre.pos[0]) * sign0;		// 線分の距離(面積) マイナスの場合も有り得る
							float d_area = std::fabs(cur.pos[0] - pre.pos[0]);
							// ---- calc spring ----
							// forceN(spring) = average(h0,h1) * area * spring_coeff
							cur.forceN = (pre.height + cur.height) * 0.5f * area;
							p_lin +=  cur.forceN * _coeff.spring;
							// torqueN(spring) = 1/3 * (p1h1 + (p1h2)/2 + (h1p2)/2 + p2h2) * area * spring_coeff
							for(int j=0 ; j<2 ; j++)
								p_tor[j] -= (1.f/3) * (pre.pos[j]*pre.height + (pre.pos[j]*cur.height + cur.pos[j]*pre.height)*0.5f + cur.pos[j]*cur.height) * area * _coeff.spring;

							// ---- calc dumper ----
							// forceN(dumper) = average(vn0, vn1) * d_area * dump_coeff
							float tmp = (pre.velN + cur.velN) * 0.5f * d_area;
							cur.forceN += tmp * _coeff.dumper;
							p_lin += tmp * _coeff.dumper;
							// torqueN(dumper) = 1/3 * (p1vn1 + p1vn2/2 + p2vn1/2 + p2vn2) * area * dump_coeff
							const RPose* rpp[2] = {&_rp0, &_rp1};
							for(int j=0 ; j<2 ; j++)
								p_tor[j] -= (1.f/3) * (pre.pos[j]*pre.velN + (pre.pos[j]*cur.velN + cur.pos[j]*pre.velN)*0.5f + cur.pos[j]*cur.velN) * d_area * _coeff.dumper;

							// ---- calc dynamic-friction ----
							// forceN(fricD) = average(fd0, fd1) * d_area * fricD_coeff
							float fd[2] = {pre.velT * cur.forceN,
											cur.velT * cur.forceN};
							p_fdLin += (fd[0] + fd[1]) * 0.5f * _coeff.fricD;
							// torqueN(fricD)
							for(int j=0 ; j<2 ; j++) {
								float dr[2] = {_nml.dir.cw(rpp[j]->getOffset() - pre.posv[j]) > 0 ? 1.f : -1.f,
												_nml.dir.cw(rpp[j]->getOffset() - cur.posv[j]) > 0 ? 1.f : -1.f};
								p_fdTor[j] += (1.f/3) * (pre.pos[j]*fd[0]*dr[0] + (pre.pos[j]*fd[1]*dr[1] + cur.pos[j]*fd[0]*dr[0])*0.5f + cur.pos[j]*fd[1]*dr[1]) * d_area * _coeff.fricD * sign0;
							}
							// TODO: calc static-friction
						}
						constexpr float sign[2] = {1,-1};
						for(int i=0 ; i<2 ; i++) {
							auto& fc = _force(i);

							fc.sdump.linear += _nml.dir * p_lin * sign[i];
							fc.sdump.torque +=  p_tor[i] * sign[i];
							fc.fricD.linear += _div * p_fdLin * sign[i];
							fc.fricD.torque += p_fdTor[i] * 0.1f;
						}
					}

				public:
					Tmp(const RPose &r0, const RPose &r1, const ConvexCore& cv, const StLineCore& nml, const RCoeff &coeff): _swI(0), _rp0(r0), _rp1(r1), _coeff(coeff), _nml(nml) {
						// 衝突ライン(2D)の法線
						_div = nml.dir * cs_mRot90[0];
						_nml.pos += _nml.dir * -10;
						_calcForce(cv, 1.f);
					}
					const DualRForce& getForce() const { return _force; }
			};
			//! 凸包を三角形に分割して抗力を計算
			DualRForce CalcRF_Convex(const RPose& rp0, const RPose& rp1, const RCoeff& coeff, const ConvexCore& cv, const StLineCore& div) {
				// 直線方向に向かって左側が表で、右が裏
				// st.dot(line)がプラスなら足し、逆なら引く
				Tmp tmp(rp0, rp1, cv, div, coeff);
				return tmp.getForce();
			}
			//! 円同士
			DualRForce CalcOV_Circle2(const Rigid& r0, const Rigid& r1, const Vec2& inner, const RCoeff& coeff, const StLineCore& div) {
				// 境界線分を計算して共通領域を2つに分けてそれぞれ積分
				// 中身がCircleだと分かっているのでポインタの読み替え
				AssertT(Trap, false, (std::domain_error)(const char*), "not implemented yet")
			}
			DualRForce CalcOV_Convex2(const Rigid& r0, const Rigid& r1, const Vec2& inner, const RCoeff& coeff, const StLineCore& div) {
				// 領域算出
				Convex cnv = Convex::GetOverlappingConvex(r0, r1, inner);
				cnv.adjustLoop();
				if(cnv.point.empty())
					return DualRForce();
				return CalcRF_Convex(r0, r1, coeff, cnv, div);
			}
			//! 円とBox含む多角形
			DualRForce CalcOV_CircleConvex[[noreturn]](const Rigid& r0, const Rigid& r1, const Vec2& inner, const RCoeff& coeff, const StLineCore& div) {
				PointL _pts;
				// pts[0]とpts[nV-1]の間は円弧を表す
				// 円弧部分は弓部で分けて凸包部分は三角形で計算し、残りは独自式で積分
				AssertT(Trap, false, (std::domain_error)(const char*), "not implemented yet")
			}
		}

		DualRForce CalcForce(const Rigid& r0, const Rigid& r1, const Vec2& inner, const RCoeff& coeff, const StLineCore& div) {
			// 実質Convex, Box, Circle専用
			if(r0.getCID() == Circle::GetCID()) {
				if(r1.getCID() == Circle::GetCID())
					return CalcOV_Circle2(r0, r1, inner, coeff, div);
				return CalcOV_CircleConvex(r0, r1, inner, coeff, div);
			}
			return CalcOV_Convex2(r0, r1, inner, coeff, div);
		}
	}
}