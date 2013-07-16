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

		void Rigid::addR(const SPResist& sp) {
			for(int i=0 ; i<4 ; i++) {
				if(!_resist[i]) {
					_resist[i] = sp;
					return;
				}
			}
			assert(false);
		}
		RForce::F Rigid::resist(int index, const RigidCR& cr) const {
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
			void Gravity::resist(RForce::F& acc, const Rigid& /*r*/, int /*index*/, const RigidCR& /*cr*/) const {
				acc.linear += _grav;
			}

			void Impact::resist(RForce::F& acc, const Rigid& r, int index, const RigidCR& cr) const {
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
			void Eular::advance(int pass, RItr itr, RItr itrE, const RigidCR& cr, float dt) {
				for(int i=0 ; itr!=itrE ; ++i,++itr) {
					auto* ptr = (*itr)->get();
					auto st = ptr->refPose().refValue();
					// 現フレームの加速度
					auto acc = ptr->resist(i, cr);
					st.acc += acc.linear / ptr->getArea(false);
					st.rotAcc += acc.torque / ptr->getInertia(false);
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
			void ImpEular::advance(int pass, RItr itr, RItr itrE, const RigidCR& cr, float dt) {
				auto *tv0 = &_tvalue[0];
				float dth2 = dt/2;
				int cur = 0;
				if(pass == 0) {
					// value = 1つ前の(計算上の)状態
					while(itr != itrE) {
						auto* ptr = (*itr)->get();
						auto dat = ptr->refPose().refValue();
						auto& ps = tv0[cur];

						// 内部のメモリに書き込むと同時に出力
						ps = dat;
						auto acc = ptr->resist(cur, cr);
						dat.ofs += dat.vel * dt;
						dat.vel += dat.acc * dt;
						dat.ang += dat.rotVel * dt;
						dat.rotVel += dat.rotAcc * dt;
						// (衝突判定結果は1フレーム遅れて出る為)
						ps.acc = acc.linear / ptr->getArea(false);
						ps.rotAcc = acc.torque / ptr->getInertia(false);

						++itr;
						++cur;
					}
				} else {
					while(itr != itrE) {
						auto* ptr = (*itr)->get();
						auto dat = ptr->refPose().refValue();
						auto& ps0 = tv0[cur];

						dat.ofs = ps0.ofs + (ps0.vel + dat.vel) * dth2;
						dat.vel = ps0.vel + (ps0.acc + dat.acc) * dth2;
						dat.ang = ps0.ang + (ps0.rotVel + dat.rotVel) * dth2;
						dat.rotVel = ps0.rotVel + (ps0.rotAcc + dat.rotAcc) * dth2;
						auto acc = ptr->resist(cur, cr);
						dat.acc = acc.linear / ptr->getArea(false);
						dat.rotAcc = acc.torque / ptr->getInertia(false);

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
			void RK4::advance(int pass, RItr itr, RItr itrE, const RigidCR& cr, float dt) {
				float dt2 = dt/2,
						dt6 = dt/6;
				int cur = 0,
					nR = itrE-itr;
				auto *tv0 = &_tvalue[0],
					*tv1 = &_tvalue[nR],
					*tv2 = &_tvalue[2*nR],
					*tv3 = &_tvalue[3*nR];

				switch(pass) {
					case 0:
						while(itr != itrE) {
							auto* ptr = (*itr)->get();
							auto dat = ptr->refPose().refValue();
							auto& ps = tv0[cur];

							ps = dat;
							auto acc = ptr->resist(cur, cr);		// 処理前の加速度
							ps.acc = acc.linear / ptr->getArea(false);
							ps.rotAcc = acc.torque / ptr->getInertia(false);
							dat.ofs += dat.vel * dt2;
							dat.vel += ps.acc * dt2;
							dat.ang += dat.rotVel * dt2;
							dat.rotVel += ps.rotAcc * dt2;

							++itr;
							++cur;
						}
						break;
					case 1:
						while(itr != itrE) {
							auto* ptr = (*itr)->get();
							auto dat = ptr->refPose().refValue();
							auto &ps0 = tv0[cur],
								&ps1 = tv1[cur];

							ps1 = dat;
							auto acc = ptr->resist(cur, cr);		// ps1の加速度
							ps1.acc = acc.linear / ptr->getArea(false);
							ps1.rotAcc = acc.torque / ptr->getInertia(false);
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
							auto* ptr = (*itr)->get();
							auto dat = ptr->refPose().refValue();
							auto &ps0 = tv0[cur],
								&ps2 = tv2[cur];

							ps2 = dat;
							auto acc = ptr->resist(cur, cr);		// ps2の加速度
							ps2.acc = acc.linear / ptr->getArea(false);
							ps2.rotAcc = acc.torque / ptr->getInertia(false);
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
							auto* ptr = (*itr)->get();
							auto dat = ptr->refPose().refValue();
							auto &ps0 = tv0[cur],
								&ps1 = tv1[cur],
								&ps2 = tv2[cur],
								&ps3 = tv3[cur];

							ps3 = dat;
							auto acc = ptr->resist(cur, cr);
							ps3.acc = acc.linear / ptr->getArea(false);
							ps3.rotAcc = acc.torque / ptr->getInertia(false);

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
		RigidMgr::RigidMgr(IItg::csptr itg, const RCoeff& coeff): RigidCR(coeff), _itg(itg) {}
		void RigidMgr::_checkCollision() {
			RigidCR::clear();
			RigidCR::setNumObjects(_broadC.getNumObj());
			_broadC.collision(*this);
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
		void RigidMgr::simulate(float dt) {
			IItg* itg = _itg.get();
			int nR = _broadC.getNumObj(),
				nItr = itg->numOfIteration();

			itg->beginIteration(nR);
			for(int i=0 ; i<nItr ; i++) {
				// 当たり判定結果は一回分だけとっておけば良い
				_checkCollision();
				for(int j=0 ; j<BroadC::NUM_TYPE ; j++) {
					auto typ = static_cast<BroadC::TYPE>(j);
					itg->advance(i, _broadC.cbegin(typ), _broadC.cend(typ), *this, dt);
				}
			}
			itg->endIteration();
		}

		// ----- 領域の積分計算関数 -----
		namespace {
			class Tmp {
				public:
					struct TmpIn {
						float height;		// 深度
						float velN, velT;	// 衝突断面の相対速度の垂直、水平成分
						float pos;			// 直線に射影した2D頂点は1次元の数値になる
						float fricD;
					};

				private:
					TmpIn	_tmp[2];
					int		_swI = 0;
					const RPose			&_rp0,
					&_rp1;
					const StLineCore	&_div;
					const RCoeff		&_coeff;
					Vec2				_nml;

					void _advance(const Vec2& p) {
						_doSwitch();
						auto& cur = _current();
						cur.height = std::fabs(_div.dir.ccw(p-_div.pos));

						// 物体Bから見た物体Aの相対速度
						Vec2 vel = _rp1.getVelocAt(p) - _rp0.getVelocAt(p);
						cur.velN = _nml.dot(vel);
						cur.velT = _div.dir.dot(vel);
						// 物体Aの重心からの相対座標 (直線方向に対して)
						cur.pos = _div.dir.dot(p) - _div.dir.dot(_rp0.getOffset());
						cur.fricD = 0;
					}
					void _doSwitch() { _swI ^= 1; }
					TmpIn& _current() { return _tmp[_swI]; }
					TmpIn& _prev() { return _tmp[_swI^1]; }

				public:
					Tmp(const RPose& r0, const RPose& r1, const StLineCore& div, const RCoeff& coeff): _swI(0), _rp0(r0), _rp1(r1), _div(div), _coeff(coeff) {
						// 衝突ライン(2D)の法線
						_nml = div.dir * cs_mRot90[0];
						if(_nml.dot(_rp0.getOffset() - div.pos) > 0)
							_nml *= -1.f;
					}
					void calcForce(RForce& dst, const ConvexCore& c, float sign0) {
						const auto& pts = c.point;
						int nV = pts.size();
						if(nV < 3)
							return;

						float p_lin = 0,
						p_tor = 0,
						p_fdLin = 0,
						p_fdTor = 0;
						_advance(pts[0]);
						for(int i=1 ; i<=nV ; i++) {
							int idx = spn::CndSub(i,nV);
							_advance(pts[idx]);
							auto& cur = _current();
							auto& pre = _prev();
							float area = (cur.pos - pre.pos) * sign0;		// マイナスの場合も有り得る
							// calc spring
							cur.fricD = (pre.height + cur.height) * 0.5f * area * _coeff.spring;
							p_tor += (-1.f/3) * (pre.pos*pre.height + (pre.pos*cur.height + cur.pos*pre.height)*0.5f + cur.pos*cur.height) * area * _coeff.spring;
							// calc dumper
							cur.fricD += (pre.velN + cur.velN) * 0.5f * area * _coeff.dumper;
							p_lin += cur.fricD;
							p_tor += (-1.f/3) * (pre.pos*pre.velN + (pre.pos*cur.velN + cur.pos*pre.velN)*0.5f + cur.pos*cur.velN) * area * _coeff.dumper;
							// dynamic-friction
							cur.fricD = (pre.velT + cur.velT) * 0.5f * cur.fricD;
							cur.fricD *= _coeff.fricD;
							p_fdLin += cur.fricD;
							p_fdTor += (-1.f/3) * (pre.pos*pre.fricD + (pre.pos*cur.fricD + cur.pos*pre.fricD)*0.5f + cur.pos*cur.fricD) * area * _coeff.fricD;
							// TODO: calc static-friction
						}
						dst.sdump.linear += _nml * p_lin;
						dst.sdump.torque += p_tor;
						dst.fricD.linear += _div.dir * p_fdLin;
						dst.fricD.torque += p_fdTor;
					}
			};
			//! 凸包を三角形に分割して抗力を計算
			RForce CalcRF_Convex(const RPose& rp0, const RPose& rp1, const RCoeff& coeff, const ConvexCore& cv, const StLineCore& div) {
				// 直線方向に向かって左側が表で、右が裏
				// st.dot(line)がプラスなら足し、逆なら引く
				Tmp tmp(rp0, rp1, div, coeff);
				RForce rf = {};
				// 直線(断面)で2つに切ってそれぞれ計算
				auto cvt = cv.splitTwo(div);
				tmp.calcForce(rf, cvt.first, 1.f);
				tmp.calcForce(rf, cvt.second, -1.f);
				return rf;
			}
			//! 円同士
			RForce CalcOV_Circle2(const Rigid& r0, const Rigid& r1, const Vec2& inner, const RCoeff& coeff, const StLineCore& div) {
				// 境界線分を計算して共通領域を2つに分けてそれぞれ積分
				// 中身がCircleだと分かっているのでポインタの読み替え
				throw std::runtime_error("not implemented yet");
			}
			RForce CalcOV_Convex2(const Rigid& r0, const Rigid& r1, const Vec2& inner, const RCoeff& coeff, const StLineCore& div) {
				// 領域算出
				Convex cnv = Convex::GetOverlappingConvex(r0, r1, inner);
				std::cout << "OverlappingConvex:" << std::endl << cnv << std::endl;
				return CalcRF_Convex(r0.getPose(), r1.getPose(), coeff, cnv, div);
			}
			//! 円とBox含む多角形
			RForce CalcOV_CircleConvex(const Rigid& r0, const Rigid& r1, const Vec2& inner, const RCoeff& coeff, const StLineCore& div) {
				PointL _pts;
				// pts[0]とpts[nV-1]の間は円弧を表す
				// 円弧部分は弓部で分けて凸包部分は三角形で計算し、残りは独自式で積分
				throw std::runtime_error("not implemented yet");
			}
		}

		RForce CalcForce(const Rigid& r0, const Rigid& r1, const Vec2& inner, const RCoeff& coeff, const StLineCore& div) {
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