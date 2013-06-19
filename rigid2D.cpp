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

		// -------------------------- IResist --------------------------
		void IResist::addNext(csptr sp) {
			if(_spNext)
				_spNext->addNext(sp);
			else
				_spNext = sp;
		}
		IResist::Accel IResist::resist(const RPose::Value& pose) const {
			Accel acc{Vec2(0),0};
			resist(acc, pose);
			return acc;
		}
		void IResist::_callNext(Accel& acc, const RPose::Value& pose) const {
			if(_spNext)
				_spNext->resist(acc, pose);
		}
		void IResist::resist(Accel& acc, const RPose::Value& pose) const {
			_callNext(acc, pose);
		}

		// -------------------------- Itg_Eular --------------------------
		void Itg_Eular::advance(RPose::Value& st, const IResist* pres, float dt) {
			// 次のフレームの位置
			st.ofs += st.vel * dt;
			st.ang += st.rotVel * dt;
			// 次のフレームの速度
			st.vel += st.acc * dt;
			st.rotVel += st.rotAcc * dt;
			// 次のフレームの加速度
			auto acc = pres->resist(st);
			st.acc += acc.linear;
			st.rotAcc += acc.rot;
		}

		// -------------------------- RigidMgr --------------------------
		RigidMgr::RigidMgr(IItg::csptr itg): _itg(itg), _res(new IResist()) {}
		void RigidMgr::add(Rigid::csptr sp) {
			_rlist.push_back(sp);
		}
		void RigidMgr::add(IResist::csptr sp) {
			_res->addNext(sp);
		}
		void RigidMgr::simulate(float dt) {
			IItg* itg = _itg.get();
			IResist* res = _res.get();
			for(auto& r : _rlist) {
				auto st = r->refPose().refValue();
				itg->advance(st, res, dt);
			}
		}
	}
}