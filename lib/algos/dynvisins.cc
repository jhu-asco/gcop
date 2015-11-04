#include "dynvisins.h"
#include "utils.h"

using namespace gcop;

/**
 * Standard stereo camera residual error
 */
struct StereoError {
  StereoError(const DynVisIns &vi, const Vector3d &z)
    : vi(vi), z(z) {}
  
  /**
   * @param o_ 3-dim rotation exp coordinates
   * @param p_ 3-dim position
   * @param l 3-dim feature position in 3d
   * @param res 2-dim vector residual
   */
  bool operator()(const double* o_,
                  const double* p_,
                  const double* l_,
                  double* res) const 
  {
    Matrix3d dR;
    if(vi.useCay) 
    {
      SO3::Instance().cay(dR, Vector3d(o_));
    }
    else
    {
      SO3::Instance().exp(dR, Vector3d(o_));
    }
    Matrix3d R = dR*vi.Ric;  // camera rotation

    Vector3d p(p_);
    Vector3d r = R.transpose()*(Vector3d(l_) - p);

    Vector3d e = (r - z)/vi.stereoStd;

    res[0] = e[0]; res[1] = e[1]; res[2] = e[2];
    return true;
  }
  
  static ceres::CostFunction* Create(const DynVisIns &vi, const Vector3d &z) {
    return (new ceres::NumericDiffCostFunction<StereoError, ceres::CENTRAL, 3, 3, 3, 3>(
                                                                                       new StereoError(vi, z)));
  }
  
  const DynVisIns &vi;
  Vector3d z;
};
/**
 * Standard perspective projection residual error
 */
struct PerspError {
  PerspError(const DynVisIns &vi, const Vector2d &z)
    : vi(vi), z(z) {}
  
  /**
   * @param o_ 3-dim rotation exp coordinates
   * @param p_ 3-dim position
   * @param l 3-dim feature position in 3d
   * @param res 2-dim vector residual
   */
  bool operator()(const double* o_,
                  const double* p_,
                  const double* l_,
                  double* res) const 
  {
    Matrix3d dR;
    if(vi.useCay) 
    {
      SO3::Instance().cay(dR, Vector3d(o_));
    }
    else
    {
      SO3::Instance().exp(dR, Vector3d(o_));
    }
    Matrix3d R = dR*vi.Ric;  // camera rotation

    Vector3d p(p_);
    Vector3d r = R.transpose()*(Vector3d(l_) - p);

    Vector2d e = (vi.K*(r/r[2]) - z)/vi.pxStd;
    //std::cout << "persp e=" << e.transpose() << std::endl;
    //std::cout << "persp z=" << z.transpose() << std::endl;
    //std::cout << "persp r=" << r.transpose() << std::endl;
    //std::cout << "persp l=" << Vector3d(l_).transpose() << std::endl;
    //std::cout << "persp p=" << p.transpose() << std::endl;

    res[0] = e[0]; res[1] = e[1];
    return true;
  }
  
  static ceres::CostFunction* Create(const DynVisIns &vi, const Vector2d &z) {
    return (new ceres::NumericDiffCostFunction<PerspError, ceres::CENTRAL, 2, 3, 3, 3>(
                                                                                       new PerspError(vi, z)));
  }
  
  const DynVisIns &vi;
  Vector2d z;
};

/**
 * Standard perspective projection residual error.  Uses analytic Jacobian.
 */
struct AnalyticPerspError : public ceres::SizedCostFunction<2, 3, 3, 3>{
  AnalyticPerspError(const DynVisIns &vi, const Vector2d &z)
    : vi(vi), z(z) {}
  
  virtual ~AnalyticPerspError() {}

  /**
   * @param parameters parameters to cost function: 3-dim rotation exp coordinates, 3-dim position, 3-dim feature position in 3d
   * @param res 2-dim vector residual
   * @param jacs 2x9-dim matrix jacobian
   */
  virtual bool Evaluate(double const* const* parameters,
                        double* res,
                        double** jacs) const 
  {
    Vector3d p(parameters[1]);
    Vector3d l(parameters[2]);
    
    Matrix3d dR;
    if(vi.useCay)
    {
      SO3::Instance().cay(dR, Vector3d(parameters[0]));
    }
    else
    {
      SO3::Instance().exp(dR, Vector3d(parameters[0]));
    }
    Matrix3d R = dR*vi.Ric;  // camera rotation

    Vector3d r = R.transpose()*(l - p);

    Vector2d e = (vi.K*(r/r[2]) - z)/vi.pxStd;

    res[0] = e[0]; res[1] = e[1];

    // Compute the Jacobian if asked for.
    if (jacs != NULL && jacs[0] != NULL && jacs[1] != NULL && jacs[2] != NULL) {
      const double fx = vi.K(0,0);
      const double fy = vi.K(1,1);
      const double fx_div_r2 = fx/r[2];
      const double fy_div_r2 = fy/r[2];
      const double r2sq = r[2]*r[2];
      const double fxr0_div_r2sq = fx*r[0]/r2sq;
      const double fyr1_div_r2sq = fy*r[1]/r2sq;
      Vector3d Rx = R.col(0);
      Vector3d Ry = R.col(1);
      Vector3d Rz = R.col(2);
   
      // de/dl
      Vector3d dexdl = (fx_div_r2*Rx - fxr0_div_r2sq * Rz)/vi.pxStd;  
      Vector3d deydl = (fy_div_r2*Ry - fyr1_div_r2sq * Rz)/vi.pxStd;  
      jacs[2][0] = dexdl[0]; jacs[2][1] = dexdl[1]; jacs[2][2] = dexdl[2];
      jacs[2][3] = deydl[0]; jacs[2][4] = deydl[1]; jacs[2][5] = deydl[2];

      //de_x/dp is -de_x/dl
      jacs[1][0] = -dexdl[0]; jacs[1][1] = -dexdl[1]; jacs[1][2] = -dexdl[2];
      jacs[1][3] = -deydl[0]; jacs[1][4] = -deydl[1]; jacs[1][5] = -deydl[2];
 
      // de/dR
      Matrix3d drhat;
      SO3::Instance().hat(drhat, dR.transpose()*(l-p));
      Matrix3d drdR = vi.Ric.transpose()*drhat;

      Vector3d dexdR = (fx_div_r2*drdR.row(0) - fxr0_div_r2sq * drdR.row(2))/vi.pxStd; 
      Vector3d deydR = (fy_div_r2*drdR.row(1) - fyr1_div_r2sq * drdR.row(2))/vi.pxStd; 
      jacs[0][0] = dexdR[0]; jacs[0][1] = dexdR[1]; jacs[0][2] = dexdR[2];
      jacs[0][3] = deydR[0]; jacs[0][4] = deydR[1]; jacs[0][5] = deydR[2];
    }
    return true;
  }

  const DynVisIns &vi;
  Vector2d z;
};

/**
 * Unit-spherical projection residual error, 
 * with simplified constant spherical unprojected covariance
 */
struct SphError {
  SphError(const DynVisIns &vi, const Vector3d &z)
    : vi(vi), z(z) {}
  
  /**
   * @param o_ 3-dim rotation exp coordinates
   * @param p_ 3-dim position
   * @param l_ 3-dim feature position in 3d
   * @param res 3-dim vector residual
   */
  bool operator()(const double* o_,
                  const double* p_,                  
                  const double* l_,
                  double* res) const 
  {
    Matrix3d dR;  
    if(vi.useCay)
    {
      SO3::Instance().cay(dR, Vector3d(o_));
    }
    else
    {
      SO3::Instance().exp(dR, Vector3d(o_));
    }
    Matrix3d R = dR*vi.Ric;  // camera rotation

    Vector3d r = R.transpose()*(Vector3d(l_) - Vector3d(p_));
    Vector3d e = (r/r.norm()- z)/vi.sphStd;
    //std::cout << "sph e=" << e.transpose() << std::endl;
    res[0] = e[0]; res[1] = e[1]; res[2] = e[2];
    return true;
  }
  
  static ceres::CostFunction* Create(const DynVisIns &vi, const Vector3d &z) {
    return (new ceres::NumericDiffCostFunction<SphError, ceres::CENTRAL, 3, 3, 3, 3>(
                                                                                  new SphError(vi, z)));
  }
  const DynVisIns &vi;
  Vector3d z;
};


/**
 * A basic cubic interpolator
 *
 * Author: Marin Kobilarov
 */
class Cubic {  
public:
  Cubic(const Vector3d &p0, const Vector3d &v0,
        const Vector3d &p1, const Vector3d &v1,
        double dt) : p0(p0), v0(v0), dt(dt) {
    double dt2 = dt*dt;
    Vector3d q1 = p1 - p0 - dt*v0;
    Vector3d q2 = v1 - v0;
    b = 6/dt2*q1 + -2/dt*q2;
    c = -6/(dt2*dt)*q1 + 3/dt2*q2;    
  }
  
  Cubic(const Vector3d &p0, const Vector3d &v0,
        double dt) : p0(p0), v0(v0), dt(dt) {
    b.setZero();
    c.setZero();
  }
  

  bool GetPos(Vector3d &p, double t) {
    if (t > dt)
      return false;
    double t2 = t*t;
    p = p0 + t*v0 + t2/2*b + t2*t/3*c;
    return true;
  }

  bool GetVel(Vector3d &v, double t) {
    if (t > dt)
      return false;
    v = v0 + t*b + t*t*c;
    return true;
  }

  bool GetAcc(Vector3d &a, double t) {
    if (t > dt)
      return false;
    a = b + 2*t*c;
    return true;
  }

  void GetDPosDp0(Matrix3d& DpDp0, double t)
  {
    Matrix3d I3 = Eigen::Matrix3d::Identity(3,3);
    double dt2 = dt*dt;
    DpDp0 = (1+ (t*t/2)*(-6/dt2) + (t*t*t/3)*(-6/(dt2*dt)*(-1)))*I3;
  }

  void GetDPosDp1(Matrix3d& DvDp1, double t)
  {
    Matrix3d I3 = Eigen::Matrix3d::Identity(3,3);
    double dt2 = dt*dt;
    DvDp1 = ((t*t/2)*(6/dt2) + (t*t*t/3)*(-6/(dt2*dt)))*I3;
  }

  void GetDPosDv0(Matrix3d& DvDv0, double t)
  {
    Matrix3d I3 = Eigen::Matrix3d::Identity(3,3);
    double dt2 = dt*dt;
    DvDv0 = (t + (t*t/2)*(6/dt2*(-dt)+(2/dt)) + (t*t*t/3)*(-6/(dt2*dt)*(-dt)+ 3/dt2*(-1) ))*I3;
  }

  void GetDPosDv1(Matrix3d& DvDv1, double t)
  {
    Matrix3d I3 = Eigen::Matrix3d::Identity(3,3);
    double dt2 = dt*dt;
    DvDv1 = ((t*t/2)*(-2/dt) + (t*t*t/3)*(3/dt2))*I3;
  }
 
  void GetDVelDp0(Matrix3d& DvDp0, double t)
  {
    Matrix3d I3 = Eigen::Matrix3d::Identity(3,3);
    double dt2 = dt*dt;
    DvDp0 = (-6*t/dt2 + t*t*(-6/(dt2*dt)*(-1)))*I3;
  }

  void GetDVelDp1(Matrix3d& DvDp1, double t)
  {
    Matrix3d I3 = Eigen::Matrix3d::Identity(3,3);
    double dt2 = dt*dt;
    DvDp1 = (6*t/dt2 + t*t*(-6/(dt2*dt)))*I3;
  }

  void GetDVelDv0(Matrix3d& DvDv0, double t)
  {
    Matrix3d I3 = Eigen::Matrix3d::Identity(3,3);
    double dt2 = dt*dt;
    DvDv0 = (1 + t*(6/dt2*(-dt)+(2/dt)) + t*t*(-6/(dt2*dt)*(-dt)+ 3/dt2*(-1) ))*I3;
  }

  void GetDVelDv1(Matrix3d& DvDv1, double t)
  {
    Matrix3d I3 = Eigen::Matrix3d::Identity(3,3);
    double dt2 = dt*dt;
    DvDv1 = (t*(-2/dt) + t*t*(3/dt2))*I3;
  }

  void GetDAccDp0(Matrix3d& DaDp0, double t)
  {
    Matrix3d I3 = Eigen::Matrix3d::Identity(3,3);
    double dt2 = dt*dt;
    DaDp0 = (-6/dt2 + 2*t*(-6/(dt2*dt)*(-1)))*I3;
  }

  void GetDAccDp1(Matrix3d& DaDp1, double t)
  {
    Matrix3d I3 = Eigen::Matrix3d::Identity(3,3);
    double dt2 = dt*dt;
    DaDp1 = (6/dt2 + 2*t*(-6/(dt2*dt)))*I3;
  }

  void GetDAccDv0(Matrix3d& DaDv0, double t)
  {
    Matrix3d I3 = Eigen::Matrix3d::Identity(3,3);
    double dt2 = dt*dt;
    DaDv0 = ((6/dt2*(-dt)+(2/dt)) + 2*t*(-6/(dt2*dt)*(-dt)+ 3/dt2*(-1) ))*I3;
  }

  void GetDAccDv1(Matrix3d& DaDv1, double t)
  {
    Matrix3d I3 = Eigen::Matrix3d::Identity(3,3);
    double dt2 = dt*dt;
    DaDv1 = ((-2/dt) + 2*t*(3/dt2))*I3;
  }
  
  /**
   * Get body-fixed angular velocity at time t, assuming exponential parametrization
   * @param body-fixed angular velocity
   * @param t time
   * @return true if success
   */
  bool GetExpVel(Vector3d &w, double t) {
    if (!GetVel(w, t))
      return false;
    Matrix3d D;
    SO3::Instance().dexp(D, -w);
    w = D*w;
    return true;
  }
  
  Vector3d p0, v0;
  double dt;
  Vector3d b, c;
};
        

/**
 * Gyro error, assuming segment is parametrized as a cubic spline
 */
struct GyroCubicError {
  /**
   * @param vi 
   * @param dt delta-t for this segment
   * @param ts local times for each measurement
   * @param ws the sequence of measurements in this segment
   */
  GyroCubicError(const DynVisIns &vi, 
                 double dt,
                 const vector<double> &ts,
                 const vector<Vector3d> &ws)
    : vi(vi), dt(dt), ts(ts), ws(ws)
  {
    assert(dt > 0);
    assert(ts.size() == ws.size());
  }
  
  /**
   * Computes IMU error between two states xa and xb
   * @param wa_ angular velocity at start of segment
   * @param wb_ angular velocity at end of segment
   * @param res residual for all gyro measurements within segment
   */
  bool operator()(const double* ra_,
                  const double* dra_,
                  const double* rb_,
                  const double* drb_,
                  double* res) const 
  {
    Vector3d ra(ra_);
    Vector3d dra(dra_);
    Vector3d rb(rb_);
    Vector3d drb(drb_);

    Cubic cub(ra, dra, rb, drb, dt);
    Vector3d r, dr;
    Matrix3d D;

    for (int i = 0; i < ws.size(); ++i) {      
    
      const double &t = ts[i];
      cub.GetPos(r, t);  // get exp coord
      cub.GetVel(dr, t); // get exp coord vel

      if(vi.useCay)
      {
        SO3::Instance().dcay(D, -r);
      }
      else
      {
        SO3::Instance().dexp(D, -r);
      }
      Vector3d w = D*dr; // the body-fixed angular velocity

      // for now just assume noise is spherical and defined in wStd
      Vector3d e = (w - ws[i] + vi.bg)/vi.wStd;
      //std::cout << "gyro e=" << e.transpose() << std::endl;
      memcpy(res + 3*i, e.data(), 3*sizeof(double));
    }
    return true;
  }
  
  static ceres::CostFunction* Create(const DynVisIns &vi, 
                                     double dt,
                                     const vector<double> &ts,
                                     const vector<Vector3d> &ws) {
    return new ceres::NumericDiffCostFunction<GyroCubicError, ceres::CENTRAL, ceres::DYNAMIC, 3, 3, 3, 3>(
                                                                                                          new GyroCubicError(vi, dt, ts, ws), ceres::TAKE_OWNERSHIP, ws.size()*3);
  }
  
  const DynVisIns &vi;
  double dt;                  ///< total time for this segment
  const vector<double> ts;    ///< sequence of relative times at which gyro measurements arrived
  const vector<Vector3d> ws;  ///< sequence of angular measurements  
};

/**
 * Gyro error, assuming segment is parametrized as a cubic spline.  Produces analytic jacobians.
 */
struct AnalyticGyroCubicError : public ceres::CostFunction{
  /**
   * @param vi 
   * @param dt delta-t for this segment
   * @param ts local times for each measurement
   * @param ws the sequence of measurements in this segment
   */
  AnalyticGyroCubicError(const DynVisIns &vi, 
                 double dt,
                 const vector<double> &ts,
                 const vector<Vector3d> &ws)
    : vi(vi), dt(dt), ts(ts), ws(ws)
  {
    assert(dt > 0);
    assert(ts.size() == ws.size());
    set_num_residuals(3*ws.size());
    vector<ceres::int32>* p_block_sizes = mutable_parameter_block_sizes();
    p_block_sizes->resize(4);
    for(int i = 0; i < 4; i++)
      (*p_block_sizes)[i] = 3;
  }

  /**
   * Computes IMU error between two states xa and xb
   * @param parameters exponential rotation, angular velocity at start of segment and exponential rotation, angular velocity at end of segment
   * @param res residual for all gyro measurements within segment
   * @param jacs jacobian of residuals for each parameter block
   */
  virtual bool Evaluate(double const* const* parameters,
                        double* res,
                        double** jacs) const 
  {
    Vector3d ra(parameters[0]);
    Vector3d dra(parameters[1]);
    Vector3d rb(parameters[2]);
    Vector3d drb(parameters[3]);

    Cubic cub(ra, dra, rb, drb, dt);
    Vector3d r, dr;
    Matrix3d D;

    for (int i = 0; i < ws.size(); ++i) {      
    
      const double &t = ts[i];
      cub.GetPos(r, t);  // get exp coord
      cub.GetVel(dr, t); // get exp coord vel

      if(vi.useCay)
      {
        SO3::Instance().dcay(D, -r);
      }
      else
      {
        SO3::Instance().dexp(D, -r);
      }
      Vector3d w = D*dr; // the body-fixed angular velocity

      // for now just assume noise is spherical and defined in wStd
      Vector3d e = (w - ws[i] + vi.bg)/vi.wStd;
      //std::cout << "gyro e=" << e.transpose() << std::endl;
      memcpy(res + 3*i, e.data(), 3*sizeof(double));

      if(jacs != NULL && jacs[0] != NULL && jacs[1] != NULL && jacs[2] != NULL && jacs[3] != NULL)
      {
        Eigen::Matrix3d dD_dr1, dD_dr2, dD_dr3;
        SO3::Instance().ddcay(dD_dr1, dD_dr2, dD_dr3, -r); 
        Eigen::Matrix3d dD1_dr, dD2_dr, dD3_dr;
        dD1_dr.row(0) = dD_dr1.row(0);
        dD1_dr.row(1) = dD_dr2.row(0);
        dD1_dr.row(2) = dD_dr3.row(0);

        dD2_dr.row(0) = dD_dr1.row(1);
        dD2_dr.row(1) = dD_dr2.row(1); 
        dD2_dr.row(2) = dD_dr3.row(1);

        dD3_dr.row(0) = dD_dr1.row(2);
        dD3_dr.row(1) = dD_dr2.row(2);
        dD3_dr.row(2) = dD_dr3.row(2);

        //de/dra
        {
          Eigen::Matrix3d dr_dra, ddr_dra;
          cub.GetDPosDp0(dr_dra, t);
          cub.GetDVelDp0(ddr_dra, t);

          Eigen::Matrix3d de_dra;
          de_dra.row(0) = (-dD1_dr*dr_dra*dr).transpose() + Vector3d(1,0,0).transpose()*D*ddr_dra;
          de_dra.row(1) = (-dD2_dr*dr_dra*dr).transpose() + Vector3d(0,1,0).transpose()*D*ddr_dra;
          de_dra.row(2) = (-dD3_dr*dr_dra*dr).transpose() + Vector3d(0,0,1).transpose()*D*ddr_dra;
          de_dra = de_dra/vi.wStd; 

          for(int j = 0; j < 3; j++)
          {
            for(int k = 0; k < 3; k++)
            {
              jacs[0][(3*i+j)*3+k] = de_dra(j,k); //de_dra has res on the rows and param on the cols
            }
          }
        }
        //de/ddra
        {
          Eigen::Matrix3d dr_ddra, ddr_ddra;
          cub.GetDPosDv0(dr_ddra, t);
          cub.GetDVelDv0(ddr_ddra, t);

          Eigen::Matrix3d de_ddra;
          de_ddra.row(0) = (-dD1_dr*dr_ddra*dr).transpose() + Vector3d(1,0,0).transpose()*D*ddr_ddra;
          de_ddra.row(1) = (-dD2_dr*dr_ddra*dr).transpose() + Vector3d(0,1,0).transpose()*D*ddr_ddra;
          de_ddra.row(2) = (-dD3_dr*dr_ddra*dr).transpose() + Vector3d(0,0,1).transpose()*D*ddr_ddra;
          de_ddra = de_ddra/vi.wStd; 

          for(int j = 0; j < 3; j++)
          {
            for(int k = 0; k < 3; k++)
            {
              jacs[1][(3*i+j)*3+k] = de_ddra(j,k); //de_ddra has res on the rows and param on the cols
            }
          }
        }
        //de/drb
        {
          Eigen::Matrix3d dr_drb, ddr_drb;
          cub.GetDPosDp1(dr_drb, t);
          cub.GetDVelDp1(ddr_drb, t);

          Eigen::Matrix3d de_drb;
          de_drb.row(0) = (-dD1_dr*dr_drb*dr).transpose() + Vector3d(1,0,0).transpose()*D*ddr_drb;
          de_drb.row(1) = (-dD2_dr*dr_drb*dr).transpose() + Vector3d(0,1,0).transpose()*D*ddr_drb;
          de_drb.row(2) = (-dD3_dr*dr_drb*dr).transpose() + Vector3d(0,0,1).transpose()*D*ddr_drb;
          de_drb = de_drb/vi.wStd; 

          for(int j = 0; j < 3; j++)
          {
            for(int k = 0; k < 3; k++)
            {
              jacs[2][(3*i+j)*3+k] = de_drb(j,k); //de_drb has res on the rows and param on the cols
            }
          }
        }
        //de/ddrb
        {
          Eigen::Matrix3d dr_ddrb, ddr_ddrb;
          cub.GetDPosDv1(dr_ddrb, t);
          cub.GetDVelDv1(ddr_ddrb, t);

          Eigen::Matrix3d de_ddrb;
          de_ddrb.row(0) = (-dD1_dr*dr_ddrb*dr).transpose() + Vector3d(1,0,0).transpose()*D*ddr_ddrb;
          de_ddrb.row(1) = (-dD2_dr*dr_ddrb*dr).transpose() + Vector3d(0,1,0).transpose()*D*ddr_ddrb;
          de_ddrb.row(2) = (-dD3_dr*dr_ddrb*dr).transpose() + Vector3d(0,0,1).transpose()*D*ddr_ddrb;
          de_ddrb = de_ddrb/vi.wStd; 

          for(int j = 0; j < 3; j++)
          {
            for(int k = 0; k < 3; k++)
            {
              jacs[3][(3*i+j)*3+k] = de_ddrb(j,k); //de_ddrb has res on the rows and param on the cols
            }
          }
        }
      }
    }
    return true;
  }
  
  const DynVisIns &vi;
  double dt;                  ///< total time for this segment
  const vector<double> ts;    ///< sequence of relative times at which gyro measurements arrived
  const vector<Vector3d> ws;  ///< sequence of angular measurements  
};

/**
 * Accelerometer error, assuming segment is parametrized as a cubic spline
 */
struct AccCubicError {

  /**
   * @param vi 
   * @param dt delta-t for this segment
   * @param ts local times for each measurement (should be in [0,dt])
   * @param as the sequence of measurements in this segment
   */
  AccCubicError(const DynVisIns &vi, 
                double dt,
                const vector<double> &ts,
                const vector<Vector3d> &as)
    : vi(vi), dt(dt), ts(ts), as(as)
  {
    assert(dt > 0);
    assert(ts.size() == as.size());
  }
  
  /**
   * Computes IMU error between two states xa and xb
   * @param wa_ angular velocity at start of segment
   * @param wb_ angular velocity at end of segment
   * @param res residual for all gyro measurements within segment
   */
  bool operator()(const double* ra_,
                  const double* pa_,
                  const double* dra_,
                  const double* va_,
                  const double* rb_,
                  const double* pb_,
                  const double* drb_,
                  const double* vb_,
                  double* res) const 
  {

    Vector3d ra(ra_);
    Vector3d pa(pa_);
    Vector3d dra(dra_);
    Vector3d va(va_);
    Vector3d rb(rb_);
    Vector3d pb(pb_);
    Vector3d drb(drb_);
    Vector3d vb(vb_);

    Cubic cr(ra, dra, rb, drb, dt);
    Cubic cp(pa, va, pb, vb, dt);

    Vector3d r, a;
    Matrix3d R;
    
    for (int i = 0; i < as.size(); ++i) {      
      const double &t = ts[i];

      cr.GetPos(r, t);
      if(vi.useCay)
      {
        SO3::Instance().cay(R, r);
      }
      else
      {
        SO3::Instance().exp(R, r);
      }
      cp.GetAcc(a, t);
            
      // for now just assume noise is spherical and defined in vi.aStd
      Vector3d e = (a - R*(as[i] - vi.ba) + vi.g0)/vi.aStd;
      //      cout << a.transpose() << as[i].transpose() << (R*as[i]).transpose() << e.transpose() << endl;
      //std::cout << "acc e=" << e.transpose() << std::endl;
      memcpy(res + 3*i, e.data(), 3*sizeof(double));
    }
    return true;
  }
  
  static ceres::CostFunction* Create(const DynVisIns &vi, 
                                     double dt,
                                     const vector<double> &ts,
                                     const vector<Vector3d> &as) {
    return new ceres::NumericDiffCostFunction<AccCubicError, ceres::CENTRAL, ceres::DYNAMIC, 3, 3, 3, 3, 3, 3, 3, 3>(
                                                                                                        new AccCubicError(vi, dt, ts, as), ceres::TAKE_OWNERSHIP, as.size()*3);
  }
  
  const DynVisIns &vi;
  double dt;                  ///< total time for this segment
  const vector<double> ts;    ///< sequence of relative times at which gyro measurements arrived
  const vector<Vector3d> as;  ///< sequence of angular measurements
};

/**
 * Accelerometer error, assuming segment is parametrized as a cubic spline.  
 * Provides analytic jacobians.
 */
struct AnalyticAccCubicError : public ceres::CostFunction {

  /**
   * @param vi 
   * @param dt delta-t for this segment
   * @param ts local times for each measurement (should be in [0,dt])
   * @param as the sequence of measurements in this segment
   */
  AnalyticAccCubicError(const DynVisIns &vi, 
                double dt,
                const vector<double> &ts,
                const vector<Vector3d> &as)
    : vi(vi), dt(dt), ts(ts), as(as)
  {
    assert(dt > 0);
    assert(ts.size() == as.size());
    set_num_residuals(3*as.size());
    vector<ceres::int32>* p_block_sizes = mutable_parameter_block_sizes();
    p_block_sizes->resize(8);
    for(int i = 0; i < p_block_sizes->size(); i++)
      (*p_block_sizes)[i] = 3;
  }
  
  /**
   * Computes accelrometer error between two states xa and xb.  Provides residual jacobians.
   * @param parameters ra, pa, dra, va, rb, pb, drb, vb 
   * @param res residual for all accelerometer measurements within segment
   * @param jacs jacobian for residuals with respect to parameters
   */
  virtual bool Evaluate(double const* const* parameters,
                        double* res,
                        double** jacs) const 
  {

    Vector3d ra(parameters[0]);
    Vector3d pa(parameters[1]);
    Vector3d dra(parameters[2]);
    Vector3d va(parameters[3]);
    Vector3d rb(parameters[4]);
    Vector3d pb(parameters[5]);
    Vector3d drb(parameters[6]);
    Vector3d vb(parameters[7]);

    Cubic cr(ra, dra, rb, drb, dt);
    Cubic cp(pa, va, pb, vb, dt);

    Vector3d r, a;
    Matrix3d R;
    
    for (int i = 0; i < as.size(); ++i) {      
      const double &t = ts[i];

      cr.GetPos(r, t);
      if(vi.useCay)
      {
        SO3::Instance().cay(R, r);
      }
      else
      {
        SO3::Instance().exp(R, r);
      }
      cp.GetAcc(a, t);
      Vector3d az = R*(as[i] - vi.ba);
      // for now just assume noise is spherical and defined in vi.aStd
      Vector3d e = (a - az + vi.g0)/vi.aStd;
      //      cout << a.transpose() << as[i].transpose() << (R*as[i]).transpose() << e.transpose() << endl;
      //std::cout << "acc e=" << e.transpose() << std::endl;
      memcpy(res + 3*i, e.data(), 3*sizeof(double));

      if(jacs != NULL && jacs[0] != NULL && jacs[1] != NULL && jacs[2] != NULL && jacs[3] != NULL 
        && jacs[4] != NULL && jacs[5] != NULL && jacs[6] != NULL && jacs[7] != NULL)
      {
        Matrix3d azh;
        SO3::Instance().hat(azh, az);
        Matrix3d dR;
        SO3::Instance().dcay(dR,r);
        // de/dra
        {
          Matrix3d DrDra;
          cr.GetDPosDp0(DrDra, t);
          Matrix3d de_dra = azh*dR*DrDra/vi.aStd;
          for(int j = 0; j < 3; j++)
          {
            for(int k = 0; k < 3; k++)
            {
              jacs[0][(3*i+j)*3+k] = de_dra(j,k); 
            }
          }
        }
        // de/dpa
        {
          Matrix3d DaDp0;
          cp.GetDAccDp0(DaDp0, t);
          Matrix3d de_dpa = DaDp0/vi.aStd;
          for(int j = 0; j < 3; j++)
          {
            for(int k = 0; k < 3; k++)
            {
              jacs[1][(3*i+j)*3+k] = de_dpa(j,k); 
            }
          }
        }
        
        // de/ddra
        {
          Matrix3d DrDdra;
          cr.GetDPosDv0(DrDdra, t);
          Matrix3d de_ddra = azh*dR*DrDdra/vi.aStd;
          for(int j = 0; j < 3; j++)
          {
            for(int k = 0; k < 3; k++)
            {
              jacs[2][(3*i+j)*3+k] = de_ddra(j,k); 
            }
          }
        }
        // de/dva
        {
          Matrix3d DaDv0;
          cp.GetDAccDv0(DaDv0, t);
          Matrix3d de_dva = DaDv0/vi.aStd;
          for(int j = 0; j < 3; j++)
          {
            for(int k = 0; k < 3; k++)
            {
              jacs[3][(3*i+j)*3+k] = de_dva(j,k); 
            }
          }
        }
        // de/drb
        {
          Matrix3d DrDrb;
          cr.GetDPosDp1(DrDrb, t);
          Matrix3d de_drb = azh*dR*DrDrb/vi.aStd;
          for(int j = 0; j < 3; j++)
          {
            for(int k = 0; k < 3; k++)
            {
              jacs[4][(3*i+j)*3+k] = de_drb(j,k); 
            }
          }
        }
        // de/dpb
        {
          Matrix3d DaDp1;
          cp.GetDAccDp1(DaDp1, t);
          Matrix3d de_dpb = DaDp1/vi.aStd;
          for(int j = 0; j < 3; j++)
          {
            for(int k = 0; k < 3; k++)
            {
              jacs[5][(3*i+j)*3+k] = de_dpb(j,k); 
            }
          }
        }
        // de/ddrb
        {
          Matrix3d DrDdrb;
          cr.GetDPosDv1(DrDdrb, t);
          Matrix3d de_ddrb = azh*dR*DrDdrb/vi.aStd;
          for(int j = 0; j < 3; j++)
          {
            for(int k = 0; k < 3; k++)
            {
              jacs[6][(3*i+j)*3+k] = de_ddrb(j,k); 
            }
          }
        }
        // de/dvb
        {
          Matrix3d DaDv1;
          cp.GetDAccDv1(DaDv1, t);
          Matrix3d de_dvb = DaDv1/vi.aStd;
          for(int j = 0; j < 3; j++)
          {
            for(int k = 0; k < 3; k++)
            {
              jacs[7][(3*i+j)*3+k] = de_dvb(j,k); 
            }
          }
        }
      
      }
    }
    return true;
  }
  
  const DynVisIns &vi;
  double dt;                  ///< total time for this segment
  const vector<double> ts;    ///< sequence of relative times at which gyro measurements arrived
  const vector<Vector3d> as;  ///< sequence of angular measurements
};


/**
 * Constant Velocity Rotational Error
 */
struct CvCubicRotError {
  /**
   * @param vi  visual-inertial estimaor
   * @param dt segment delta t
   */
  CvCubicRotError(DynVisIns &vi, 
                  double dt)
    : vi(vi), dt(dt)
  { 
    assert(dt > 0);
  }
  
  /**
   * Computes dyn error between two states
   */
  bool operator()(const double* ra_,
                  const double* dra_,
                  const double* rb_,
                  const double* drb_,
                  double* res) const 
  {
    Vector3d ra(ra_);
    Vector3d dra(dra_);
    Vector3d rb(rb_);
    Vector3d drb(drb_);

    Cubic cub(ra, dra, rb, drb, dt);

    Vector3d e1 = (sqrt(dt)/vi.dwStd)*(cub.b + dt*cub.c);
    Vector3d e2 = (sqrt(dt/3)/vi.dwStd)*(dt*cub.c);
    memcpy(res, e1.data(), 3*sizeof(double));
    memcpy(res + 3, e2.data(), 3*sizeof(double));

    return true;
  }
  
  static ceres::CostFunction* Create(DynVisIns &vi, 
                                     double dt) {
    return (new ceres::NumericDiffCostFunction<CvCubicRotError, ceres::CENTRAL, 6, 3, 3, 3, 3>(
                                                                                          new CvCubicRotError(vi, dt)));
  }
  DynVisIns &vi;
  double dt;
};



/**
 * Constant Velocity Rotational Error
 */
struct CvCubicPosError {
  /**
   * @param vi  visual-inertial estimaor
   * @param dt segment delta t
   */
  CvCubicPosError(DynVisIns &vi, 
                  double dt)
    : vi(vi), dt(dt)
  { 
    assert(dt > 0);
  }
  
  /**
   * Computes dyn error between two states
   */
  bool operator()(const double* pa_,
                  const double* va_,
                  const double* pb_,
                  const double* vb_,
                  double* res) const 
  {
    Vector3d pa(pa_);
    Vector3d va(va_);
    Vector3d pb(pb_);
    Vector3d vb(vb_);

    Cubic cub(pa, va, pb, vb, dt);

    Vector3d e1 = (sqrt(dt)/vi.dvStd)*(cub.b + dt*cub.c);
    Vector3d e2 = (sqrt(dt/3)/vi.dvStd)*(dt*cub.c);
    memcpy(res, e1.data(), 3*sizeof(double));
    memcpy(res + 3, e2.data(), 3*sizeof(double));

    return true;
  }
  
  static ceres::CostFunction* Create(DynVisIns &vi, 
                                     double dt) {
    return (new ceres::NumericDiffCostFunction<CvCubicPosError, ceres::CENTRAL, 6, 3, 3, 3, 3>(
                                                                                          new CvCubicPosError(vi, dt)));
  }
  DynVisIns &vi;
  double dt;
};



/** 
 * Prior residual on the state
 */
struct StatePrior {
  /**
   * @param vi visual-inertial estimator
   * @param x0 prior mean
   * @param P0 prior covariance
   */
  StatePrior(DynVisIns &vi, 
             const Body3dState &x0)
    : vi(vi), x0(x0)
  {    
    LLT<Matrix12d> llt(x0.P.inverse());   // assume P0>0
    this->W = llt.matrixU();   // the weight matrix W such that W'*W=inv(P0)
  }
  
  /**
   * Computes IMU state prior error
   * @param x_q 6-dim vector containing orientation and position
   * @param x_v 3-dim vector containing velocity v
   */
  bool operator()(const double* r_,
                  const double* p_,
                  const double* dr_,
                  const double* v_,
                  double* res) const 
  {
    Matrix3d R;
    Matrix3d D;

    Vector3d r(r_);
    if(vi.useCay)
    {
      SO3::Instance().cay(R, r);
    }
    else
    {
      SO3::Instance().exp(R, r);
    }

    Vector3d er;
    if(vi.useCay)
    {
      SO3::Instance().cayinv(er, x0.R.transpose()*R); 
    }
    else
    {
      SO3::Instance().log(er, x0.R.transpose()*R); 
    }
   
    Vector12d e;
    e.head<3>() = er;
    e.segment<3>(3) = Vector3d(p_) - x0.p;
    
    // @mk: this could be optimized to perform the dexpinv to the x0 prior and just 
    // take the diff in coordinates here:
    if(vi.useCay)
    {
      SO3::Instance().dcay(D, -r);
    }
    else
    {
      SO3::Instance().dexp(D, -r);
    }
    Vector3d w = D*Vector3d(dr_); // the body-fixed angular velocity

    e.segment<3>(6) = w - x0.w;
    e.tail<3>() = Vector3d(v_) - x0.v;
    e = W*e;

    //std::cout << "prior e=" << e.transpose() << std::endl;
    memcpy(res, e.data(), 12*sizeof(double));

    return true;
  }
  
  static ceres::CostFunction* Create(DynVisIns &vi, 
                                     const Body3dState &x0) {
    return (new ceres::NumericDiffCostFunction<StatePrior, ceres::CENTRAL, 12, 3, 3, 3, 3>(
                                                                                           new StatePrior(vi, x0)));
    
  }

  DynVisIns &vi;       ///< visual-inertial estimator
  Body3dState x0;      ///< prior state x0
  Matrix12d W;         ///< residual weight matrix W is such that W'*W=inv(P0)
};

/** 
 * Prior residual on a feature
 */
struct FeaturePrior {
  /**
   * @param vi visual-inertial estimator
   * @param x0 prior mean
   * @param P0 prior covariance
   */
  FeaturePrior(DynVisIns &vi, 
             const DynVisIns::Point &p)
    : vi(vi), p(p)
  {    
    LLT<Matrix3d> llt(p.P.inverse());   // assume P0>0
    this->W = llt.matrixU();   // the weight matrix W such that W'*W=inv(P0)
  }
  
  /**
   * Computes feature prior error
   * @param l_ 3-dim vector containing feature coordinates
   * @param res 3-dim vector containing residual
   */
  bool operator()(const double* l_,
                  double* res) const 
  {
    Vector3d l(l_);

    Vector3d e;
    e = W*(l-p.l);

    memcpy(res, e.data(), 3*sizeof(double));

    return true;
  }
  
  static ceres::CostFunction* Create(DynVisIns &vi, 
                                     const DynVisIns::Point &p) {
    return (new ceres::NumericDiffCostFunction<FeaturePrior, ceres::CENTRAL, 3, 3>(
                                                                                           new FeaturePrior(vi, p)));
    
  }

  DynVisIns &vi;       ///< visual-inertial estimator
  DynVisIns::Point p;      ///< prior feature position
  Matrix3d W;         ///< residual weight matrix W is such that W'*W=inv(P0)
};

DynVisIns::DynVisIns() : t(-1), tc(-1), problem(NULL), l_opti_map(NULL), n_good_pnts(0) {
  
  maxIterations = 50;
  maxCams = 0;
  minPnts = 20;
  ceresActive = false;
  useHuberLoss = true;
  checkPtActiveFlag = false;
  v = 0;
  
  // initial state / prior info
  x0.Clear();  
  x0.P.topLeftCorner<3,3>().diagonal().setConstant(.0001);  // R
  x0.P.block<3,3>(3,3).diagonal().setConstant(.0001);      // p
  x0.P.block<3,3>(6,6).diagonal().setConstant(.0001);    // w
  x0.P.block<3,3>(9,9).diagonal().setConstant(.0001);    // v
  
  // from IMU to Cam rotation: first -90Z  then -90X
  Ric << 0, 0, 1, -1, 0, 0, 0, -1, 0;      

  useImu = true;
  useCam = true;
  useDyn = true;
  usePrior = true;
  useFeatPrior = false;
  useCay = false;
  useAnalyticJacs = false;

  optBias = false;

  sphMeas = false;


  pxStd = 1;   // pixel measurement error standard deviation
  sphStd = 1;  // corresponding spherial error -- this will be set to the right value later
  
  dwStd = 1;   // angular acc white nose (rad/s^2)
  dvStd = 5;   // linear acc white noise (m/s^2)

  wStd = .001;    // gyro noise
  aStd = .02;     // acc noise

  g0 << 0, 0, 9.81;

  camId = -1;
  camId0 = 0;
}

DynVisIns::~DynVisIns()
{
}



bool DynVisIns::ProcessImu(double t, const Vector3d &w, const Vector3d &a) {

  // initialize time if this is the first IMU measurement
  if (this->t < 0) {
    this->t = t;
    return true;
  }

  
  double dt = t - this->t;

  if (dt <= 0) {
    cout <<"[W] IMU data out of sync dt=" << dt << endl;
    return false;
  }

  // for now camera is required
  assert(useCam);
      
  // accumulate measurements for this last camera segment 
  // only if a camera frame has alrady been added
  if (useCam && camId >= 0) {
    assert(t > tc);
    //    assert(tss.size());
    //    assert(wss.size());
    //    assert(ass.size());

    Camera &cam = cams[camId];

    cam.ts.push_back(t - tc); // add local time
    cam.ws.push_back(w);
    cam.as.push_back(a);
  }

  // update current time
  this->t = t;
  return true;
}


/**
 * Process feature data
 * @param t time
 * @param zcs current measured points
 * @param zcInds current measured point indices
 * @return true on success
 */
bool DynVisIns::ProcessCam(double t, const vector<Vector2d> &zcs, const vector<int> &pntIds) 
{
  assert(useCam);  
  //  assert(this->t >= 0); // assume at least one IMU measurement has arrived 

  //  if (useImu) {
    //    cout << "dt=" << dt << endl;
    if (t - this->t < 0) {
      cout << "[W] DynVisIns::ProcessCam: frame out of sync dt=" << t - this->t << endl;
      return false;
    }
    //  }

  Camera cam; // new camera to be added
  ++camId;    // increase global id (it was initialized to -1)

  // update global time and camera time
  this->t = t;  

  // if not first camera then set delta-t b/n last and this camera
  // and init state to previous cam state
  if (camId > 0) {
    cams[camId - 1].dt = t - tc;
    cam.x = cams[camId - 1].x;
  } else {
    cam.dt = 0;
    cam.x = x0;
  }
   // above one could use the propagated state x instead of x0 to initialize using IMU dead-reconing -- only a good idea if initial pose is correct, otherwise accelerometer-based odometry will be off
  
  this->tc = t;  // update last camera time to to current time
  
  // add observations
  for (int i = 0; i < zcs.size(); ++i) {        
    
    int pntId = pntIds[i];
    const Vector2d &z = zcs[i];
    
    // if this is a new point, then add it to point map
    if (pnts.find(pntId) == pnts.end()) {
      
      Point pnt;
      pnt.usePrior = false;
      pnt.active = false;
      //      pnt.l = Vector3d(1,0,0);
      // generate a spherical measurement
      Vector3d lu((zcs[i][0] - K(0,2))/K(0,0),
                  (zcs[i][1] - K(1,2))/K(1,1),
                  1);    
      lu = lu/lu.norm(); // unit normal in camera frame
      
      pnt.l = cam.x.p + cam.x.R*Ric*lu; // unit normal in spatial frame
      
      pnts[pntId] = pnt;
    } 

    pnts[pntId].zs[camId] = z;
    
    /* covariance of points
       static Vector3d e3(0, 0, 1);
       Vector3d b = e3.cross(lc);
       b = b/b.norm();
       
       Matrix3d Rl;
       Rl.col(0) = lc;
       Rl.col(1) = b;
       Rl.col(2) = lc.cross(b);
       
       Matrix3d Pl;
       Pl << 25, 0, 0, 0, .1, 0, 0, 0, .1;
       
       P.block(15 + i*3, 15 + i*3, 3, 3) = Rl*Pl*Rl.transpose();
    */
  
    cam.pntIds.push_back(pntId);
  }

  // add camera to map
  cams[camId] = cam;

  // This is expensive if we add multiple cameras between calls to compute.  Instead, 
  // only remove cameras once compute is called and reset the prior only once.  Also, set a prior
  // on features from the removed cameras that are still present in the optimization.
  // Foreach camera removed
  //   Foreach feature seen by camera
  //     Remove camera meaurement from feature measurement list
  //     If last measurement
  //       Remove feature from set
  //     Else
  //       Add feature to set

  // check if within current window
  /*
  if (maxCams > 0 && cams.size() > maxCams) {
    assert(cams.size() == maxCams + 1);

    // if prior is used then reset it to the second oldest cam
    if (usePrior)
      ResetPrior(camId0 + 1);

    RemoveCamera(camId0);
    camId0++;
    cout << "[I] DynVisIns::ProcessCam: maxCams=" << maxCams << " window reached. Removing first cam." << endl;
  }
  */

  return true;
}

/**
 * Process stereo feature data
 * @param t time
 * @param zcs current measured points
 * @param zcInds current measured point indices
 * @return true on success
 */
bool DynVisIns::ProcessStereoCam(double t, const vector<Vector3d> &zcs, const vector<int> &pntIds) 
{
  assert(useCam);  

  if (t - this->t < 0) {
    cout << "[W] DynVisIns::ProcessStereoCam: frame out of sync dt=" << t - this->t << endl;
    return false;
  }

  Camera cam; // new camera to be added
  ++camId;    // increase global id (it was initialized to -1)

  // update global time and camera time
  this->t = t;  

  // if not first camera then set delta-t b/n last and this camera
  // and init state to previous cam state
  if (camId > 0) {
    cams[camId - 1].dt = t - tc;
    cam.x = cams[camId - 1].x;
  } else {
    cam.dt = 0;
    cam.x = x0;
  }
   // above one could use the propagated state x instead of x0 to initialize using IMU dead-reconing -- only a good idea if initial pose is correct, otherwise accelerometer-based odometry will be off
  
  this->tc = t;  // update last camera time to to current time
  
  // add observations
  for (int i = 0; i < zcs.size(); ++i) {        
    
    int pntId = pntIds[i];
    const Vector3d &z = zcs[i];
    
    // if this is a new point, then add it to point map
    if (pnts.find(pntId) == pnts.end()) {
      
      Point pnt;
      pnt.usePrior = false;
      pnt.active = false;
      
      pnt.l = cam.x.p + cam.x.R*Ric*z; // unit normal in spatial frame
      
      pnts[pntId] = pnt;
    } 

    pnts[pntId].z3ds[camId] = z;
    cam.pntIds.push_back(pntId);
  }

  // add camera to map
  cams[camId] = cam;

  return true;
}

// pnt_zs_removed keeps a set of all the points which had measurements removed yet are still
//   present in the optimization.
bool DynVisIns::RemoveCamera(int id, std::set<int>* pnt_zs_removed) 
{
  map<int, Camera>::iterator camIter = cams.find(id);
  if (camIter == cams.end()) {
    cout << "[W] DynVisIns::RemoveCamera: cannot find camera with id#" << id << endl;
    return false;
  }
  Camera &cam = camIter->second;

  // go through all points seen by the camera
  vector<int>::iterator iter;
  for (iter = cam.pntIds.begin(); iter != cam.pntIds.end(); ++iter) {
    int pntId = *iter;
    Point &pnt = pnts[pntId];
    pnt.zs.erase(id);  // erase the feature measurement of this point by cam
    if(pnt_zs_removed)
      pnt_zs_removed->insert(pntId);
    if (!pnt.zs.size()) { // if this point now has no measurements then delete it
      pnts.erase(pntId); 
      if(pnt_zs_removed)
        pnt_zs_removed->erase(pntId);
      cout << "[I] DynVisIns::RemoveCam: pnt id#" << pntId << " is invisible and removed." << endl;
    }
  }
  // remove the camera
  cams.erase(id);
}

bool DynVisIns::RemovePoint(int id) 
{
  // TODO
  // currently points are being removed dynamically only within ProcessCam if a camera is being removed
}

bool DynVisIns::ResetPrior(int id, std::set<int>* pnt_ids)
{
  if (!v) {
    cout << "[W] DynVisIns::ResetPrior: optimization vector not set!" << endl;
    return false;
  }

  if (!ceresActive) {
    cout << "[W] DynVisIns::ResetPrior: ceres is not active!" << endl;
    return false;
  }
  
  map<int, Camera>::iterator camIter = cams.find(id);
  if (camIter == cams.end()) {
    cout << "[W] DynVisIns::ResetPrior: camera id#" << id << " not found." << endl;
    return false;
  }
  Camera &cam = camIter->second;
  x0 = cam.x;
  
  ceres::Covariance::Options options;
  // Options below are way too slow...but they take care of rank-deficient Jacobian.
  //options.algorithm_type = ceres::DENSE_SVD;
  //options.null_space_rank = -1;
  ceres::Covariance covariance(options);
  
  vector<pair<const double*, const double*> > covariance_blocks;

  int ind = (id - camId0)*12;    // index into the CERES opt vector

  //  double *r = this->v + ind;
  //  double *p = this->v + ind + 3;
  //  double *dr = this->v + ind + 6;
  //  double *dp = this->v + ind + 9;
  
  // go through pairs of (r,p,dr,dp)
  assert(useDyn || useImu);
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j) {
      double *v1 = this->v + ind + 3*i;
      double *v2 = this->v + ind + 3*j;
      covariance_blocks.push_back(make_pair(v1, v2));
    }
  }

  std::vector<int> v_pt_idxs;
  if(pnt_ids)
  {
    set<int>::iterator pntIter;
    for(pntIter = pnt_ids->begin(); pntIter != pnt_ids->end(); pntIter++)
    {
      // TODO: could have an inverse map available to avoid this search...
      int v_pt_idx;
      for(v_pt_idx = 0; v_pt_idx < l_opti_map->size(); ++v_pt_idx)
      {
        if(l_opti_map->at(v_pt_idx) == *pntIter)
          break;
      }
      if(v_pt_idx == l_opti_map->size())
      {
        cout << "[W] DynVisIns::ResetPrior: failed to find ptId " << *pntIter << " in previous optimization vector (this could be because it was never set to active)...skipping." << endl;
        continue;
      }
      double *vpt = this->v + 12*num_opti_cams + v_pt_idx*3;
      covariance_blocks.push_back(make_pair(vpt,vpt));
      v_pt_idxs.push_back(v_pt_idx);
    }
  }  
  
  assert(problem);
  //CHECK(covariance.Compute(covariance_blocks, problem));
  if(!covariance.Compute(covariance_blocks, problem))
  {
    cout << "[W] DynVisIns::ResetPrior: jacobian is rank deficcient." << endl;
  } 
 
  // CERES works in row major, while Eigen is column major
  Matrix<double, 3, 3, RowMajor> M;
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j) {
      double *v1 = this->v + ind + 3*i;
      double *v2 = this->v + ind + 3*j;
      covariance.GetCovarianceBlock(v1, v2, M.data());
      x0.P.block<3,3>(3*i, 3*j) = M;
    }
  }  
  
  if(pnt_ids)
  {
    Matrix<double, 3, 3, RowMajor> M;
    for (int i = 0; i < v_pt_idxs.size(); ++i) {
      double *vpt = this->v + 12*num_opti_cams + v_pt_idxs[i]*3;
      covariance.GetCovarianceBlock(vpt, vpt, M.data());
      //cout << "[I] DynVisIns::ResetPrior: set pt " << l_opti_map->at(v_pt_idxs[i]) 
      //  << " prior to " << std::endl << M << std::endl;
      pnts[l_opti_map->at(v_pt_idxs[i])].P = M;
      pnts[l_opti_map->at(v_pt_idxs[i])].usePrior = true;
    }  
  }
}


bool DynVisIns::Compute() {

  if(useAnalyticJacs && !useCay)
  {
    cout << "[E] DynVisIns::Compute: must use cayley map if using analytic jacobians.  Jacobians are taken with respect to cayley map." << endl;
    assert(useCay);
  }

  // check if cams within current window
  // reset feature and state priors if removing cams
  if  (maxCams > 0 && v && cams.size() > maxCams) {
    //assert(cams.size() == maxCams + 1);
    std::set<int> pnt_zs_removed;
    int newCamId0 = cams.size() - maxCams + camId0;
    cout << "[I] DynVisIns::Compute: maxCams=" << maxCams << " window reached. Removing " 
      << cams.size() - maxCams << " cams. New cam id0=" << newCamId0 << endl;
    for(int i = camId0; i < newCamId0; ++i)
    {
      if(useFeatPrior)
      {
        RemoveCamera(i, &pnt_zs_removed);
      }
      else
      {
        RemoveCamera(i);
      }
    }

    // if prior is used then reset it to cam at the end of the new window
    if (usePrior)
    {
      if(useFeatPrior)
      {
        ResetPrior(newCamId0, &pnt_zs_removed);
      }
      else
      {
        ResetPrior(newCamId0);
      }
    }
    camId0 = newCamId0;
  }
 
      
  n_good_pnts = 0;
  vector<pair<int,Point&>> gpoints;
  map<int, Point>::iterator gpntIter;
  for (gpntIter = pnts.begin(); gpntIter != pnts.end(); ++gpntIter) {
    // only consider points with >1 observations
    if(gpntIter->second.zs.size() > 1 && (!checkPtActiveFlag || gpntIter->second.active))
    {
      gpoints.push_back(pair<int,Point&>(gpntIter->first, gpntIter->second));
    }
  }
  n_good_pnts = gpoints.size();

  if(n_good_pnts < minPnts)
  {
    cout << "[E] DynVisIns::Compute: fewer than " << minPnts 
      << " good points...aborting optimization." << endl;
    return false;
  }

  if(problem)
    delete problem;

  problem = new ceres::Problem();

  if (v)
    delete[] v;
  if (l_opti_map)
    delete l_opti_map;

  v = new double[12*cams.size() + 3*n_good_pnts + (optBias ? 6 : 0)];
  l_opti_map = new std::vector<int>(n_good_pnts);  
  num_opti_cams = cams.size();

  ToVec(v, l_opti_map);
  if (useCam) {
    // for efficiency, instead of computing a projected covariance on the tangent 
    // of the unit sphere, that needs to be recomputed for every measurement
    // since the projection depends on the measurement,
    // we assume a constant ball of radius sphStd, averaged on the u-v plane
    sphStd = pxStd/sqrt(fx*fx + fy*fy)/2;
    
    assert(this->pnts.size());
    
    
    // first iterate through points 
    //map<int, Point>::iterator pntIter;
    //int i = 0;
    //for (pntIter = pnts.begin(); pntIter != pnts.end(); ++pntIter) {
    for (int i = 0; i < gpoints.size(); i++) {
      int pntId = gpoints[i].first;
      Point &pnt = gpoints[i].second;
      if(pnt.zs.size() <= 1)
      {
        cout << "[I] DynVisIns::Compute: pntId=" << pntId << " only has one observation...skipping." << endl;
        continue;
      }
      if(checkPtActiveFlag && !pnt.active)
      {
        cout << "[I] DynVisIns::Compute: pntId=" << pntId << " not active...skipping." << endl;
        continue;
      }
      double *l = v + 12*cams.size() + 3*i;

      if(useFeatPrior && pnt.usePrior)
      {
        ceres::CostFunction* cost = FeaturePrior::Create(*this, pnt);
        problem->AddResidualBlock(cost,
                                 NULL,
                                 l);
      }

      //std::cout << "pnt #" << pntId << " obs=" << pnt.zs.size() << std::endl;
      map<int, Vector2d>::iterator zIter;
      for (zIter = pnt.zs.begin(); zIter != pnt.zs.end(); ++zIter) {//int i = 0; i < pnt.zs.size(); ++i) {
        int camId = zIter->first;
        Vector2d &z = zIter->second;
        ceres::CostFunction* cost_function;
        if(useAnalyticJacs)
        {
          cost_function = new AnalyticPerspError(*this, z);
        }
        else
        {
          cost_function =
            //        sphMeas ?
            //        SphError::Create(*this, lus[i]) :
            PerspError::Create(*this, z);
        }
        
        double *x = v + 12*(camId - camId0);
        
        //        assert(zCamInds[i] < xs.size());
        //        assert(zInds[i] < ls.size());
        
        
        //      ceres::LossFunction *loss_function = new ceres::HuberLoss(100.0);
        //      ceres::LossFunctionWrapper* loss_function(new ceres::HuberLoss(1.0), ceres::TAKE_OWNERSHIP);
        if(useHuberLoss)
        {  
          problem->AddResidualBlock(cost_function,
                                   new ceres::HuberLoss(2.0),
                                   x, x + 3, l);
        }
        else
        {
          problem->AddResidualBlock(cost_function,
                                   NULL,
                                   x, x + 3, l);
        }

        // for now restrict point coordinates to [-20,20] meters, assuming we're in a small room
        problem->SetParameterLowerBound(l, 0, -20);
        problem->SetParameterLowerBound(l, 1, -20);
        problem->SetParameterLowerBound(l, 2, -20);
        problem->SetParameterUpperBound(l, 0, 20);
        problem->SetParameterUpperBound(l, 1, 20);
        problem->SetParameterUpperBound(l, 2, 20);      
      }

      map<int, Vector3d>::iterator z3dIter;
      for (z3dIter = pnt.z3ds.begin(); z3dIter != pnt.z3ds.end(); ++z3dIter) 
      {
        int camId = z3dIter->first;
        Vector3d &z = z3dIter->second;
        ceres::CostFunction* cost_function;
        //if(useAnalyticJacs)
        //{
        //  cost_function = new AnalyticPerspError(*this, z);
        //}
        //else
        //{
          cost_function =
            //        sphMeas ?
            //        SphError::Create(*this, lus[i]) :
            StereoError::Create(*this, z);
        //}
        
        double *x = v + 12*(camId - camId0);
        
        //        assert(zCamInds[i] < xs.size());
        //        assert(zInds[i] < ls.size());
        
        
        //      ceres::LossFunction *loss_function = new ceres::HuberLoss(100.0);
        //      ceres::LossFunctionWrapper* loss_function(new ceres::HuberLoss(1.0), ceres::TAKE_OWNERSHIP);
        if(useHuberLoss)
        {  
          problem->AddResidualBlock(cost_function,
                                   new ceres::HuberLoss(2.0),
                                   x, x + 3, l);
        }
        else
        {
          problem->AddResidualBlock(cost_function,
                                   NULL,
                                   x, x + 3, l);
        }

        // for now restrict point coordinates to [-20,20] meters, assuming we're in a small room
        problem->SetParameterLowerBound(l, 0, -20);
        problem->SetParameterLowerBound(l, 1, -20);
        problem->SetParameterLowerBound(l, 2, -20);
        problem->SetParameterUpperBound(l, 0, 20);
        problem->SetParameterUpperBound(l, 1, 20);
        problem->SetParameterUpperBound(l, 2, 20);      
      }
      //++i;
    }
  }

  if (useImu) {
    //    assert(xs.size() >= tss.size());
    //    assert(dts.size() == tss.size());
    
    map<int, Camera>::iterator camIter;
    for (camIter = cams.begin(); camIter != cams.end(); ++camIter) {

      int camId = camIter->first;
      Camera &cam = camIter->second;

      // ignore last frame
      if (camId == this->camId)
        continue;
      
      assert(cam.ts.size());
      assert(cam.dt > 0);
      //    for (int i = 0; i < tss.size(); ++i) {
      //      vector<double> &ts = tss[i];
      //      vector<Vector3d> &ws = wss[i];
      //      vector<Vector3d> &as = ass[i];
      //      assert(ts.size() == ws.size());
      //      assert(ts.size() == as.size());
      
      ceres::CostFunction* gyroCost;
      if(useAnalyticJacs)
      {
        gyroCost = new AnalyticGyroCubicError(*this, cam.dt, cam.ts, cam.ws);
      }
      else
      {
        gyroCost = GyroCubicError::Create(*this, cam.dt, cam.ts, cam.ws);
      }
      
      double *xa = v + 12*(camId - camId0);
      double *xb = v + 12*(camId + 1 - camId0);
      
      problem->AddResidualBlock(gyroCost,
                               NULL /* squared loss */,
                               xa, xa + 6, xb, xb + 6);
      
      ceres::CostFunction* accCost;
      if(useAnalyticJacs)
      {
        accCost = new AnalyticAccCubicError(*this, cam.dt, cam.ts, cam.as);
      }
      else
      {
        accCost = AccCubicError::Create(*this, cam.dt, cam.ts, cam.as);
      }
      
      problem->AddResidualBlock(accCost,
                               NULL /* squared loss */,
                               xa, xa + 3, xa + 6, xa + 9,
                               xb, xb + 3, xb + 6, xb + 9);        
      // could also add some box constraints on the state?
    }
  }
  
  if (useDyn) {
    
    //    assert(xs.size() == dts.size() + 1);
    map<int, Camera>::iterator camIter;
    for (camIter = cams.begin(); camIter != cams.end(); ++camIter) {

      
      int camId = camIter->first;
      Camera &cam = camIter->second;

      // ignore last frame
      if (camId == this->camId)
        continue;
      
      assert(cam.dt > 0);

      //    for (int i = 0; i < tss.size(); ++i) {
      //      vector<double> &ts = tss[i];
      //      vector<Vector3d> &ws = wss[i];
      //      vector<Vector3d> &as = ass[i];
      //      assert(ts.size() == ws.size());
      //      assert(ts.size() == as.size());
    
        //    for (int i = 0; i < dts.size(); ++i) {
        //      assert(dts[i] > 0);
      
      double *xa = v + 12*(camId - camId0);
      double *xb = v + 12*(camId + 1 - camId0);
      
      ceres::CostFunction* rotCost = CvCubicRotError::Create(*this, cam.dt);    
      problem->AddResidualBlock(rotCost,
                               NULL /* squared loss */,
                               xa, xa + 6, xb, xb + 6);
      
      ceres::CostFunction* posCost = CvCubicPosError::Create(*this, cam.dt);    
      problem->AddResidualBlock(posCost,
                               NULL /* squared loss */,
                               xa + 3, xa + 9, xb + 3, xb + 9);
      
      // could also add some box constraints on the state?
    }
  }

  if (usePrior) {
    ceres::CostFunction* cost = StatePrior::Create(*this, x0);
    problem->AddResidualBlock(cost,
                             NULL,
                             v, v + 3, v + 6, v + 9);
  }
    
  ceres::Solver::Options options;
  // options.linear_solver_type = ceres::DENSE_SCHUR;
  options.linear_solver_type = ceres::SPARSE_NORMAL_CHOLESKY;
  options.minimizer_progress_to_stdout = true;
  options.max_num_iterations = maxIterations;

  ceres::Solver::Summary summary;
  ceres::Solve(options, problem, &summary);
  std::cout << summary.FullReport() << "\n";  
  
  FromVec(v, l_opti_map);
  ceresActive = true;

  //  delete[] v;
  //  v = 0;

  return true;
  
}

/*
void hermite3(Vector3d &c1, Vector3d &c2, Vector3d &c3,
              const Vector3d &pa, const Vector3d *va,
              const Vector3d &pb, const Vector3d *vb)
{
  Vector3d d = pb - pa;
  c1 = d*va.norm()/d.norm();
  c2 = 3*(pb - pa) - c1 - 2*va;
  c3 = -2*(pb - pa) + c1 + va;
}
*/

// generate synthetic data
  bool DynVisIns::GenData(DynVisIns &tvi, int ns, int np, int ni) 
{
  fx = 453.23520207;
  fy = 453.72298392;
  cx = 391.85891497;
  cy = 282.24403976;
  K(0,0) = fx; K(0,1) = 0; K(0,2) = cx;
  K(1,0) = 0; K(1,1) = fy; K(1,2) = cy;

  // tvi is the "true" VI estimator
  map<int, Camera> &cams = tvi.cams;
  map<int, Point> &pnts = tvi.pnts;

  //  this->ls.resize(ls.size());

  // generate a grid of features on a plane 3 meters ahead
  //  int n1 = sqrt(ls.size());
  //  assert(ls.size() == n1*n1);
  int n1 = sqrt(np);

  for (int i = 0; i < n1; ++i) {
    for (int j = 0; j < n1; ++j) {
      int id = i*n1 + j;        
      // true points are at a vertical plane at distance 3 meters forward
      Point pnt;
      pnt.l = Vector3d( 3, ((double)(j-n1/2.0))/n1, ((double)(i-n1/2.0))/n1);

      pnts[id] = pnt;
      
      pnt.l /= pnt.l.norm();

      // initialize estimator points using unit vectors pointing towards the points
      this->pnts[id] = pnt;
    }
  }    
  cout << "Generated " << pnts.size() << " points" << endl;

  // #of segments
  //  int ns = xs.size() - 1;

  double t = 0;
  double dt = 1.0/ns;  // time-step of each segment
  

  //  this->xs.resize(xs.size());
  //  this->tss.resize(xs.size() - 1);
  //  this->wss.resize(xs.size() - 1);
  //  this->ass.resize(xs.size() - 1);  

  // initialize first state using the prior
  Camera cam0;
  cam0.x = this->x0;  
  cam0.dt = dt;
  tvi.camId++;   // is -1 initially
  cams[tvi.camId] = cam0;
  this->cams[tvi.camId] = cam0;

  //  tvi.xs.push_back(this->x0);

  // initialize simulated path using x0
  Vector3d r, p, dr, v;
  FromState(r, p, dr, v, tvi.x0);
  
  for (int i = 0; i < ns; ++i) {

    Camera cam;
    cam.dt = dt;
            
    // random angular acceleration
    Vector3d aw = dwStd*Vector3d(randn(), randn(), randn());
    Vector3d av = dvStd*Vector3d(randn(), randn(), randn());
    Vector3d jerk(0,0,0);  // zero jerk
    
    Cubic cw(r, dr, dt); cw.b = aw; cw.c = jerk;
    Cubic cv(p, v, dt); cv.b = av; cv.c = jerk;
    
    // generate IMU measurements
    if (useImu) {                        
      Matrix3d D;
      
      for (int j=1; j <=ni; j++) {
        double ti = j*dt/(ni+1); // relative IMU time
        Body3dState xt;
        Vector3d rt, pt, drt, vt, at;
        cw.GetPos(rt, ti);
        cw.GetVel(drt, ti);           
        cv.GetPos(pt, ti); 
        cv.GetVel(vt, ti);
        cv.GetAcc(at, ti);
        
        ToState(xt, rt, pt, drt, vt);
        cam.ts.push_back(ti);   // add an IMU measurement          
        cam.ws.push_back(xt.w);   // gyro reading in body frame
        cam.as.push_back(xt.R.transpose()*(at + g0));  // accel reading in body frame
      }
      
      //      tss[i] = ts;
      //      wss[i] = ws;
      //      ass[i] = as;
    }
    
    // update pos and vel
    cw.GetPos(r, dt);
    cw.GetVel(dr, dt);
    cv.GetPos(p, dt);
    cv.GetVel(v, dt);      


    // true state
    ToState(cam.x, r, p, dr, v);

    tvi.camId++;
    cams[tvi.camId] = cam;
    // init estimator data to first state
    this->camId = tvi.camId;
    this->cams[this->camId] = cam0;
  }
  
  // generate feature meas
  if (useCam)  {
    map<int, Camera>::iterator camIter;

    for (camIter = cams.begin(); camIter != cams.end(); ++camIter){
      int camId = camIter->first;
      Camera &cam = camIter->second;
      Body3dState &x = cam.x;
      map<int, Point>::iterator pntIter;
      for (pntIter = pnts.begin(); pntIter != pnts.end(); ++pntIter) {
        int pntId = pntIter->first;
        Point &pnt = pntIter->second;
        Matrix3d R = x.R*Ric;  // camera rotation
        Vector3d r = R.transpose()*(pnt.l - x.p);
        
        this->pnts[pntId].zs[camId] = K*(r/r[2]);

        cam.pntIds.push_back(pntId);
        // spherical measurements
        //        this->lus.push_back(r/r.norm());      
        
          // pixel measurements
        // this->zs.push_back();
        
        // this->zInds.push_back(j);
        // this->zCamInds.push_back(k);
      }
    }
  }

  /*
  map<int, Camera>::iterator camIter;
  for (camIter = this->cams.begin(); camIter != this->cams.end(); ++camIter){
    int camId = camIter->first;
    Camera &cam = camIter->second;
    Body3dState &x = cam.x;
    cout << "cam id=" << camId << " at pos p=" << cam.x.p.transpose() << endl;
    map<int, Point>::iterator pntIter;
    for (pntIter = this->pnts.begin(); pntIter != this->pnts.end(); ++pntIter) {
      int pntId = pntIter->first;
      Point &pnt = pntIter->second;
      cout << "pnt id=" << pntId << " pos=" << pnt.l.transpose() << " ";
      map<int, Vector2d>::iterator zIter;
      for (zIter = pnt.zs.begin(); zIter != pnt.zs.end(); ++zIter)
        cout << " zid= " << zIter->first << " z=" << zIter->second.transpose() << " ";
    }
    cout << endl;
  }
  */
  
  // useImu = false;
  
  //  cout << "Generated " << pnts.size() << " feature measurements" << endl;
  //  cout << "Generated " << tss.size() << " IMU segments" << endl;

  return true;
}


bool DynVisIns::MakeFeatures(vector<Vector2d> &zs, 
                             vector<int> &pntIds,
                             const Body3dState &x, 
                             const map<int, Point> &pnts)
{
  int l = 780;
  int h = 560;
  zs.clear();
  pntIds.clear();
  map<int, Point>::const_iterator pntIter;
  for (pntIter = pnts.begin(); pntIter != pnts.end(); ++pntIter) {
    int pntId = pntIter->first;
    const Point &pnt = pntIter->second;
    Vector2d z;
    MakeFeature(z, x, pnt.l);

    // Check if feature is in camera window
    if(z(0) >=0 && z(0) <= l && z(1) >= 0 && z(1) <= h)
    {
      zs.push_back(z);
      pntIds.push_back(pntId);
    }
  }
  return true;
}


bool DynVisIns::MakeFeature(Vector2d &z, const Body3dState &x, const Vector3d &l)
{
  Matrix3d R = x.R*Ric;  // camera rotation
  Vector3d r = R.transpose()*(l - x.p);
  z = K*(r/r[2]);
  return true;
}


bool DynVisIns::SimData(DynVisIns &tvi, int ns, int np, int ni, double dt) 
{
  // set camera data
  fx = 453.23520207;
  fy = 453.72298392;
  cx = 391.85891497;
  cy = 282.24403976;
  K(0,0) = fx; K(0,1) = 0; K(0,2) = cx;
  K(1,0) = 0; K(1,1) = fy; K(1,2) = cy;

  // tvi is the "true" VI estimator
  map<int, Camera> &cams = tvi.cams;
  map<int, Point> &pnts = tvi.pnts;

  // number of features on each side of grid
  int n1 = sqrt(np);

  for (int i = 0; i < n1; ++i) {
    for (int j = 0; j < n1; ++j) {
      int id = i*n1 + j;        
      // true points are at a vertical plane at distance 3 meters forward
      Point pnt;
      pnt.l = Vector3d( 3, 9*((double)(j-n1/2.0))/n1, 9*((double)(i-n1/2.0))/n1);
      pnts[id] = pnt;
    }
  }    
  cout << "Generated " << pnts.size() << " points" << endl;

  //  double dt = tf/ns;  // time-step of each segment

  // initialize true first state using the true prior
  Camera cam;
  cam.dt = dt;

  vector<Vector2d> zs;
  vector<int> pntIds;
  if(tvi.camId == -1)
  {
    tvi.t = 0;
    cam.x = tvi.x0;
    tvi.camId++;   // is -1 initially
    cams[tvi.camId] = cam;

    MakeFeatures(zs, pntIds, cam.x, pnts);
    ProcessCam(tvi.t, zs, pntIds);
    for(int i = 0; i < pntIds.size(); i++)
    {
      std::cout << "ID["<<i<<"]="<<pntIds[i] << std::endl;
    }
    std::cout <<"init sim" << std::endl;

  }
  else
  {
    cam = cams[tvi.camId];
  }
  double t = tvi.t;
  std::cout<<"t="<<t<<std::endl;
  //  this->cams[tvi.camId] = cam0;



  // initialize simulated path using x0
  Vector3d r, p, dr, v;
  FromState(r, p, dr, v, cam.x);
    
  for (int i = 0; i < ns; ++i) {
    // random angular acceleration
    Vector3d aw = dwStd*Vector3d(randn(), randn(), randn());
    Vector3d av = dvStd*Vector3d(randn(), randn(), randn());
    Vector3d jerk(0,0,0);  // zero jerk
    
    Cubic cw(r, dr, dt); cw.b = aw; cw.c = jerk;
    Cubic cv(p, v, dt); cv.b = av; cv.c = jerk;
    
    // local IMU time-step
    double dti = dt/(ni+1);

    // generate IMU measurements
    if (useImu) {                        
      Matrix3d D;
      
      for (int j = 1; j <= ni; j++) {
        double ti = j*dti; // relative IMU time
        Body3dState xt;
        Vector3d rt, pt, drt, vt, at;
        cw.GetPos(rt, ti);
        cw.GetVel(drt, ti); 
        cv.GetPos(pt, ti);
        cv.GetVel(vt, ti);
        cv.GetAcc(at, ti);
        
        ToState(xt, rt, pt, drt, vt);

        //ProcessImu(tvi.t + ti, xt.w, xt.R.transpose()*(at + g0));
        ProcessImu(t + ti, xt.w, xt.R.transpose()*(at + g0));
      }
    }
    // update cam time, pos and vel
    //tvi.t += dt;
    t += dt;
    cw.GetPos(r, dt);
    cw.GetVel(dr, dt);
    cv.GetPos(p, dt);
    cv.GetVel(v, dt);      
    
    // true state
    Camera cam;
    cam.dt = dt;            
    ToState(cam.x, r, p, dr, v);
    tvi.camId++;
    cams[tvi.camId] = cam;

    // process cam
    MakeFeatures(zs, pntIds, cam.x, pnts);
    for(int j = 0; j < pntIds.size(); j++)
    {
      std::cout << "ID["<<j<<"]="<<pntIds[j] << std::endl;
    }
    //ProcessCam(tvi.t, zs, pntIds);
    ProcessCam(t, zs, pntIds);
  }

  tvi.t = t;
/*
  map<int, Camera>::iterator camIter;
  for (camIter = this->cams.begin(); camIter != this->cams.end(); ++camIter){
    int camId = camIter->first;
    Camera &cam = camIter->second;
    Body3dState &x = cam.x;
    cout << "cam id=" << camId << " at pos p=" << cam.x.p.transpose() << endl;
    map<int, Point>::iterator pntIter;
    for (pntIter = this->pnts.begin(); pntIter != this->pnts.end(); ++pntIter) {
      int pntId = pntIter->first;
      Point &pnt = pntIter->second;
      cout << "pnt id=" << pntId << " pos=" << pnt.l.transpose() << " ";
      map<int, Vector2d>::iterator zIter;
      for (zIter = pnt.zs.begin(); zIter != pnt.zs.end(); ++zIter)
        cout << " zid= " << zIter->first << " z=" << zIter->second.transpose() << " ";
    }
    cout << endl;
  }
  */
  
  // useImu = false;
  
  //  cout << "Generated " << pnts.size() << " feature measurements" << endl;
  //  cout << "Generated " << tss.size() << " IMU segments" << endl;

  return true;
}



bool DynVisIns::LoadFile(const char* filename) {
  ifstream file(filename);    
  
  if (!file.is_open()) {
    std::cerr << "[E] unable to open file " << filename << std::endl;
    return false;
  }
  
  double dummy;
  double t; // time
  int n;    // nof observations at time t
  int k = 0;
  int msgType;

  do {
    file >> msgType >> t;
    if (msgType == 3) {
      file >> fx >> fy >> cx >> cy >> n;
      
      // params are fx, fy, cx, cy
      K(0,0) = fx; K(0,1) = 0; K(0,2) = cx;
      K(1,0) = 0; K(1,1) = fy; K(1,2) = cy;
      
      if (!n) 
        continue;
      
      vector<Vector2d> zcs(n);
      vector<int> zcInds(n);
      
      for (int i = 0; i < n; ++i) {
        file >> zcs[i][0] >> zcs[i][1];
      }
      for (int i = 0; i < n; ++i) {
        file >> zcInds[i];
      }        
      //      if (x.v.norm() > 0.1)
      if (useCam)
        ProcessCam(t, zcs, zcInds);     
    }

    if (msgType == 1) {
      Vector3d a0, w, a;
      file >> a0[0] >> a0[1] >> a0[2] >> w[0] >> w[1] >> w[2] >> a[0] >> a[1] >> a[2];

      //      imu.a0 = a0;
      //      imu.a0 << 0, 0, 9.81;
      g0 << 0, 0, 9.81;
      cout << a.transpose() << endl;

      if (useImu)
        ProcessImu(t, w, a);
    }
    if (msgType == 2) {
      for (int i = 0; i <6; ++i)
        file >> dummy;
    }      
    k++;
    //    if (dtss.size()>55)
    //      break;
  } while(!file.eof());
  
  file.close();    
  cout << "Added " << cams.size() << " frames and " << pnts.size() << " points" << endl;
  
  return true;
}
