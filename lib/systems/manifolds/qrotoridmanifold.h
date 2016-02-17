#ifndef QROTORIDMANIFOLD_H
#define QROTORIDMANIFOLD_H

#include "body3dmanifold.h"

namespace gcop {
  using namespace std;
  using namespace Eigen;

  typedef Matrix<double, 12, 1> Vector12d;
  typedef Matrix<double, 15, 1> Vector15d;
  typedef Matrix<double, 15, 15> Matrix15d;

  /** Quadrotor state for system ID model. It is an extension of Body3dState. It adds the commanded rpy as a virtual state.
    */

  class QRotorIDState : public Body3dState {
  public:
      Vector3d u; ///< Commanded rpy
      // IF needed add the covariance matrix for full state
      void Clear()
      {
          Body3dState::Clear();
          u.setZero();
      }
  };

  class QRotorIDManifold : public Manifold<QRotorIDState, 15> {
  public:
    static QRotorIDManifold& Instance()
    {
        static QRotorIDManifold instance;
        return instance;
    }

    void Lift(Vector15d &v,
              const QRotorIDState &xa,
              const QRotorIDState &xb)
    {
        Vector12d vref;
        body_manifold.Lift(vref, xa, xb);
        v.head<12>() = vref;
        Vector3d v_tail;
        v_tail = xb.u - xa.u;
        for(int i = 0; i < 3; i++)
        {
          v_tail(i) = (v_tail(i) > M_PI)?(v_tail(i) - 2*M_PI):(v_tail(i) < -M_PI)?(v_tail(i)+2*M_PI):v_tail(i);
        }
        v.tail<3>() = v_tail;
    }

    void Retract(QRotorIDState &xb,
                 const QRotorIDState &xa,
                 const Vector15d &v)
    {
        body_manifold.Retract(xb,xa,v.head(12));
        xb.u = xa.u + v.tail<3>();
        for(int i = 0; i < 3; i++)
        {
          xb.u(i) = (xb.u(i) > M_PI)?(xb.u(i) - 2*M_PI):(xb.u(i) < -M_PI)?(xb.u(i)+2*M_PI):xb.u(i);
        }
    }

    void dtau(Matrix15d &M, const Vector15d &v)
    {
        M.setIdentity();
        Matrix12d mtopcorner = M.topLeftCorner<12,12>();
        Vector12d vtop = v.head<12>();
        body_manifold.dtau(mtopcorner,vtop);
        M.topLeftCorner<12,12>() = mtopcorner;
    }

    void dtauinv(Matrix15d &M, const Vector15d &v)
    {
        M.setIdentity();
        Matrix12d mtopcorner = M.topLeftCorner<12,12>();
        Vector12d vtop = v.head<12>();
        body_manifold.dtauinv(mtopcorner, vtop);
        M.topLeftCorner<12,12>() = mtopcorner;
    }

    void Adtau(Matrix15d &M, const Vector15d &v)
    {
        M.setIdentity();
        Matrix12d mtopcorner = M.topLeftCorner<12,12>();
        Vector12d vtop = v.head<12>();
        body_manifold.Adtau(mtopcorner, vtop);
        M.topLeftCorner<12,12>() = mtopcorner;
    }

    bool useCay;   ///< whether to use the Cayley map instead of exponential

  private:
    QRotorIDManifold(): body_manifold(Body3dManifold::Instance()), useCay(false)
    {}
    Body3dManifold &body_manifold;
  };
}

#endif // QROTORIDMANIFOLD_H
