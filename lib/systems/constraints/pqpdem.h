#ifndef GCOP_PQPDEM_H
#define GCOP_PQPDEM_H

#include "constraint.h"
#include "dem.h"
#include "PQP/PQP.h"

namespace gcop {
  
  template <typename T = VectorXd,
    int _nx = Dynamic, 
    int _nu = Dynamic, 
    int _np = Dynamic>
    class PqpDem : public Constraint<T, _nx, _nu, _np, 1> {
  public:

  typedef Matrix<double, 1, 1> Vectorgd;
  typedef Matrix<double, 1, _nx> Matrixgxd;
  typedef Matrix<double, 1, _nu> Matrixgud;
  typedef Matrix<double, 1, _np> Matrixgpd;

  typedef Matrix<double, _nx, 1> Vectornd;
  typedef Matrix<double, _nu, 1> Vectorcd;
  typedef Matrix<double, _np, 1> Vectormd;


    PqpDem(const Dem& dem, double cr = 0.1, double sd = 0.0);
    virtual ~PqpDem();
    
    bool operator()(Vectorgd &g,
                    double t, const T &x, const Vectorcd &u,
                    const Vectormd *p = 0, 
                    Matrixgxd *dgdx = 0, Matrixgud *dgdu = 0,
                    Matrixgpd *dgdp = 0);    
    
    const Dem& dem;   ///< digital elevation map

    PQP_Model *pm;
    PQP_REAL pt[3];
    PQP_REAL pR[3][3];

    PQP_Model *dm;
    PQP_REAL dt[3];
    PQP_REAL dR[3][3];

    PQP_DistanceResult dres;
    PQP_ToleranceResult tres;

    double cr;       ///< collision radius
    double sd;       ///< additional safety distance

  };

  
  template <typename T, int _nx, int _nu, int _np> 
    PqpDem<T, _nx, _nu, _np>::PqpDem(const Dem& dem, double cr, double sd) :
    Constraint<T, _nx, _nu, _np, 1>(), dem(dem), cr(cr), sd(sd)
{
  pm = new PQP_Model;
  dm = new PQP_Model;

  PQP_REAL p1[3], p2[3], p3[3], p4[3];

  pm->BeginModel();
  p1[0] = 0; p1[1] = cr; p1[2] = 0;
  p2[0] = 0; p2[1] = -cr; p2[2] = 0;
  p3[0] = cr; p3[1] = 0; p3[2] = 0;
  pm->AddTri(p1, p2, p3, 0);  
  p1[0] = 0; p1[1] = 0; p1[2] = cr;
  p2[0] = 0; p2[1] = 0; p2[2] = -cr;
  p3[0] = cr; p3[1] = 0; p3[2] = 0;
  pm->AddTri(p1, p2, p3, 1);
  pm->EndModel();
  pm->MemUsage(1);

  dm->BeginModel();
  int count  = 0;
  double p00[3], p10[3], p11[3], p01[3];
  for (int i=0; i < dem.ni-1; ++i) {
    for (int j=0; j < dem.nj-1; ++j) {
      dem.Get(p00, i, j);
      dem.Get(p10, i+1, j);
      dem.Get(p11, i+1, j+1);
      dem.Get(p01, i, j+1);
      p1[0] = p00[0]; p1[1] = p00[1]; p1[2] = p00[2];
      p2[0] = p10[0]; p2[1] = p10[1]; p2[2] = p10[2];
      p3[0] = p11[0]; p3[1] = p11[1]; p3[2] = p11[2];
      p4[0] = p01[0]; p4[1] = p01[1]; p4[2] = p01[2];
      dm->AddTri(p1, p2, p3, count);
      dm->AddTri(p1, p3, p4, count+1);
      count += 2;
    }
  }
  dm->EndModel();
  dm->MemUsage(1);

  pt[0] = pt[1] = pt[2] = 0;
  pR[0][0] = pR[1][1] = pR[2][2] = 1.0;
  pR[0][1] = pR[1][0] = pR[2][0] = 0.0;
  pR[0][2] = pR[1][2] = pR[2][1] = 0.0;

  dt[0] = dt[1] = dt[2] = 0;
  dR[0][0] = dR[1][1] = dR[2][2] = 1.0;
  dR[0][1] = dR[1][0] = dR[2][0] = 0.0;
  dR[0][2] = dR[1][2] = dR[2][1] = 0.0;
}

  template <typename T, int _nx, int _nu, int _np> 
    PqpDem<T, _nx, _nu, _np>::~PqpDem()
  {
    delete dm;
    delete pm;
  }
  
  

template <typename T, int _nx, int _nu, int _np> 
    bool PqpDem<T, _nx, _nu, _np>::operator()(Vectorgd &g,
                                              double t, const T &x, const Vectorcd &u,
                                              const Vectormd *rho, 
                                              Matrixgxd *dgdx, Matrixgud *dgdu,
                                              Matrixgpd *dgdp)
{
  //  const double *q = s.x + oi;
  //  SET3(pt, q);
  Vector3d p = x.second.head(3); // position
  pt[0] = p[0];
  pt[1] = p[1];
  pt[2] = p[2];

  int res = PQP_Distance(&dres, pR, pt, pm, dR, dt, dm, 0.0, 0.0);

  if (res != PQP_OK)
    cerr << "Warning: PqpDem:Distance: res=" << res << endl;

  assert(res == PQP_OK);

  double d = dres.Distance();

  // subtract safety distance
  d -= sd;

  //  assert(d>=0);

  //  if (dem.Inside(dres.P1()[0], dres.P1()[1], dres.P1()[2]))

  bool in = dem.Inside(p[0], p[1], p[2]);
  if (in)
    d = -d;
  
  // cout << "d=" << d << endl;

  if (dgdx && 0) {
    Vector3d p2(dres.P2());
    Vector3d dp = p2 - p;
    dp.normalize();
    dgdx->segment(3,3) = dp;
    //    dgdx->head(3).normalize();

    //    MINUS3(g, q, dres.P2());
    //    double gn = NORM3(g);
    //    DIV3(g,g,gn);
    if (in) {
      dgdx->segment(3,3) = -dp;
    }
  }

  g[0] = -d;
  return (d > 0);
}



};



#endif
