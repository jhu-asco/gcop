#include <limits>
#include <iostream>
#include <utility>
#include "mbs.h"
#include <iomanip>


using namespace gcop;

void print(const MbsState &x)
{
  for (int i =0; i < x.gs.size(); ++i) {
    cout << "gs[" << i << "]=" << endl << x.gs[i] << endl;
  }
  for (int i =0; i < x.vs.size(); ++i) {
    cout << "vs[" << i << "]=" << endl << x.vs[i].transpose() << endl;
  }  
  for (int i =0; i < x.dgs.size(); ++i) {
    cout << "dgs[" << i << "]=" << endl << x.dgs[i] << endl;
  }

  cout << "r=" << x.r.transpose() << endl;
  cout << "dr=" << x.dr.transpose() << endl;
}


Mbs::Mbs(int nb, int c, bool fixed) : nb(nb), fixed(fixed),
                                      System(*new MbsManifold(nb, fixed), *new Rn<>(c)), 
                                      links(nb),
                                      joints(nb-1), 
                                      pis(nb), cs(nb), se3(SE3::Instance()),
                                      method(EULER),
                                      iters(2),
                                      debug(false), basetype(FLOATBASE),
                                      damping(VectorXd::Zero(nb - 1)),
                                      ag(0,0,-9.81)
                                      //                          fq(*this), 
                                      //                          fva(*this),
                                      //                          fvb(*this), 
                                      //                          fu(*this)
{
  //ag << 0, 0, -9.81;
}
 
  
Mbs::~Mbs()
{
  delete &U;
  delete &X;
}

void Mbs::Force(VectorXd &f, double t, const MbsState &x, const VectorXd &u,
                MatrixXd *A, MatrixXd *B) 
{
  if(basetype == FLOATBASE) {
    assert(!fixed);
    assert(6 + nb - 1 == f.size());
    f.head(6) = u.head(6);
  } else if(basetype == AIRBASE) {
    assert(!fixed);
    f[0] = u[0];
    f[1] = u[1];
    f[2] = u[2];
    f[3] = 0;
    f[4] = 0;
    f[5] = u[3];
  }
  
  // this is ok for fixed base also 
  f.tail(nb-1) = u.tail(nb-1) - damping.cwiseProduct(x.dr);
  
  if (A)
    A->setZero();
  
  if (B) {
    if(basetype == AIRBASE) {
      B->setZero();
      (*B)(0,0) = 1;
      (*B)(1,1) = 1;
      (*B)(2,2) = 1;
      for(int count = 5;count < f.size();count++)
        (*B)(count,count-2) = 1;
    } else {
      B->setIdentity();      
    }
  }
}

void Mbs::Init() 
{
  //  Ips[0] = links[0].I.asDiagonal();
  
  for (int i = 0; i < nb-1; ++i) {    
    Joint &jnt = joints[i];
    jnt.Init();    
    //    Ips[i+1] = jnt.Ac.transpose()*links[i+1].I.asDiagonal()*jnt.Ac;
  }
  
  pis[0] = -1; // root has no parent

  // assume that parent indices pis[*] are already filled in   
  for (int i= nb-1; i > 0; --i)
    cs[pis[i]].push_back(i);
}


double Mbs::F(VectorXd &v, double t, const MbsState &x, 
              const VectorXd &u, double h,
              MatrixXd *A, MatrixXd *B)
{
  MatrixXd M(nb + 5, nb + 5);
  Mass(M, x);
  
  //  cout << M << endl;
     
  VectorXd a(nb + 5);  // accelerations
  
  VectorXd b(nb + 5); // bias
  Bias(b, t, x);
  //DBias(b, t, x, x, h);

  VectorXd f(nb + 5);
  Force(f, t, x, u);  // compute control/external forces

  //  cout << "control fu=" << fu << endl;

  if (debug)
    cout << "f=" << f << endl;
  
  //    f.setZero();
  
  // compute acceleration
  LLT<MatrixXd> llt;
  llt.compute(M);
  if (llt.info() == Eigen::Success) {
    a = llt.solve(f - b);
  } else {
    cout << "[W] Mbs::F Mass matrix not positive definite!" << endl;
    return 0;
  }
  
  // a = -M.inverse()*f; 
  v.head(6) = h*(x.vs[0] + h*a.head(6));            // updated base velocity
  v.segment(6, nb-1) = h*(x.dr + h*a.tail(nb-1));   // updated joint velocities
  v.segment(nb+5, 6) = h*a.head(6);                 
  v.segment(nb+11, nb-1) = h*a.tail(nb-1);
}


void Mbs::Mass(MatrixXd &M, const MbsState &x) const 
{
  int n = nb - 1 + 6*(!fixed);

  // assumes that FK has been called on x, i.e. that it is consistent
  M.setZero();

  vector<Matrix6d> Ics(nb);
  
  for (int i=0; i < nb; ++i) {
    Ics[i] = links[i].I.asDiagonal();
  }

  int i0 = 6*(!fixed);
  
  for(int i = nb-1; i >= 0; --i) {
    int j = pis[i];
    if (j >= 0) { // if has a parent
      Matrix6d Aij;
      Matrix4d gi;
      assert(i>0);
      
      se3.inv(gi, x.dgs[i-1]);    
      
      se3.Ad(Aij, gi);
      Ics[j] += Aij.transpose()*Ics[i]*Aij;
    }
    
    if (i > 0) { // if not the root
      Vector6d F(Ics[i]*joints[i-1].S);
      M(i0 + i - 1, i0 + i - 1) = joints[i-1].S.dot(F);
      
      j=i;
      while(pis[j] >= 0) {
        int l = pis[j];
        Matrix6d Ajl;
        Matrix4d gi;
        se3.inv(gi, x.dgs[j-1]);
        se3.Ad(Ajl, gi);
        F = Ajl.transpose()*F;
        j = l;
        if (j > 0) {
          M(i0 + i - 1, i0 + j - 1) = F.dot(joints[j - 1].S);
          M(i0 + j - 1, i0 + i - 1) = M(i0 + i - 1, i0 + j - 1);
        }
      }
      if (!fixed) {
        M.block(0, 6+i-1, 6, 1) = F;
        M.block(6+i-1, 0, 1, 6) = F.transpose();
      }
    }
  }
  if (!fixed)
    M.topLeftCorner(6,6) = Ics[0];
}

  
void Mbs::FK(MbsState &x)
{
  // given: x.gs[0], x.r 
  // compute: x.gs[*], x.dgs[*], x.Ms[*]
  
  Matrix4d dg;
  for (int i = 1; i < nb; ++i) {
    int pi = pis[i];
    se3.exp(dg, x.r[i-1]*joints[i-1].a);
    x.dgs[i-1] = joints[i-1].gp*dg;
    x.gs[i] = x.gs[pi]*x.dgs[i-1]*joints[i-1].gci;
  }
}


void Mbs::ID(VectorXd &f,
             double t, const MbsState &x, const VectorXd &u) 
{
  VectorXd b(nb + 5); // bias
  Bias(b, t, x);
  //  DBias(b, t, x, x, h);

  VectorXd fu(nb + 5);
  Force(fu, t, x, u);  // compute control/external forces

  //  cout << "control fu=" << fu << endl;

  f = b - fu;
}



void Mbs::Bias(VectorXd &b,
               double t,
               const MbsState &x) const
{
  
  // assume that Kstep was called so that all velocities are propagated
  //   VectorXd b(nb+5); // bias
    
  vector<Vector6d> ps(nb);
  Matrix6d A;
  Matrix4d gi;

  Vector6d gr;
  gr.setZero();

  int i0 = 6*(!fixed);

  for (int i = nb - 1; i >= 0; --i) {

    const Vector6d &v = x.vs[i];

    //    cout << "Bias v[ " << i << "]=" << v.transpose() << endl;

    const Vector6d &I = links[i].I;
    Vector6d mu = I.cwiseProduct(v);

    gr.tail<3>() = links[i].m*(x.gs[i].topLeftCorner<3,3>().transpose()*ag);
    
    //    ps[i].head(3) = v.head<3>().cross(mu.head<3>()) + v.tail<3>().cross(mu.tail<3>());
    //    ps[i].tail(3) = v.head<3>().cross(mu.tail<3>());    
    ps[i].head(3) = v.head<3>().cross(mu.head<3>());
    ps[i].tail(3) = v.head<3>().cross(mu.tail<3>());    

    ps[i] = ps[i] - gr;

    if (debug) {
      //cout << "i=" << i << endl;
      //      cout << "v=" << v << endl;
      //      cout << "mua=" << mua << endl;
      //      cout << "mub=" << mub << endl;      
      //      cout << "ps[" << i << "]=" << ps[i] << endl;
    }
    
    for (int ci = 0; ci < cs[i].size(); ++ci) {
      
      int j = cs[i][ci];   // index of ci-th child of i
      
      if (debug)
        cout << "CHILD: i=" << i << " j=" << j << endl;
      assert(j > 0);
      se3.inv(gi, x.dgs[j-1]);
      se3.Ad(A, gi);
      ps[i] += A.transpose()*ps[j];
    }
    
    if (debug)
      cout << "i=" << i << "ps=" << ps[i] << endl;
    
    if (i > 0)
      b[i0 + i - 1] = ps[i].dot(joints[i-1].S);
    else if (!fixed)
      b.head(i0) = ps[0];
  }
}


double Mbs::HeunStep(MbsState& xb, double t, const MbsState& xa,
                     const VectorXd &u, double h,
                     MatrixXd *A, MatrixXd *B)
{
  VectorXd a(nb+5);
  Acc(a, t, xa, u);

  MbsState xn(nb);
  xn.vs[0] = xa.vs[0] + h*a.head(6);
  xn.dr = xa.dr + h*a.tail(nb-1);
  KStep(xn, xa, h);
  VectorXd an(nb+5);
  Acc(an, t+h, xn, u);  

  xb.vs[0] = xa.vs[0] + h/2*(a.head(6) + an.head(6));
  xb.dr = xa.dr + h/2*(a.tail(nb - 1) + an.tail(nb - 1));
  KStep(xb, xa, h);
  return 0;
}


double Mbs::EulerStep(MbsState& xb, double t, const MbsState& xa,
                      const VectorXd &u, double h,
                      MatrixXd *A, MatrixXd *B)
{
  int n = nb - 1 + 6*(!fixed);

  VectorXd a(n);
  Acc(a, t, xa, u);

  if (!fixed)
    xb.vs[0] = xa.vs[0] + h*a.head(6);

  xb.dr = xa.dr + h*a.tail(nb - 1);

  ClampVelocity(xb);

  KStep(xb, xa, h);

  return 0;

}


void Mbs::ClampPose(MbsState &x, int i) const
{
  assert(i >= 0 && i < nb);
  for (int j = 0; j < 3; ++j) {
    if (x.gs[i](j,3) > X.ub.gs[i](j,3))
      x.gs[i](j,3) = X.ub.gs[i](j,3);
    else if (x.gs[i](j,3) < X.lb.gs[i](j,3))
      x.gs[i](j,3) = X.lb.gs[i](j,3);
  }
}


void Mbs::ClampShape(MbsState &x, int i) const
{
  assert(i >= 0 && i < x.r.size());
  if (x.r[i] > X.ub.r[i])
    x.r[i] = X.ub.r[i];
  else if (x.r[i] < X.lb.r[i])
    x.r[i] = X.lb.r[i];
}


void Mbs::ClampVelocity(MbsState &x) const
{
  if (!fixed) {
    Vector6d &v = x.vs[0];
    const Vector6d &vu = X.ub.vs[0];
    const Vector6d &vl = X.lb.vs[0];
    
    for (int i = 0; i < 6; ++i)
      if (v[i] > vu[i])
        v[i] = vu[i];
      else if (v[i] < vl[i])
        v[i] = vl[i];
  }

  for (int i = 0; i < nb-1; ++i)
    if (x.dr[i] > X.ub.dr[i])
      x.dr[i] = X.ub.dr[i];
    else if (x.dr[i] < X.lb.dr[i])
      x.dr[i] = X.lb.dr[i];
}


void Mbs::Acc(VectorXd &a, double t, const MbsState &x, const VectorXd &u)
{
  int n = nb - 1 + 6*(!fixed);

  MatrixXd M(n, n);
  Mass(M, x);
  //  M.setIdentity();

  VectorXd b(n); // bias
  Bias(b, t, x);
  //  b.setZero();
  
  //  cout << "b=" << b.transpose() << endl;

  VectorXd f(n);
  Force(f, t, x, u);  // compute control/external forces
  
  if (debug)
    cout << "f=" << f << endl;  
  
  // compute acceleration
  LLT<MatrixXd> llt;
  llt.compute(M);
  if (llt.info() == Eigen::Success) {
    a = llt.solve(f - b);
  } else {
    cout << "[W] Mbs::Acc Mass matrix not positive definite!" << endl;
  }  
}



double Mbs::Step(MbsState& xb, double t, const MbsState& xa,
                 const VectorXd &u, double h,
                 MatrixXd *A, MatrixXd *B)
{
  //  if (method == EULER)
  //    return System::Step(xb, t, xa, u, h, A, B);

  
  if (method == EULER)
    return EulerStep(xb, t, xa, u, h, A, B);
  else if (method == HEUN)
    return HeunStep(xb, t, xa, u, h, A, B);
  else if (method == TRAP)
    return TrapStep(xb, t, xa, u, h, A, B);
  else 
    cout << "[W] Mbs::Step: unsupported method " << method << endl;
  return 0;
}


double Mbs::TrapStep(MbsState &xb, double t, const MbsState& xa,
                     const VectorXd &u, double h,
                     MatrixXd *A, MatrixXd *B)
{
  // initialize new velocity with old velocity
  VectorXd vdr(nb + 5);
  vdr.head(6) = xa.vs[0];
  vdr.tail(nb-1) = xa.dr;

  // Newton-Euler error
  VectorXd e(nb + 5);

  // Jacobian
  VectorXd em(nb + 5);
  VectorXd ep(nb + 5);
  VectorXd dvdr(nb + 5); 
  MatrixXd De(nb + 5, nb + 5);

  // search direction
  VectorXd d(nb + 5);

  double eps = 1e-3;

  for (int j = 0; j < iters; ++j) {    
    NE(e, vdr, xb, t, xa, u, h);    
    NewtonEulerJacobian(De, xb, xa, h);

    vdr = vdr - De.lu().solve(e);
    continue;

    // finite differences
    for (int i = 0; i < nb + 5; ++i) {
      dvdr.setZero();
      dvdr[i] = eps;
      //      NE(em, vdr - dvdr, xb, t, xa, u, h);
      NE(ep, vdr + dvdr, xb, t, xa, u, h);
      
      //      De.col(i) = (ep - em)/(2*eps);
      De.col(i) = (ep - e)/(eps);
      //      std::cout << std::setprecision(15) ;
      if (0 && i==nb+2) {
        cout << "e=" << e.transpose() << endl;       
        cout << "ep=" << ep.transpose() << endl;
        cout << "em=" << em.transpose() << endl;
        cout << "KUR" << endl;
        print(xa);
        print(xb);
        cout << "PUT" << endl;
        //        exit(0);
      }
      //      assert(De.col(i).norm() > 1e-10);    
    }  

    vdr = vdr - De.lu().solve(e);
    continue;
    
    FullPivHouseholderQR<MatrixXd> qr;
    qr.compute(De);


    
    if(qr.rank() == nb + 5) {
      
      // take a Newton step
      //    vdr = vdr - De.inverse()*e;
      //    De.lu().solve(e, &d);
      //vdr = vdr - d;
      
      //    vdr = vdr - De.lu().solve(e);
      
      /*
        FullPivHouseholderQR<MatrixXd> qr1;
        qr1.compute(De.block(0,0,6,6));
        cout << "J rank=" << qr1.rank() << endl;
        
        FullPivHouseholderQR<MatrixXd> qr2;
        qr2.compute(De.block(6,6,nb-1,nb-1));
        cout << "m rank=" << qr2.rank() << endl;
        
        FullPivHouseholderQR<MatrixXd> qr3;
        qr3.compute(De.block(6,6,3,3));
        cout << "m rank=" << qr3.rank() << endl;
        
        FullPivHouseholderQR<MatrixXd> qr4;
        qr4.compute(De.block(9,9,3,3));
        cout << "m rank=" << qr4.rank() << endl;
      */
      
      //    cout << "rank=" << qr.rank() << endl;
      
      VectorXd a = qr.solve(e);
      vdr = vdr - a;
      //    cout << "De=" << endl << De << endl;
      //    assert(e.isApprox(De*a));    
    }
  }
  
  if (e.norm() > 1e-3) {
    cout << "[W] Mbs::Step: residual e seems high e=" << e << endl;
  }
  
  xb.vs[0] = vdr.head(6);
  xb.dr = vdr.tail(nb-1);
  KStep(xb, xa, h);

  if (B) {
    assert(B->rows() == 2*(nb + 5));
    assert(B->cols() == c);
    MatrixXd Bv(nb + 5, c);
    VectorXd f(nb + 5);
    Force(f, t, xa, u, 0, &Bv);  // compute control/external forces
    B->topRows(nb + 5).setZero();
    B->bottomRows(nb + 5) = -h*Bv;
  }
}



void Mbs::NE(VectorXd &e, const VectorXd &vdr,
             MbsState &xb,
             double t, const MbsState &xa,
             const VectorXd &u, double h)
{
  xb.vs[0] = vdr.head(6);
  xb.dr = vdr.tail(nb-1);
  KStep(xb, xa, h);
  VectorXd bd(nb + 5);

  DBias(bd, t, xb, xa, h);

  VectorXd f(nb + 5);
  Force(f, t, xa, u);  // compute control/external forces

  e = bd - h*f;
}


void Mbs::DBias(VectorXd &b,
                double t,
                const MbsState &xb,
                const MbsState &xa, double h)
{
  
  // assume that Kstep was called so that all velocities are propagated
  //   VectorXd b(nb+5); // bias
    
  vector<Vector6d> ps(nb);
  Matrix6d Da, Db;
  Matrix6d A;
  Matrix4d gi;

  Vector6d gr;
  gr.setZero();
  
  for (int i = nb - 1; i >= 0; --i) {
    const Vector6d &va = xa.vs[i];
    const Vector6d &vb = xb.vs[i];
    const Vector6d &I = links[i].I;
    
    gr.tail<3>() = links[i].m*(xa.gs[i].topLeftCorner<3,3>().transpose()*ag);

    /*    
    se3.dcayinv(Db, h*vb);
    se3.dcayinv(Da, -h*va);
    ps[i] = Db.transpose()*(I.cwiseProduct(vb)) -  Da.transpose()*(I.cwiseProduct(va)) - h*gr;
    */
    Vector6d pb, pa;
    se3.tlnmu(pb, h*vb, I.cwiseProduct(vb));
    se3.tlnmu(pa, -h*va, I.cwiseProduct(va));
    ps[i] = pb - pa - h*gr;
    //    ps[i] = I.cwiseProduct(vb) - I.cwiseProduct(va) - h*gr;

    if (debug) {
      cout << "i=" << i << endl;
      cout << "va=" << va << endl;
      cout << "vb=" << vb << endl;
      //      cout << "mua=" << mua << endl;
      //      cout << "mub=" << mub << endl;      
      cout << "ps[" << i << "]=" << ps[i] << endl;
    }
    
    for (int ci = 0; ci < cs[i].size(); ++ci) {
      
      int j = cs[i][ci];   // index of ci-th child of i
      
      if (debug)
        cout << "CHILD: i=" << i << " j=" << j << endl;
      assert(j > 0);
      se3.inv(gi, xa.dgs[j-1]);
      se3.Ad(A, gi);
      ps[i] += A.transpose()*ps[j];
    }
    
    //      if (i > 0)
    //        ps[i] = joints[i-1].Ac.transpose()*ps[i];
    
    if (debug)
      cout << "i=" << i << "ps=" << ps[i].transpose() << endl;
    
    if (i > 0)
      b[6+i-1] = ps[i].dot(joints[i-1].S);
    else
      b.head(6) = ps[0];
  }
}




void Mbs::Rec(MbsState &x, double h)
{
  // given reduced state: x.gs[0], x.r, x.vs[0], x.dr
  // compute: x.gs[*], x.dgs[*], x.vs[*]

  // kinematics
  FK(x);
  
  assert(x.gs.size() == nb);
  assert(x.vs.size() == nb);

  // previous state
  MbsState xp(x);
  Matrix4d g;

  //  se3.cay(g, -h*x.vs[0]); // get previous pose
  se3.exp(g, -h*x.vs[0]); // get previous pose    
  xp.gs[0] = x.gs[0]*g;

  xp.r = x.r - h*x.dr;    // get previous joint angles
  FK(xp);                 // expand previous configuration
  //std::cout << std::setprecision(5) ;

  
  for (int i = 1; i < nb; ++i) { // take the difference and divide by time for each body
    se3.inv(g, xp.gs[i]);
		//cout<<"xp.gs["<<i<<"]: "<<endl<<xp.gs[i]<<endl;
		//cout<<"g"<<endl<<g<<endl;
		//cout<<"g*x.gs["<<i<<"]: "<<endl<<g*x.gs[i]<<endl;

    //    se3.cayinv(x.vs[i], g*x.gs[i]);
    se3.log(x.vs[i], g*x.gs[i]);

    //    cout<<"x.vs["<<i<<"]: "<<endl<<x.vs[i]<<endl;

    x.vs[i] /= h;
		//cout<<"h "<<h<<endl;
		//cout<<"x.vs["<<i<<"]: "<<endl<<x.vs[i]<<endl;
  }
}



void Mbs::KStep(MbsState &xb, const MbsState &xa, double h, bool impl) {
  
  // given: xa, xb.vs[0], xb.dr[*]
  // compute: xb.dgs[*], xb.gs[*],  xb.vs[*], xb.r[*]    

  Matrix4d dg;
  Matrix4d gi;
    
  if (!fixed) {
    se3.exp(dg, impl ? h*xb.vs[0] : h*xa.vs[0]);
    xb.gs[0] = xa.gs[0]*dg;
    ClampPose(xb, 0);
  }
  
  for (int i = 1; i < nb; ++i) {
    Joint &jnt = joints[i-1];
    xb.r[i-1] = xa.r[i-1] + (impl ? h*xb.dr[i-1] : h*xa.dr[i-1]);
    ClampShape(xb, i-1);

    //    cout << xb.r[i-1]*jnt.a << endl;

    se3.exp(dg, xb.r[i-1]*jnt.a);
    
    assert(!std::isnan(dg(0,0)));
    xb.dgs[i-1] = jnt.gp*dg*jnt.gci;  // this is also in FK
    
    int pi = pis[i];
    xb.gs[i] = xb.gs[pi]*xb.dgs[i-1];

    //      se3.cayinv(xb.vs[i], gi*xb.gs[pi]*xb.dgs[i-1]);
    se3.inv(gi, xa.gs[i]);

    se3.log(xb.vs[i], gi*xb.gs[i]);
    assert(!std::isnan(xb.vs[i][0]));
    
    xb.vs[i] /= h;
  }
}


void Mbs::NewtonEulerJacobian(MatrixXd &De, const MbsState &xb, const MbsState &xa, double h) 
{  
  // given: xa, xb.vs[0], xb.dr[*]
  // compute: xb.dgs[*], xb.gs[*],  xb.vs[*], xb.r[*]    

  Matrix6d Dv[nb];
  MatrixXd Dr[nb];

  Matrix6d Dp; 
  Matrix6d Dm; 
  se3.tln(Dp, h*xb.vs[0]);  // this must be dtau(-h*v[0])

  Dv[0] = Dp;
  Dr[0] = MatrixXd::Zero(6, nb-1);
  
  Matrix4d dg;
  Matrix4d gi;
  Matrix6d A;
    

  for (int i = 1; i < nb; ++i) {
    Joint &jnt = joints[i-1];    
    int pi = pis[i];

    se3.inv(gi, xa.dgs[i-1]);
    se3.Ad(A, gi);

    Dv[i] = A*Dv[pi];
    Dr[i] = A*Dr[pi];
    Dr[i].col(i-1) += jnt.S;
  }

  MatrixXd P[nb];

  for (int i = nb - 1; i >= 0; --i) {
    const Vector6d &v = xb.vs[i];
    const Vector6d &I = links[i].I;
    
    se3.adt(A, -I.cwiseProduct(v)/2);
    se3.tln(Dp, h*v);
    A = A + Dp.transpose()*I.asDiagonal();
    se3.tln(Dm, -h*v);
    A = A*Dm;
    P[i].resize(6, nb + 5);
    P[i].leftCols(6) = A*Dv[i];
    P[i].rightCols(nb-1) = A*Dr[i];

    for (int ci = 0; ci < cs[i].size(); ++ci) {      
      int j = cs[i][ci];   // index of ci-th child of i      
      assert(j > 0);

      se3.inv(gi, xa.dgs[j-1]);
      se3.Ad(A, gi);

      P[i] += A.transpose()*P[j];
    }
    
    if (i > 0) {
      De.row(6 + i - 1) = joints[i-1].S.transpose()*P[i];
    } else {
      De.topRows(6) = P[0];
    }
  }
}
