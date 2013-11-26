#include <limits>
#include <iostream>
#include <utility>
#include "mbs.h"

using namespace gcop;

Mbs::Mbs(int nb, int c) : nb(nb), 
                          System(*new MbsManifold(nb), *new Rn<>(c)), 
                          links(nb),
                          joints(nb-1), Ips(nb),
                          pis(nb), cs(nb), se3(SE3::Instance()),
                          debug(false)
                          //                          fq(*this), 
                          //                          fva(*this),
                          //                          fvb(*this), 
                          //                          fu(*this)
{
}
  
  
Mbs::~Mbs()
{
  delete &U;
  delete &X;
}

void Mbs::Init() 
{
  
  Ips[0] = links[0].I.asDiagonal();
  
  for (int i = 0; i < nb-1; ++i) {
    
    Joint &jnt = joints[i];
    
    se3.Ad(jnt.Ac, jnt.gc);
    se3.inv(jnt.gpi, jnt.gp);
    se3.inv(jnt.gci, jnt.gc);
    
    jnt.S = jnt.Ac*jnt.a;
    
    Ips[i+1] = jnt.Ac.transpose()*links[i+1].I.asDiagonal()*jnt.Ac;
  }
  
  pis[0] = -1;
  
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
  VectorXd f(nb + 5);  // bias forces
  
  // compute bias
  ID(f, t, x, x, u, h);
  
  if (debug)
    cout << "f=" << f << endl;
  
  //    f.setZero();
  
  // compute acceleration
  LLT<MatrixXd> llt;
  llt.compute(M);
  if (llt.info() == Eigen::Success) {
    a = -llt.solve(f);
  } else {
    cout << "[W] Mbs::F Mass matrix not positive definite!" << endl;
    return 0;
  }
  
  // a = -M.inverse()*f;
  v.head(6) = h*(x.vs[0] + h*a.head(6));
  v.segment(6, nb-1) = h*(x.dr + h*a.tail(nb-1));
  v.segment(nb+5, 6) = h*a.head(6);
  v.segment(nb+11, nb-1) = h*a.tail(nb-1);
}


void Mbs::Mass(MatrixXd &M, const MbsState &x) 
{
    
  // assumes that FK has been called on x
  M.setZero();

  vector<Matrix6d> Ics(nb);
  
  for (int i=0; i < nb; ++i) {
    Ics[i] = links[i].I.asDiagonal();
  }
  
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
    
    if (i > 0) {
      Vector6d F(Ics[i]*joints[i-1].S);
      M(6+i-1, 6+i-1) = joints[i-1].S.dot(F);
      
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
          M(6+i-1, 6+j-1) = F.dot(joints[j-1].S);
          M(6+j-1, 6+i-1) = M(6+i-1, 6+j-1);
        }
      }
      //        M.block<6,1>(0, 6+i-1) = F;
      //        M.block<1,6>(6+i-1, 0) = F.transpose();
      M.block(0, 6+i-1, 6, 1) = F;
      M.block(6+i-1, 0, 1, 6) = F.transpose();
    }
  }
  // M.topLeftCorner<6,6>() = Ics[0];
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
             double t, const MbsState &xb, const MbsState &xa,
             const VectorXd &u, double h) 
{
  
  // assume that Kstep was called so that all velocities are propagated
  
  vector<Vector6d> ps(nb);
  Matrix6d Da, Db;
  Matrix6d A;
  Matrix4d gi;

  VectorXd b(nb+5); // bias
  
  Vector3d ag(0, 0, -9.81);
  Vector6d gr;
  gr.setZero();

  for (int i = nb - 1; i >= 0; --i) {
    const Vector6d &va = xa.vs[i];
    const Vector6d &vb = xb.vs[i];
    const Vector6d &I = links[i].I;

    gr.tail<3>() = links[i].m*(xa.gs[i].topLeftCorner<3,3>().transpose()*ag);
    
    se3.dcayinv(Db, h*vb);
    se3.dcayinv(Da, -h*va);
    ps[i] = (Db.transpose()*(I.cwiseProduct(vb)) -  Da.transpose()*(I.cwiseProduct(va)))/h - gr;

    /*
    Vector6d mua, mub;

    se3.tln(mub, h*vb, I.cwiseProduct(vb));
    se3.tln(mua, -h*va, I.cwiseProduct(va));

    //    mub = I.cwiseProduct(vb);
    //    mua = I.cwiseProduct(va);    

    ps[i] = (mub - mua)/h;    
    */

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
      cout << "i=" << i << "ps=" << ps[i] << endl;
    
    if (i > 0)
      b[6+i-1] = ps[i].dot(joints[i-1].S);
    else
      b.head(6) = ps[0];
  }
  
  //  cout << "bias b=" << b << endl;

  VectorXd fu(nb + 5);
  Force(fu, t, xa, u);  // compute control/external forces

  //  cout << "control fu=" << fu << endl;

  f = b - fu;

  /*
  
  MatrixXd M(nb+5,nb+5);
  Mass(M, xa);

  VectorXd vdra(nb+5);
  VectorXd vdrb(nb+5);

  vdra.head(6) = xa.vs[0];  
  vdra.tail(nb-1) = xa.dr;  
  vdrb.head(6) = xb.vs[0];  
  vdrb.tail(nb-1) = xb.dr;  

  cout << "M=" << M << endl;

  f = M*(vdrb - vdra)/h + b - fu;

  */

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
  
  for (int i = 1; i < nb; ++i) { // take the difference and divide by time for each body
    se3.inv(g, xp.gs[i]);

    //    se3.cayinv(x.vs[i], g*x.gs[i]);
    se3.log(x.vs[i], g*x.gs[i]);

    x.vs[i] /= h;
  }
}


void Mbs::KStep(MbsState &xb, const MbsState &xa, double h) {
  
  // given: xa, xb.vs[0], xb.dr[*]
  // compute: xb.dgs[*], xb.gs[*],  xb.vs[*], xb.r[*]
  
  Matrix4d dg;
  Matrix4d gi;
  
  //    xb.Ms[0].setIdentity();
  
  for (int i = 0; i < nb; ++i) {
    if (i > 0) {
      xb.r[i-1] = xa.r[i-1] + h*xb.dr[i-1];
      
      Joint &jnt = joints[i-1];
      
      int pi = pis[i];
      se3.inv(gi, xa.gs[i]);
      se3.exp(dg, xb.r[i-1]*jnt.a);
      xb.dgs[i-1] = jnt.gp*dg*jnt.gci;  // this is also in FK
      
      //      se3.cayinv(xb.vs[i], gi*xb.gs[pi]*xb.dgs[i-1]);
      se3.log(xb.vs[i], gi*xb.gs[pi]*xb.dgs[i-1]);

      xb.vs[i] /= h;
    } 

    //    se3.cay(dg, h*xb.vs[i]);
    se3.exp(dg, h*xb.vs[i]);

    xb.gs[i] = xa.gs[i]*dg;
    
    //TODO: add else above
  }
}
