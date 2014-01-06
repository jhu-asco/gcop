#include <limits>
#include <iostream>
#include <utility>
#include "mbs.h"
#include <iomanip>


using namespace gcop;

Mbs::Mbs(int nb, int c) : nb(nb), 
                          System(*new MbsManifold(nb), *new Rn<>(c)), 
                          links(nb),
                          joints(nb-1), Ips(nb),
                          pis(nb), cs(nb), se3(SE3::Instance()),
                          method(TRAP),
                          iters(2),
                          debug(false), basetype("chainbase"),
                          damping(VectorXd::Zero(nb-1)),
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
	if(basetype == "chainbase")	 
	{
		assert(6 + nb-1 == f.size());
		f.head(6) = u.head(6);
	}
	else if(basetype == "airbase")
	{
		f[0] = u[0];
		f[1] = u[1];
		f[2] = u[2];
		f[3] = 0;
		f[4] = 0;
		f[5] = u[3];
	}

  f.tail(nb-1) = u.tail(nb-1) - damping.cwiseProduct(x.dr);

  if (A)
    A->setZero();
  if (B)
	{
		if(basetype == "chainbase")
			B->setIdentity();
		else if(basetype == "airbase")
		{
			B->setZero();
			(*B)(0,0) = 1;
			(*B)(1,1) = 1;
			(*B)(2,2) = 1;
			for(int count = 5;count < f.size();count++)
				(*B)(count,count-2) = 1;
		}
	}
}

void Mbs::Init() 
{
  Ips[0] = links[0].I.asDiagonal();
  
  for (int i = 0; i < nb-1; ++i) {    
    Joint &jnt = joints[i];
    jnt.Init();    
    Ips[i+1] = jnt.Ac.transpose()*links[i+1].I.asDiagonal()*jnt.Ac;
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
  VectorXd f(nb + 5);  // bias forces
  
  // compute bias
  ID(f, t, x, u);
  
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
  v.head(6) = h*(x.vs[0] + h*a.head(6));            // updated base velocity
  v.segment(6, nb-1) = h*(x.dr + h*a.tail(nb-1));   // updated joint velocities
  v.segment(nb+5, 6) = h*a.head(6);                 
  v.segment(nb+11, nb-1) = h*a.tail(nb-1);
}


void Mbs::Mass(MatrixXd &M, const MbsState &x) const 
{
  // assumes that FK has been called on x, i.e. that it is consistent
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
             double t, const MbsState &x, const VectorXd &u) 
{
  VectorXd b(nb + 5); // bias
  Bias(b, t, x);

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

  for (int i = nb - 1; i >= 0; --i) {
    const Vector6d &v = x.vs[i];
    const Vector6d &I = links[i].I;
    Vector6d mu = I.cwiseProduct(v);

    gr.tail<3>() = links[i].m*(x.gs[i].topLeftCorner<3,3>().transpose()*ag);
    
    ps[i].head(3) = v.head<3>().cross(mu.head<3>());
    ps[i].tail(3) = v.tail<3>().cross(mu.head<3>()) + v.head<3>().cross(mu.tail<3>());    
    ps[i] = ps[i] - gr;

    if (debug) {
      cout << "i=" << i << endl;
      cout << "v=" << v << endl;
      //      cout << "mua=" << mua << endl;
      //      cout << "mub=" << mub << endl;      
      cout << "ps[" << i << "]=" << ps[i] << endl;
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
    
    //      if (i > 0)
    //        ps[i] = joints[i-1].Ac.transpose()*ps[i];
    
    if (debug)
      cout << "i=" << i << "ps=" << ps[i] << endl;
    
    if (i > 0)
      b[6+i-1] = ps[i].dot(joints[i-1].S);
    else
      b.head(6) = ps[0];
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


void Mbs::Acc(VectorXd &a, double t, const MbsState& x, const VectorXd &u)
{
  MatrixXd M(nb + 5, nb + 5);
  Mass(M, x);

  VectorXd b(nb + 5); // bias
  Bias(b, t, x);

  VectorXd f(nb + 5);
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
  if (method == EULER)
    return System::Step(xb, t, xa, u, h, A, B);

  if (method == HEUN)
    return HeunStep(xb, t, xa, u, h, A, B);

  // initialize new velocity with old velocity
  VectorXd vdr(nb + 5);
  vdr.head(6) = xa.vs[0];
  vdr.tail(nb-1) = xa.dr;

  // Newton-Euler error
  VectorXd e(nb + 5);

  // Jacobian
  VectorXd en(nb + 5);
  VectorXd dvdr(nb + 5); 
  MatrixXd De(nb + 5, nb + 5);

  // search direction
  VectorXd d(nb + 5);

  double eps = 1e-6;
  for (int j = 0; j < iters; ++j) {    
    NE(e, vdr, xb, t, xa, u, h);      
    for (int i = 0; i < nb + 5; ++i) {
      dvdr.setZero();
      dvdr[i] = eps;
      NE(en, vdr + dvdr, xb, t, xa, u, h);    
      De.col(i) = (en - e)/eps;
    }    
    // take a Newton step
    //    vdr = vdr - De.inverse()*e;
    //    De.lu().solve(e, &d);
    //vdr = vdr - d;
    vdr = vdr - De.lu().solve(e);
  }  

  if (e.norm() > 1e-3) {
    cout << "[W] Mbs::Step: residual e seems high e=" << e << endl;
  }
  
  xb.vs[0] = vdr.head(6);
  xb.dr = vdr.tail(nb-1);
  KStep(xb, xa, h);
}



void Mbs::NE(VectorXd &e, const VectorXd &vdr, 
             MbsState &xb,               
             double t, const MbsState &xa, 
             const VectorXd &u, double h)
{
  xb.vs[0] = vdr.head(6);
  xb.dr = vdr.tail(nb-1);
  KStep(xb, xa, h);
  VectorXd b(nb + 5);

  DBias(b, t, xb, xa, h);

  VectorXd f(nb + 5);
  Force(f, t, xa, u);  // compute control/external forces

  e = b - f;
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
    
    se3.dcayinv(Db, h*vb);
    se3.dcayinv(Da, -h*va);
    ps[i] = (Db.transpose()*(I.cwiseProduct(vb)) -  Da.transpose()*(I.cwiseProduct(va)))/h - gr;

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
    
    //TODO: add else above for efficiency
  }
}
