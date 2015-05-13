#ifndef GCOP_SDDP_H
#define GCOP_SDDP_H
/** Creates a sampling based ddp. The algorithm provides a method to linearize the system based on sample trajectories around provided us.
*/

#include "docp.h"
#include <exception>      // std::exception
#include <stdexcept>
#include <random>

namespace gcop {
  
  using namespace std;
  using namespace Eigen;
  
  template <typename T, int nx = Dynamic, int nu = Dynamic, int np = Dynamic> class SDdp :
    public Docp<T, nx, nu, np> {
    
    typedef Matrix<double, nx, 1> Vectornd;
    typedef Matrix<double, nu, 1> Vectorcd;
    typedef Matrix<double, np, 1> Vectorpd;

    typedef Matrix<double, nx, nx> Matrixnd;
    typedef Matrix<double, nx, nu> Matrixncd;
    typedef Matrix<double, nu, nx> Matrixcnd;
    typedef Matrix<double, nu, nu> Matrixcd;

    //External render Function for trajectories:
    typedef void(RenderFunc)(int, vector<T>&);
    
  public:
    /**
     * Create an optimal control problem using a system, a cost, and 
     * a trajectory given by a sequence of times, states, and controls. 
     * The times ts must be given, the initial state xs[0] must be set, and
     * the controls us will be used as an initial guess for the optimization.
     *
     * After initialization, every call to Iterate() will optimize the 
     * controls us and states xs and modify them accordingly. Problems involving
     * time-optimization will also modify the sequence of times ts.
     * 
     * @param sys system
     * @param cost cost
     * @param ts (N+1) sequence of discrete times
     * @param xs (N+1) sequence of discrete states
     * @param us (N) sequence of control inputs
     * @param update whether to update trajectory xs using initial state xs[0] and inputs us.
     *               This is necessary only if xs was not already generated from us.
     */
    SDdp(System<T, nx, nu, np> &sys, Cost<T, nx, nu, np> &cost, 
         vector<double> &ts, vector<T> &xs, vector<Vectorcd> &us, 
         Vectorpd *p = 0, bool update = true);
    
    virtual ~SDdp();


    /**
     * Perform one DDP iteration. Internally calls:
     *
     * Backward -> Forward -> Update. The controls us and trajectory xs
     * are updated. 
     */
    void Iterate();

    /**
     *  Forward pass
     */
    void Forward();

    /**
     *  Backward pass
     */
    void Backward();

    /**
     * Linearize around existing us and xs by collecting samples
     */
    void Linearize();
        
    int N;        ///< number of discrete trajectory segments
    
    double mu;    ///< current regularization factor mu
    double mu0;   ///< minimum regularization factor mu
    double dmu0;  ///< regularization factor modification step-size
    double mumax; ///< maximum regularization factor mu
    double a;     ///< step-size
    double current_a; ///< current step-size
    
    std::vector<Vectorcd> dus;  ///< computed control change
    std::vector<T> xss;             ///< Sample states (N+1) vector
    std::vector<T> xsprev;          ///< Mean of Sample states (N+1) vector
    
    std::vector<Vectorcd> kus;
    std::vector<Vectorcd> du_sigma;///< Stores the stdeviation at each segment
    std::vector<Matrixcnd> Kuxs;
    std::vector<Vectorcd> Qud; ///< Inverse Variance for sampling du

		std::default_random_engine randgenerator; ///< Default random engine

    std::normal_distribution<double> normal_dist;///< Creates a normal distribution
    //std::uniform_real_distribution<double> uniform_dist;///< Creates a normal distribution

    Vectornd Lx;
    Matrixnd Lxx;
    Vectorcd Lu;
    Matrixcd Luu;
    
    Matrixnd P;
    Vectornd v;

    Vectorcd duscale;///<Scales sampled du according to this vector
    Vectornd dxscale;///<Scales sampled dx0 according to this vector
    RenderFunc *external_render;///<RenderFunction for rendering samples
    
    double V;
    Vector2d dV;

    double s1;   ///< Armijo/Bertsekas step-size control factor s1
    double s2;   ///< Armijo/Bertsekas step-size control factor s2
    double b1;   ///< Armijo/Bertsekas step-size control factor b1
    double b2;   ///< Armijo/Bertsekas step-size control factor b2
    
    char type;   ///< type of algorithm (choices are PURE, DDP, LQS), LQS is default. In the current implementation second-order expansion of the dynamics is ignored, so DDP and LQS are identical. Both LQS and DDP are implemented using the true nonlinear dynamics integration in the Forward step, while PURE uses the linearized version in order to match exactly the true Newton step. 
    int Ns; ///< Number of samples for linearizing
    
    static const char PURE = 0;     ///< PURE version of algorithm (i.e. stage-wise Newton)
    static const char DDP = 1;      ///< DDP version of algorithm
    static const char LQS = 2;      ///< Linear-Quadratic Subproblem version of algorithm due to Dreyfus / Dunn&Bertsekas / Pantoja
    
    /*
    class Fx : public Function<T, n, n> {
    public:
    Fx(System &sys) :
      Function<T>(sys.U.nspace), sys(sys) {
      }
      
      void F(Vectornd &f, const T &x) {

        
        sys.X.L
        this->xa = xa;
        sys.FK(this->xa);
        sys.ID(f, t, *xb, this->xa, *u, h);
      }
      
      Mbs &sys;
      
      double t;
      const MbsState *xb;
      const Vectorcd *u;
      double h;
    };

    class Fu : public Function<VectorXd> {
    public:
    Fu(Mbs &sys) :
      Function<VectorXd>(sys.U), sys(sys) {
      }
      
      void F(VectorXd &f, const VectorXd &u) {
        sys.ID(f, t, *xb, *xa, u, h);
      }
      
      Mbs &sys;
      
      double t;
      const MbsState *xa;
      const MbsState *xb;
      double h;
    };

    */

    bool pd(const Matrixnd &P) {
      LLT<Matrixnd> llt;     
      llt.compute(P);
      return llt.info() == Eigen::Success;
    }

    bool pdX(const MatrixXd &P) {
      LLT<MatrixXd> llt;     
      llt.compute(P);
      return llt.info() == Eigen::Success;
    }  

    void Reset()
    {
      this->current_a = 1.0;
      this->mu = (1e-3);
      this-> mu0 = (1e-3); 
      this->dmu0 = (2); 
      this->a = (1); 
      this->s1 = (0.1); 
      this->s2 = (0.5); 
      this->b1 = (0.25); 
      this->b2 = (2);
      for(int count = 0; count < N; count++)
      {
        Kuxs[count].setZero();
        kus[count].setZero();
      }
    }

  };

  using namespace std;
  using namespace Eigen;
  
  template <typename T, int nx, int nu, int np> 
    SDdp<T, nx, nu, np>::SDdp(System<T, nx, nu, np> &sys, 
                            Cost<T, nx, nu, np> &cost, 
                            vector<double> &ts, 
                            vector<T> &xs, 
                            vector<Matrix<double, nu, 1> > &us,
                            Matrix<double, np, 1> *p,
                            bool update) : Docp<T, nx, nu, np>(sys, cost, ts, xs, us, p, false),
    N(us.size()), 
    mu(1e-3), mu0(1e-3), dmu0(2), mumax(1e6), a(1), current_a(1),
    dus(N), kus(N), Kuxs(N), Qud(N), du_sigma(N),
    s1(0.1), s2(0.5), b1(0.25), b2(2),
    type(LQS),
    normal_dist(0,1),
    //uniform_dist(-1,1),
    external_render(0),
    xss(xs),
    xsprev(xs),
    Ns(30)
   // normal_dist(0,0.02)
    //normal_dist(0,0.001)
    {
      assert(N > 0);
      assert(sys.X.n > 0);
      assert(sys.U.n > 0);

      assert(ts.size() == N+1);
      assert(xs.size() == N+1);
      assert(xss.size() == N+1);
      assert(xsprev.size() == N+1);
      assert(us.size() == N);

      if (nx == Dynamic || nu == Dynamic) {
        for (int i = 0; i < N; ++i) {
          dus[i].resize(sys.U.n);
          //          As[i].resize(sys.X.n, sys.X.n);
          //          Bs[i].resize(sys.X.n, sys.U.n);
          kus[i].resize(sys.U.n);
          Kuxs[i].resize(sys.U.n, sys.X.n);
        }

        Lx.resize(sys.X.n);
        Lxx.resize(sys.X.n, sys.X.n);
        Lu.resize(sys.U.n);
        Luu.resize(sys.U.n, sys.U.n);

        P.resize(sys.X.n, sys.X.n);
        v.resize(sys.X.n);       
      }
      duscale = Vectorcd::Constant(0.02);//Default initialization of scale
      for(int i =0; i<N; ++i)
      {
        dus[i].setZero();
        Kuxs[i].setZero();//Set initial feedback matrices to zero
        kus[i].setZero();
        Qud[i] = Vectorcd::Constant(1.0/0.0004);//Set the inverse variance to the given value in the beginning
        du_sigma[i] = duscale;// Set the sigma to small value to begin with
      }
      dxscale = Vectornd::Constant(0.02);//Default initialization of scale
      //duscale.setZero();
      cout<<"duscale: "<<duscale<<endl;
      if (update) {
        this->Update(false);

     ++(this->nofevaluations);
        xss = xs;//Copy xs to xss in the beginning
        //this->Update(true);//#DEBUG
        //this->Linearize();
      }
    }
  
  template <typename T, int nx, int nu, int np> 
    SDdp<T, nx, nu, np>::~SDdp()
    {
    }
    
  template <typename T, int nx, int nu, int np> 
    void SDdp<T, nx, nu, np>::Backward() {

    typedef Matrix<double, nx, 1> Vectornd;
    typedef Matrix<double, nu, 1> Vectorcd;
    typedef Matrix<double, nx, nx> Matrixnd;
    typedef Matrix<double, nx, nu> Matrixncd;
    typedef Matrix<double, nu, nx> Matrixcnd;
    typedef Matrix<double, nu, nu> Matrixcd;  
    
    double t = this->ts.back();
    const T &x = this->xs.back();
    const Vectorcd &u = this->us.back();
    double L = this->cost.L(t, x, u, 0, 0, &Lx, &Lxx);
    
    V = L;
    dV.setZero();
    
    v = Lx;
    P = Lxx;
    
    Vectornd Qx;
    Vectorcd Qu;    
    Matrixnd Qxx;
    Matrixcd Quu;
    Matrixcd Quum;
    Matrixcnd Qux;   
    
    Matrixnd At;
    Matrixcnd Bt;
    
    Matrixcd Ic = MatrixXd::Identity(this->sys.U.n, this->sys.U.n);
    
    for (int k = N-1; k >=0; --k) {
      
      t = this->ts[k];
      double h = this->ts[k+1] - t;

      const T &x = this->xs[k];
      const Vectorcd &u = this->us[k];
      double L = this->cost.L(t, x, u, h, 0, &Lx, &Lxx, &Lu, &Luu);

      if (std::isnan(L)) {
        cout << "NAN" << " k=" << k << " Lx=" << Lx << " Lxx=" << Lxx << endl;
        getchar();
      }
      
      V += L;
      
      const Matrixnd &A = this->As[k];
      const Matrixncd &B = this->Bs[k];
      Vectorcd &ku = kus[k];
      Matrixcnd &Kux = Kuxs[k];
      
      At = A.transpose();    
      Bt = B.transpose();
      
      Qx = Lx + At*v;
      Qu = Lu + Bt*v;
      Qxx = Lxx + At*P*A;
      Quu = Luu + Bt*P*B;
      Qux = Bt*P*A;
      
      double mu = this->mu;
      double dmu = 1;
     
      if (this->debug) {
      if (!pd(P)) {
        cout << "P[" << k << "] is not pd =" << endl << P << endl;
        cout << "Luu=" << endl << Luu << endl;
      }

      Matrixcd Pp = Bt*P*B;
      if (!pdX(Pp)) {
        cout << "B'PB is not positive definite:" << endl << Pp << endl;
      }
      }

      LLT<Matrixcd> llt;     
      
      while (1) {
        Quum = Quu + mu*Ic;
        llt.compute(Quum);
        
        // if OK, then reduce mu and break
        if (llt.info() == Eigen::Success) {
          // this is the standard quadratic rule specified by Tassa and Todorov
          dmu = min(1/this->dmu0, dmu/this->dmu0);
          if (mu*dmu > this->mu0)
            mu = mu*dmu;
          else
            mu = this->mu0;
          
        if (this->debug) 
          cout << "[I] SDdp::Backward: reduced mu=" << mu << " at k=" << k << endl;

        //Store Diag value for inverse Variance:
        Qud[k] = Quum.diagonal();
        cout<<"Qud: "<<Qud[k].transpose()<<endl;


          break;
        }
        
        // if negative then increase mu
        dmu = max(this->dmu0, dmu*this->dmu0);
        mu = max(this->mu0, mu*dmu);
        
        if (this->debug) {
          cout << "[I] SDdp::Backward: increased mu=" << mu << " at k=" << k << endl;
          cout << "[I] SDdp::Backward: Quu=" << Quu << endl;
        }

        if (mu > mumax) {
          cout << "[W] SDdp::Backward: mu= " << mu << " exceeded maximum !" << endl;          
          if (this->debug)
            getchar();
          break;
        }
      }

      if (mu > mumax)
        break;
      
      ku = -llt.solve(Qu);
      Kux = -llt.solve(Qux);
      
      assert(!std::isnan(ku[0]));
      assert(!std::isnan(Kux(0,0)));

      v = Qx + Kux.transpose()*Qu;
      P = Qxx + Kux.transpose()*Qux;
      
      dV[0] += ku.dot(Qu);
      dV[1] += ku.dot(Quu*ku/2);
    }
    
    //    if (debug)
    cout << "[I] SDdp::Backward: current V=" << V << endl;    
  }

	template <typename T, int nx, int nu, int np> 
    void SDdp<T, nx, nu, np>::Forward() {
    //static int count_iterate = 0;
    //count_iterate++;

    typedef Matrix<double, nx, 1> Vectornd;
    typedef Matrix<double, nu, 1> Vectorcd;
    typedef Matrix<double, nx, nx> Matrixnd;
    typedef Matrix<double, nx, nu> Matrixncd;
    typedef Matrix<double, nu, nx> Matrixcnd;
    typedef Matrix<double, nu, nu> Matrixcd;  
    
    // measured change in V
    double dVm = 1;
    
    //double a = this->a;
    this->current_a = this->a;//Start with big step size
    
    while (dVm > 0) {

      Vectornd dx = VectorXd::Zero(this->sys.X.n);
			dx.setZero();//Redundancy
      T xn = this->xs[0];
      this->sys.Reset(xn,this->ts[0]);//Reset to initial state
      Vectorcd un;
      
      double Vm = 0;
      
      for (int k = 0; k < N; ++k) {
        const Vectorcd &u = this->us[k];
        Vectorcd &du = dus[k];
        
        const Vectorcd &ku = kus[k];
        const Matrixcnd &Kux = Kuxs[k];
        
        du = current_a*ku + Kux*dx;
        un = u + du; 
        
        Rn<nu> &U = (Rn<nu>&)this->sys.U;
        if (U.bnd) {
          for (int j = 0; j < u.size(); ++j) 
            if (un[j] < U.lb[j]) {
							//cout<<"I hit bound"<<endl;//[DEBUG]
              un[j] = U.lb[j];
              du[j] = un[j] - u[j];
            } else
              if (un[j] > U.ub[j]) {
                un[j] = U.ub[j];
                du[j] = un[j] - u[j];
              }
        }

        /*if(count_iterate == 1)
        {
          //[DEBUG]:
          cout<<"du[ "<<k<<"]:\t"<<du.transpose()<<endl;
          cout<<"dx[ "<<k<<"]:\t"<<dx.transpose()<<endl;
          cout<<"ku[ "<<k<<"]:\t"<<ku.transpose()<<endl;
          cout<<"Kux[ "<<k<<"]:"<<endl<<Kux<<endl;
        }
        */
        
        const double &t = this->ts[k];
        double h = this->ts[k+1] - t;

        double L = this->cost.L(t, xn, un, h);
        Vm += L;
        
        // update dx
        if (type == PURE) {
          dx = this->As[k]*dx + this->Bs[k]*du;
          this->sys.X.Retract(xn, xn, dx);
        } else {
          double h = this->ts[k+1] - t;
					//Adding nan catching :
					try
					{
						this->sys.Step(xn, un, h, this->p);
					}
					catch(std::exception &e)
					{
						//[DEBUG] Statement
						std::cerr << "exception caught: " << e.what() << '\n';
						std::cout<<" Iteration counter: "<<k<<endl;
						std::cout<<" u: "<<u.transpose()<<endl;
						std::cout<<" du: "<<du.transpose()<<endl;
						//More debug statements:
						//std::cout<<" a: "<<a<<"\t ku: "<<ku.transpose()<<"\t dx: "<<dx.transpose()<<"\n kux: \n"<<Kux<<endl;
						//[Culprit dx]

						throw std::runtime_error(std::string("Nan observed"));
						return;
					}
					//cout<<"dx: "<<dx.transpose()<<endl;//[DEBUG]//Previous iteration dx
					//std::cout<<" du: "<<du.transpose()<<endl;//[DEBUG]
          this->sys.X.Lift(dx, this->xs[k+1], xn);
					/*if((k == 29) || (k == 30) || (k == 31))
					{
						//cout dx:
						cout<<"dx[ "<<(k+1)<<"]:\t"<<dx.transpose()<<endl;
						cout<<"un[ "<<(k+1)<<"]:\t"<<un.transpose()<<endl;
						cout<<"du[ "<<k<<"]:\t"<<du.transpose()<<endl;
						cout<<"Printing states ["<<(k+1)<<"]"<<endl;
						this->sys.print(this->xs[k+1]);
						cout<<"Printing xn: "<<endl;
						this->sys.print(xn);
						if(k == 31)
						{
							//exit(0);
						}
					}
					*/
					//cout<<"xs[ "<<(k+1)<<"]:\t"<<this->xs[k+1]<<endl;
					//cout<<"un[ "<<(k+1)<<"]:\t"<<un.transpose()<<endl;
          
          //          cout << xn.gs[0] << " " << xn.r << " " << xn.vs[0] << " " << xn.dr << endl;

          //          cout << this->xs[k+1].gs[0] << " " << this->xs[k+1].r << " " << this->xs[k+1].vs[0] << " " << this->xs[k+1].dr << endl;
          assert(!std::isnan(dx[0]));
        }
      }
      ++(this->nofevaluations);
      
      double L = this->cost.L(this->ts[N], xn, un, 0);
      Vm += L;
      
      if (this->debug)
        cout << "[I] SDdp::Forward: measured V=" << Vm << endl;
      
      dVm = Vm - V;
      
      if (dVm > 0) {
       current_a *= b1;
       if(current_a == 0)
         break;
        if (current_a < 1e-12)
        {
          current_a = 0;
        }
        if (this->debug)
          cout << "[I] SDdp::Forward: step-size reduced a=" << current_a << endl;
        
        continue;
      }
      
      /*double r = dVm/(current_a*dV[0] + current_a*current_a*dV[1]);
      if (r < s1)
       current_a = b1*current_a;
      else 
        if (r >= s2) 
         current_a = b2*current_a;    
    
      if (this->debug)
        cout << "[I] SDdp::Forward: step-size a=" << current_a << endl;    
        */
    }
    this->J = V + dVm;//Set the optimal cost after one iteration
  }

  template <typename T, int nx, int nu, int np> 
    void SDdp<T, nx, nu, np>::Iterate() {
    
    Linearize();//Linearize with the current us
    Backward();
    Forward();
    for (int k = 0; k < N; ++k)
      this->us[k] += dus[k];
    this->Update(false);
    ++(this->nofevaluations);
    //this->Update(true);//#DEBUG
  }
 template <typename T, int nx, int nu, int np>
    void SDdp<T, nx, nu, np>::Linearize(){
      //static int count_iterate = 0;
			randgenerator.seed(370212);
      MatrixXd dusmatrix(nu*N,Ns);
      MatrixXd dxsmatrix(nx*(N+1),Ns);

      //Vectorcd du;
      Vectorcd us1;

      int count = 0;

      //Create dus random matrix:
//#pragma omp parallel for private(count)
      /*for(count = 0;count < Ns;count++)
      {
        for(int count1 = 0;count1 < N*nu;count1++)
        {
          //dusmatrix(count1,count)  = normal_dist(randgenerator);
          int count_u = count1%nu;//find the index corresponding to du_scale
          dusmatrix(count1,count)  = duscale(count_u)*uniform_dist(randgenerator);
        }
      }
      */
      //cout<<dusmatrix<<endl;//#DEBUG
      //getchar();

      /*Find xsprev based on dus and us:
      this->sys.Reset(this->xs[0],this->ts[0]);
      for(int count_traj = 0; count_traj< N;count_traj++)
      {
        us1 = this->us[count_traj] - this->dus[count_traj];
        this->sys.Step(xsprev[count_traj+1],us1,(this->ts[count_traj+1])-(this->ts[count_traj]), this->p);
      }
      */

      //dxsmatrix.block(0,0,nx,Ns).setZero();//Set zero first block
      Vectornd dx;
      T x0_sample;
      Rn<nu> &U = (Rn<nu>&)this->sys.U;
      //NonParallelizable code for now
      // First do some iterations to adjust the dus_scale(stdeviation of du) so that dx_stdeviation lies close to dx_scale (target stdeviation) across trajectory
      //Set all du_sigma to given default values:
      for(count = 0; count < N; count++)
      {
        du_sigma[count] = duscale;//Set to given default value
      }
      //double kp = 0.005;
      /*
      double kp = 0.002;
      int Ns1 = 15;
      vector<double> dxsquare(N+1, 0);//0 to N
      Vectornd dxscale1 = dxscale;
      for(int count_it = 0;count_it < 2;count_it++)
      {
        for(count = 0;count < Ns1; count++)
        {
          //Adding Variance to duscale to ensure common variance over the trajectory///////////////////:
          // First Sample:
          for(int count1 = 0; count1 < nx; count1++)
          {
            dx(count1) = dxscale1(count1)*normal_dist(randgenerator);//Adjust dx_scale 
          }
          //dxsmatrix.block<nx,1>(0,count) = dx; 
          this->sys.X.Retract(x0_sample, this->xs[0], dx);
          xss[0] = x0_sample;
          this->sys.Reset(x0_sample,this->ts[0]);
          for(int count_traj = 0; count_traj< N;count_traj++)
          {
            this->sys.X.Lift(dx, this->xs[count_traj], xss[count_traj]);//This is for feedback 
            dxsquare[count_traj] += dx.squaredNorm();
            //cout<<"dx: "<<dx.transpose()<<endl;
            //if(du_sigma[count_traj].maxCoeff() > 0.00011)//Only Use Feedback if variance is atleast 0.002
            {
              us1 = this->us[count_traj] + this->Kuxs[count_traj]*dx;//Before Update
            }
            */
            /*else
            {
              us1 = this->us[count_traj];//No Feedback
            }
            */
            /*for(int count_u = 0;count_u < nu; count_u++)
            {
              */
              /*if(du_sigma[count_traj][count_u] < 0.002)
                du_sigma[count_traj][count_u] = 0.002;
              else if(du_sigma[count_traj][count_u] > 0.2)
                du_sigma[count_traj][count_u] = 0.2;//Bounds ond du_sigma
                */
                /*
              us1[count_u] = us1[count_u] + du_sigma[count_traj][count_u]*normal_dist(randgenerator);
            }
            //cout<<"du_sigma: "<<du_sigma[count_traj].transpose()<<endl;
            //getchar();
            this->sys.Step(xss[count_traj+1],us1,(this->ts[count_traj+1])-(this->ts[count_traj]), this->p);
          }
          this->sys.X.Lift(dx, this->xs[N], xss[N]);//This is for feedback
          dxsquare[N] += dx.squaredNorm();
          //cout<<(dx.squaredNorm()/dxscale.squaredNorm()-1)<<endl;
        }
        cout<<(dxsquare[0]/(Ns1*dxscale1.squaredNorm()) - 1)<<endl;
        for(count = 1;count <= N; count++)
        {
          for(int count_sigma = 0;count_sigma < count; count_sigma++)
          {
            du_sigma[count_sigma] -= (kp*(count_sigma+1)/double(count))*Vectorcd::Constant(dxsquare[count]/(Ns1*dxscale.squaredNorm()) - 1);
            //cout<<"du_sigma, dx-dxscale["<<count_sigma<<"]: "<<du_sigma[count_sigma].transpose()<<endl;
            for(int count_u = 0;count_u < nu; count_u++)
            {
              if(du_sigma[count_sigma][count_u] < 0.005)
                du_sigma[count_sigma][count_u] = 0.005;
              else if(du_sigma[count_sigma][count_u] > 0.2)
                du_sigma[count_sigma][count_u] = 0.2;//Bounds ond du_sigma
            }
          }
          //dxscale1 -= (kp/double(count))*Vectornd::Constant(dxsquare[count]/(Ns1*dxscale.squaredNorm()) - 1);
          cout<<(dxsquare[count]/(Ns1*dxscale.squaredNorm()) - 1)<<endl;
          //set dxsquare back to zero:
          dxsquare[count] = 0;
        }
        dxsquare[0] = 0;
        for(count = 0; count< N; count++)
        {
          cout<<"du_sigma["<<count<<"]: "<<du_sigma[count].transpose()<<endl;
        }
        //cout<<(dx.squaredNorm()/dxscale.squaredNorm() - 1)<<endl;
        cout<<"------Iteration done----"<<endl;
      }
        */
      //getchar();
      for(count = 0;count < Ns;count++)
      {
        //Set to initial state perturbed by small amount:
        for(int count1 = 0; count1 < nx; count1++)
        {
          dx(count1) = this->dxscale(count1)*normal_dist(randgenerator);//Adjust dx_scale 
        }
        dxsmatrix.block<nx,1>(0,count) = dx; 
        //cout<<"dx0: "<<dx.transpose()<<endl;
        //cout<<dx.norm();

        this->sys.X.Retract(x0_sample, this->xs[0], dx);

        xss[0] = x0_sample;
        this->sys.Reset(x0_sample,this->ts[0]);
        for(int count1 = 0;count1 < N;count1++)
        {
          //this->sys.X.Lift(dx, xsprev[count1], xss[count1]);//This is for feedback
          //us1 = this->us[count1] - this->dus[count1] + (this->current_a)*this->kus[count1] + this->Kuxs[count1]*dx;//Before Update
          this->sys.X.Lift(dx, this->xs[count1], xss[count1]);//This is for feedback
          dxsmatrix.block<nx,1>((count1)*nx,count) = dx; 
          //if(du_sigma[count1].maxCoeff() > 0.00011)//Only Use Feedback if variance is greater than 0.0001(min)
          {
            us1 = this->us[count1] + this->Kuxs[count1]*dx;//Before Update
            //us1 = this->us[count1];//No Feedback
            //cout<<"Feedback: "<< (this->Kuxs[count1]*dx).transpose()<<endl;
          }
          /*else
          {
            us1 = this->us[count1];//No Feedback
          }
          */
         // if(count_iterate != 2)
          {
            for(int count_u = 0;count_u < nu; count_u++)
            {
              //Make sure Qud is not smaller than some threshold:
              //dusmatrix(count1,count)  = duscale(count_u)*uniform_dist(randgenerator);
              //if(Qud[count1][count_u]<(1.0/(duscale[count_u]*duscale[count_u])))
              assert((du_sigma[count1][count_u] < 0.21)&&(du_sigma[count1][count_u]>0.00009));
              us1[count_u] = us1[count_u] + du_sigma[count1][count_u]*normal_dist(randgenerator);
              //else
                //us1[count_u] = us1[count_u] + (1.0/sqrt(Qud[count1][count_u]))*normal_dist(randgenerator);
            }
            //getchar();
          }
          //The sampled control should be within the control bounds of the system !!!
          if (U.bnd) {
            for (int count_u = 0; count_u < nu; ++count_u)
              if (us1[count_u] < U.lb[count_u]) {
                //cout<<"I hit bound"<<endl;//[DEBUG]
                us1[count_u] = U.lb[count_u];
             //   du[j] = un[j] - u[j];
              } else
                if (us1[count_u] > U.ub[count_u]) {
                  us1[count_u] = U.ub[count_u];
                  //du[j] = un[j] - u[j];
                }
          }

          dusmatrix.block<nu,1>(count1*nu, count) = us1 - this->us[count1];//Verify that this is zero without randomness #DEBUG
          //dumatrix.block<nu,1>(count1*nu,count) = us1;
          /*if(count_iterate == 2)
          {
            cout<<"us1, us: "<<us1.transpose()<<"\t , "<<(this->us[count1]).transpose()<<endl;
            //cout<<"dus: "<<(this->dus[count1]).transpose()<<"\t"<<.transpose()<<endl;
            cout<<"dx[ "<<count1<<"]:\t"<<dx.transpose()<<endl;
            cout<<"ku[ "<<count1<<"]:\t"<<this->kus[count1].transpose()<<endl;
            cout<<"Kux[ "<<count1<<"]:"<<endl<<this->Kuxs[count1]<<endl;
            cout<<"current_a"<<current_a<<endl;
            //cout<<"ku, Kxu: "<<endl<<(this->kus[count1])<<endl<<(this->Kuxs[count1])<<endl;
          }
          */
          //cout<<"du, u: "<<du.transpose()<<"\t"<<(this->us[count1]).transpose()<<endl;
          this->sys.Step(xss[count1+1],us1,(this->ts[count1+1])-(this->ts[count1]), this->p);
          //this->sys.X.Lift(dx, this->xs[count1+1], xss[count1+1]);//This is for fitting Linear Model
          //cout<<" "<<dx.norm();
          /*if(count_iterate == 2)
          {
            cout<<"dx: "<<count1<<"\t"<<dx.transpose()<<endl; 
            //getchar();
          }
          */
        }
        ++(this->nofevaluations);
        //Final step:
        this->sys.X.Lift(dx, this->xs[N], xss[N]);//This is for feedback
        dxsmatrix.block<nx,1>((N)*nx,count) = dx; 
        //cout<<endl; 

        //Render trajectory samples if external rendering function is provided:
        if(external_render)
        {
          external_render(count,xss);//ID for the sample trajectory
        }
        //if(count_iterate == 2)
          //getchar();
      }

      //Matrix<double, nx, nx+nu>Abs;
      //cout<<dxsmatrix<<endl;//#DEBUG
      //getchar();
      MatrixXd Abs(nx+nu,nx);
      MatrixXd XUMatrix(nx+nu,Ns);

      MatrixXd Xverifymatrix(nx,Ns);//Verify 

      //Compute As and Bs for every k:
      for(count = 0;count < N;count++)
      {
        XUMatrix<<dxsmatrix.block(count*nx,0,nx,Ns), dusmatrix.block(count*nu,0,nu,Ns);
        Abs = (XUMatrix*XUMatrix.transpose()).ldlt().solve(XUMatrix*dxsmatrix.block((count+1)*nx,0,nx,Ns).transpose());
        // Compare As with true As #DEBUG: 
        /*cout<<"Computed: As, Bs: "<<endl;
        cout<<endl<<Abs.block<nx,nx>(0,0).transpose()<<endl;
        cout<<endl<<Abs.block<nu,nx>(nx,0).transpose()<<endl;
        cout<<"True As, bs: "<<endl;
        cout<<endl<<this->As[count]<<endl;
        cout<<endl<<this->Bs[count]<<endl;
        cout<<"Diff As, bs: "<<endl;
        cout<<endl<<(Abs.block<nx,nx>(0,0).transpose() - this->As[count])<<endl;
        cout<<endl<<(Abs.block<nu,nx>(nx,0).transpose() - this->Bs[count])<<endl;
        */
//        getchar(); 

        this->As[count] = Abs.block<nx,nx>(0,0).transpose();
        this->Bs[count] = Abs.block<nu,nx>(nx,0).transpose();

        /******VERIFY**********/
        Xverifymatrix = Abs.transpose()*XUMatrix;
        //cout<<endl<<"dXpredicted: "<<endl<<Xverifymatrix<<endl;
        //cout<<endl<<"True_dX: "<<dxsmatrix.block((count+1)*nx,0,nx,Ns)<<endl;
        cout<<endl<<"Error_predicted: "<<(Xverifymatrix - dxsmatrix.block((count+1)*nx,0,nx,Ns)).squaredNorm()<<endl;
        //cout<<endl<<Abs<<endl;
        //getchar();

        //Debug: Print:
/*        cout<<"As, Bs, Abs: "<<endl;
        cout<<endl<<Abs<<endl;
        cout<<endl<<this->As[count]<<endl;
        cout<<endl<<this->Bs[count]<<endl;
        cout<<"Comparing dxs[count+1] with As*dxs[count] + Bs[count]*du forall samples "<<endl;
        cout<<"Prediction - Actual: "<<endl;
        cout<<endl<<(Abs.transpose()*XUMatrix - dxsmatrix.block((count+1)*nx,0,nx,Ns))<<endl;
        getchar();
        */
      }
      //count_iterate++;
    }
}

#endif
