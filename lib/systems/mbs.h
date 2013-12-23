#ifndef GCOP_MBS_H
#define GCOP_MBS_H

#include "body3d.h"
#include "joint.h"
#include "mbsmanifold.h"
#include "function.h"
#include "mbscspace.h"
#include "mbstspace.h"

namespace gcop {
  
  using namespace std;
  using namespace Eigen;
  
  /**
   * The discrete mechanics is
   * based on an implicit symplectic variational integrator which is second-order accurate
   * and momentum-balance preserving. Regularity is guaranteed by choosing a time-step 
   * which does not jump out of the convex region around the initial guess from the previous
   * dynamics iteration.
   *
   * Author: Marin Kobilarov marin(at)jhu.edu
   */
  class Mbs : public System< MbsState, VectorXd> {
    
    /*
    typedef Matrix<double, c, 1> Vectorcd;
    typedef Matrix<double, c, c> Matrixcd;
    typedef Matrix<double, 6, c> Matrix6xcd;
    typedef Matrix<double, MBS_DIM(nb), c> Matrixncd;
    typedef Matrix<double, MBS_DIM(nb), MBS_DIM(nb)> Matrixnd;
    
    typedef Matrix<double, 6 + nb - 1, 1> Vectormd;
    */
    //    typedef pair<Matrix4d, Vector<double, 2*m - 6> > MbsState;
    
  public:
    /*
    class Fq : public Function<MbsState> {
    public:
    Fq(Mbs &sys) :
      Function<MbsState>(sys.cspace), sys(sys), xa(sys.nb) {
      }
      
      void F(VectorXd &f, const MbsState &xa) {
        this->xa = xa;
        sys.FK(this->xa);
        sys.ID(f, t, *xb, this->xa, *u, h);
      }
      
      Mbs &sys;
      
      double t;
      MbsState xa;
      const MbsState *xb;
      const Vectorcd *u;
      double h;
    };

    class Fva : public Function<MbsState> {
    public:
    Fva(Mbs &sys) :
      Function<MbsState>(sys.tspace), sys(sys), x(sys.nb), xa(sys.nb) {        
      }
      
      void F(VectorXd &f, const MbsState &xa) {
        this->xa = xa;
        sys.Rec(this->xa, h);        
        sys.ID(f, t, *xb, this->xa, *u, h);
      }
      
      Mbs &sys;

      MbsState x;      ///< previous state before xa
      double t;      
      MbsState xa;     ///< temporary object for xa
      const MbsState *xb;    ///< pointer to xb
      const Vectorcd *u;         ///< pointer to u
      double h;
    };    

    class Fvb : public Function<MbsState> {
    public:
    Fvb(Mbs &sys) :
      Function<MbsState>(sys.tspace), sys(sys), xb(sys.nb) {  
      }
      
      void F(VectorXd &f, const MbsState &xb) {
        this->xb = xb;
        sys.KStep(this->xb, *xa, h);
        sys.ID(f, t, this->xb, *xa, *u, h);
      }
      
      Mbs &sys;
      
      double t;
      const MbsState *xa;
      MbsState xb;
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

    Mbs(int nb, int c);

    virtual ~Mbs();

    void Init();

    /*
      double Step(MbsState &xb, double t, const MbsState &xa, 
      const VectorXd &u, double h,
      MatrixXd *A = 0, MatrixXd *B = 0);
    */

    double F(VectorXd &v, double t, const MbsState &xa, 
             const VectorXd &u, double h,
             MatrixXd *A = 0, MatrixXd *B = 0);


    /**
     * Inverse-dynamics function. Computes the residual of the
     * dynamical update from xa to xb given controls u. This is regarded as
     * the discrete-time equivalent of (b - B*u)
     * @param f residual
     * @param t time
     * @param x state
     * @param u controls
     */
    void ID(VectorXd &f,
            double t, const MbsState &x,
            const VectorXd &u);

    /**
     * Compute total bias b, i.e. such that M*a + b = B*u
     * @param b bias
     * @param t time
     * @param x state
     */
    void Bias(VectorXd &b,
              double t, const MbsState &x) const;

    /**
     * Discrete bias
     * @param b bias
     * @param t time
     * @param xb state
     * @param xa state
     * @param h time-step
     */
    void DBias(VectorXd &b,
                double t,
                const MbsState &xb, 
                const MbsState &xa, double h);    
    
    /**
     * Forward kinematics: given base pose and joint angles, update the poses of
     * all bodies in the graph.
     * @param x state
     */
    void FK(MbsState &x);

    /**
     * Kinematic-step. Given previous state xa, the next state xb is updated
     * by taking the joint angles from xa, and joint velocities from xb, updating the
     * joint angles in xb, and updating the xb body poses
     * @param xa previous state
     * @param xb next state
     * @param h time-step
     */
    void KStep(MbsState &xb, const MbsState &xa, double h);

    /**
     * Compute the mass matrix at a given state x
     * @param M mass matrix
     * @param x state
     */
    void Mass(MatrixXd &M, const MbsState &x) const;

    /**
     * Total resulting force on the system from external (e.g. gravity)
     * and internal (control) inputs
     * @param f total force f(x,u,t)
     * @param t time
     * @param x state
     * @param u controls
     * @param A jacobian Dxf
     * @param B jacobian Duf
     */
    virtual void Force(VectorXd &f, double t, const MbsState &x, const VectorXd &u,MatrixXd *A = 0, MatrixXd *B = 0);
    

    void Rec(MbsState &x, double h);

    int nb;                   ///< number of rigid bodies

    vector<Body3d<>> links;   ///< links
    vector<Joint> joints;     ///< joints

    vector<Matrix6d> Ips;     ///< A'*I*A (nb) vector
    
    vector<int> pis;          ///< parent indexes
    vector<vector<int> > cs;  ///< child lists

    Vector3d ag;             ///< acceleration due to gravity (0, 0, -9.81) by default

    
    SE3 &se3;                 ///< singleton reference for performing SE(3) operations

    bool debug;


    // below is for autodiff stuff
    //    MbsManifold X;
    //    Rn<> U;
    //    MbsCspace cspace;
    //    MbsTspace tspace;


    //    Fq fq;
    //    Fva fva;
    //    Fvb fvb;
    //    Fu fu;
  };
}

#endif
