#ifndef GCOP_GP_H
#define GCOP_GP_H

#include <Eigen/Dense>

namespace gcop {
  class GP {
  public:
    GP(int d, int n);
    
    virtual ~GP();
    
    virtual void Sample();
    
    bool Add(const Eigen::VectorXd &x,
             double f);    
    
    void Train();
    
    void Train(const Eigen::MatrixXd &Xs, 
               const Eigen::VectorXd &fs);
    
    double Predict(const Eigen::VectorXd &x, 
                   double *s = 0) const;
    
    double SqExp(const Eigen::VectorXd &xa, 
                 const Eigen::VectorXd &xb) const;
    
    double LogL(double dll[2] = 0);

    double PI(const Eigen::VectorXd &x, double fmin) const;
    

    double OptParams();
    
    int d;
    int n;
    
    Eigen::MatrixXd Xs;
    Eigen::VectorXd fs;
    
    Eigen::MatrixXd K;
    Eigen::MatrixXd Ki;
    
    Eigen::MatrixXd L;
    
    Eigen::VectorXd a;
    
    double l;
    double s;
    
    double sigma;
    
    bool cf;   ///< propagate cholesky factor L rather than K^{-1}

    bool eps;  ///< prohibit adding points that are eps-close in L_2 to existing data 
    
  };
}

#endif
