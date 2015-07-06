#ifndef GCOP_UNIFORMSPLINETPARAM_H
#define GCOP_UNIFORMSPLINETPARAM_H

#include "tparam.h"

namespace gcop {
  
  using namespace std;
  using namespace Eigen;

  /**
   * General trajectory parametrization functionality. A trajectory can be parametrizes using
   * a set of discrete controls, a continuous control paramertization, a set of discrete states, 
   * a set of discrete flat outputs, etc...
   *
   * This version implements the Uniform spline interpolation of the data. It also implements least square fitting of the given control data using spline
   * This is indigenous support which supports only upto degree = 4.
   * Also it assumes tks has atleast 3 points i.e minimum of 2 segments
   *
   * Author: Marin Kobilarov (c) 2005--2013
   * Author2: Gowtham Garimella
   */
  template <typename T, 
    int nx = Dynamic, 
    int nu = Dynamic,
    int np = Dynamic> class UniformSplineTparam : public Tparam<T, nx, nu, np> {
    
    typedef Matrix<double, nx, 1> Vectornd;
    typedef Matrix<double, nu, 1> Vectorcd;
    typedef Matrix<double, nx, nx> Matrixnd;
    typedef Matrix<double, nx, nu> Matrixncd;
    typedef Matrix<double, nu, nx> Matrixcnd;
    typedef Matrix<double, nu, nu> Matrixcd;  

    typedef Matrix<double, np, 1> Vectormd;

  public:
    UniformSplineTparam(System<T, nx, nu, np> &sys, const VectorXd &tks, int degree = 2);//By default use quadratic spline
    
    bool To(VectorXd &s, 
            const vector<double> &ts, 
            const vector<T> &xs, 
            const vector<Vectorcd> &us,
            const Vectormd *p = 0);
    
    bool From(vector<double> &ts, 
              vector<T> &xs, 
              vector<Vectorcd> &us,
              const VectorXd &s,
              Vectormd *p = 0);
    
    const VectorXd &tks;
    int degree; //Degree of the spline // p = (m - n - 1) where m is knot vector size; n is the control vector size (tks size)

    //Basis Function:
    void FindBasis(VectorXd &basis, double u)
    {
      assert((u <=1) && (u>=0));
      assert(degree <= 4 && degree>=1);//I only have closed form solution for basis from Wiki http://www.wikiwand.com/en/Irwin%E2%80%93Hall_distribution#/Special_cases
      //Can put an assert for basis size not doing right now
      //In general it is easy to use deboor algorithm to find the basis for any n.
      double x;
      switch(degree)
      {
        case 1:
          basis[1] = u;
          x = 1+u;
          basis[0] = 2-x;
          break;
        case 2:
          basis[2] = 0.5*u*u;
          x = 1+u;
          basis[1] = 0.5*(-2*x*x + 6*x -3);
          x = 2+u;
          basis[0] = 0.5*(x*x - 6*x +9);
          break;
        case 3:
          basis[3] = (1.0/6.0)*u*u*u;
          x = 1+u;
          basis[2] = (1.0/6.0)*(-3*x*x*x + 12*x*x -12*x + 4);
          x = 2+u;
          basis[1] = (1.0/6.0)*(3*x*x*x - 24*x*x +60*x - 44);
          x = 3+u;
          basis[0] = (1.0/6.0)*(-1*x*x*x + 12*x*x -48*x +64);
          break;
        case 4:
          basis[4] = (1.0/24.0)*u*u*u*u;
          x = 1+u;
          basis[3] = (1.0/24.0)*(-4*x*x*x*x +20*x*x*x - 30*x*x + 20*x -5);
          x = 2+u;
          basis[2] = (1.0/24.0)*(6*x*x*x*x -60*x*x*x + 210*x*x - 300*x +155);
          x = 3+u;
          basis[1] = (1.0/24.0)*(-4*x*x*x*x +60*x*x*x - 330*x*x + 780*x -655);
          x = 4+u;
          basis[0] = (1.0/24.0)*(x*x*x*x -20*x*x*x + 150*x*x - 500*x +625);
          break;
          ////ADD MORE AS U GO////
      }
    }
  };
  
  template <typename T, int nx, int nu, int np> 
    UniformSplineTparam<T, nx, nu, np>::UniformSplineTparam(System<T, nx, nu, np> &sys, const VectorXd &tks, int degree) :  Tparam<T, nx, nu, np>(sys, (tks.size()+degree-1)*(sys.U.n)),degree(degree), tks(tks) {

      assert(tks.size() >=2);//atleast 1 segments for which above constructor size is valid
      //cout<<"Degree: "<<degree<<endl;
      
      //cout<<"tks: "<<tks.transpose()<<endl;
  }

  template <typename T, int nx, int nu, int np> 
    bool UniformSplineTparam<T, nx, nu, np>::To(VectorXd &s, 
                                             const vector<double> &ts, 
                                             const vector<T> &xs, 
                                             const vector<Vectorcd> &us,
                                             const Vectormd *p) {
      s.resize(this->ntp);
      //Find control points(uks) given us
      int nofcontrolpoints = (degree + this->tks.size() - 1);
      int tks_index = 0;

      //Initialize variables for least square fit; Fitting each control dimension separately
      MatrixXd A(us.size(), nofcontrolpoints);//Basis Matrix
      MatrixXd Asquare(nofcontrolpoints, nofcontrolpoints);//Basis Matrix
      VectorXd c(us.size());//RHS 
      A.setZero();//Initialize A
      VectorXd basis(degree+1);

      //Find Basis Matrix to find the control points:
      for(int ind = 0;ind < us.size(); ind++)
			{
				//Find the region in which ts[i] lies:
				if(tks_index < (tks.size()-1))
					while((ts[ind] - tks[tks_index+1])>1e-17)
						tks_index +=1;
				//tks_index = (ts[ind] - tks[tks_index+1])<-1e-17?tks_index:(tks_index+1);//1e-17 so that it considers numerical accuracy
				//Find the Basis for the normalized coordinate:
				//cout<<"index: "<<tks_index<<endl;
				double u = (ts[ind] - tks[tks_index])/(tks[tks_index+1] - tks[tks_index]);
				//cout<<"tks["<<(tks_index)<<"]: "<<tks[tks_index]<<"\t"<<tks[tks_index+1]<<endl;
				//cout<<"ts: "<<ts[ind]<<endl;
				//cout<<"u: "<<u<<endl;
				FindBasis(basis, u);
				A.block(ind, tks_index, 1, degree+1) = basis.transpose();//Create Basis Matrix
			}
      //Asquare = (A.transpose()*A);
      //cout<<"A: "<<endl<<A<<endl;
      //cout<<"Asquare: "<<endl<<Asquare<<endl;
      //getchar();
      for(int ind = 0; ind < (this->sys.U.n); ind++)
      {
        for(int uind = 0; uind < us.size(); uind++)
        {
          c(uind) = us[uind](ind);
        }
        Asquare = A.transpose()*A;//Weighted matrix
        VectorXd b = Asquare.ldlt().solve(A.transpose()*c);
        cout<<"Error: "<<(A*b - c).squaredNorm()<<endl;
      //  cout<<"c: "<<c.transpose()<<endl;
        //cout<<"b: "<<b.transpose()<<endl;
        //Copy the elements back into vector s:
        for(int sind = 0; sind < nofcontrolpoints; sind++)
        {
          s(sind*(this->sys.U.n) + ind) = b(sind);
        }
      }
      cout<<"s: "<<s.transpose()<<endl;
      //getchar();
      //Verify if this s is good:
      /*
      {
        tks_index = 0;
        Vectorcd usi;
        for (int i = 0; i < us.size(); ++i) {
          //Find the region in which ts[i] lies:
          if(tks_index < (tks.size()-1))
            tks_index = (ts[i] - tks[tks_index+1])<-1e-17?tks_index:(tks_index+1);
          //Find the Basis for the normalized coordinate:
          double u = (ts[i] - tks[tks_index])/(tks[tks_index+1] - tks[tks_index]);
          FindBasis(basis, u);
          usi.setZero();
          //Weight the spline using the basis:
          for(int degree_count = 0; degree_count <= degree; degree_count++)
          {
            usi += basis[degree_count]*s.segment((tks_index+degree_count)*nu,nu);
          }
          cout<<"Basis: "<<basis.transpose()<<endl;

          cout<<"us_pred["<<i<<"]: "<<usi.transpose()<<"us_act: "<<us[i].transpose()<<endl;
        }
      }
      getchar();
      */

      //s.setZero();
      return true;
    }
  
  template <typename T, int nx, int nu, int np> 
    bool UniformSplineTparam<T, nx, nu, np>::From(vector<double> &ts, 
                                                  vector<T> &xs, 
                                                  vector<Vectorcd> &us,
                                                  const VectorXd &s,
                                                  Vectormd *p) {
    assert(s.size() == (degree + this->tks.size() - 1)*(this->sys.U.n));
    //cout<<"s.size(): "<<s.size()<<endl;
    //cout<<"tks.size(): "<<tks.size()<<endl;
    //cout<<"degree: "<<degree<<endl;
    //cout<<"sys.U.n: "<<(this->sys.U.n)<<endl;
    //cout<<"s: "<<s.transpose()<<endl;
    //getchar();

    this->sys.Reset(xs[0],ts[0]);
    int tks_index = 0;
    VectorXd basis(degree+1);
    Vectorcd usi;
    for (int i = 0; i < us.size(); ++i) {
      //Find the region in which ts[i] lies:
      if(tks_index < (tks.size()-1))
        while((ts[i] - tks[tks_index+1])>1e-17)
          tks_index +=1;
        //tks_index = (ts[i] - tks[tks_index+1])<-1e-17?tks_index:(tks_index+1);
      //Find the Basis for the normalized coordinate:
      double u = (ts[i] - tks[tks_index])/(tks[tks_index+1] - tks[tks_index]);
      FindBasis(basis, u);
      us[i].setZero();
      //Weight the spline using the basis:
      for(int degree_count = 0; degree_count <= degree; degree_count++)
      {
        //us[i] += basis[degree_count]*s.segment((tks_index+degree_count)*(this->sys.U.n), this->sys.U.n);
        //cout<<"degree_count: "<<degree_count<<endl;
        //cout<<"s segment: "<<(tks_index+degree_count)*(this->sys.U.n)<<endl;
        us[i] += basis[degree_count]*s.segment((tks_index+degree_count)*(this->sys.U.n),(this->sys.U.n));
      }

    //  cout<<"ts["<<i<<"]: "<<ts[i]<<endl;
      //cout<<"us["<<i<<"]: "<<us[i].transpose()<<"ts: "<<ts[i]<<endl;
      this->sys.Step(xs[i+1], us[i], ts[i+1] - ts[i], p);
      //cout<<"Xs["<<(i+1)<<"]"<<xs[i+1].transpose()<<endl;//#DEBUG
      //cout<<"us["<<i<<"]"<<us[i].transpose()<<endl;
    }
    //getchar();
    return true;
  }
}

#endif
