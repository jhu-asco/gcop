#include "model.h"
#include <assert.h>
#include <iostream>

using namespace Eigen;
using namespace gcop;
using namespace std;

Model::Model(int nx, int nz, int nu, int nc) :
  nx(nx), nz(nz), nu(nu), nc(nc), nr(nx - nc), quat(false),
  Q(nr, nr),
  R(nz, nz)  
{
  //  Q << MatrixXf::Identity(nr, nr);
  //  R << MatrixXf::Identity(nz, nz);
}


Model::~Model()
{
}

bool Model::f(VectorXd &xn,
              const VectorXd &x)
{
  cout << "[W] Model::f: unimplmented. xn=" << xn << " x=" << x << endl;
  return false;
}


bool Model::f(VectorXd &xn,
              const VectorXd &x,
              const VectorXd &u)
{
  cout << "[W] Model::f: unimplmented. xn=" << xn << " x=" << x << " u=" << u << endl;
  return false;
}


bool Model::h(VectorXd &z,
              const VectorXd &x)
{
  cout << "[W] Model::h: unimplmented. z=" << z << " x=" << x << endl;
  return false;
}


bool Model::dfdx(MatrixXd &A,
                 const VectorXd &x)
{
  cout << "[W] Model::dfdx: unimplmented. A=" << A << " x=" << x << endl;
  return false;
}

bool Model::dfdx(MatrixXd &A,
                 const VectorXd &x,
                 const VectorXd &u)
{
  cout << "[W] Model::dfdx: unimplmented. A=" << A << " x=" << x << " u=" << u << endl;
  return false;

}

bool Model::dfdu(MatrixXd &B,
                 const VectorXd &x,
                 const VectorXd &u)
{
  cout << "[W] Model::dfdu: unimplmented. B=" << B << " x=" << x << " u=" << u << endl;
  return false;

}

bool Model::dhdx(MatrixXd &H, 
                 const VectorXd &x)
{
  cout << "[W] Model::dhdx: unimplmented. H=" << H << " x=" << x << endl;
  return false;
}
