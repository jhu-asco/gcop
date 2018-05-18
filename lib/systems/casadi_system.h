#ifndef GCOP_CASADI_SYSTEM_H
#define GCOP_CASADI_SYSTEM_H

#include "system.h"
#include <Eigen/Dense>
#include <assert.h>
#include <casadi/casadi.hpp>
#include <string>

namespace gcop {

using namespace std;
using namespace Eigen;
namespace cs = casadi;

/**
* A wrapper around casadi function interface to implement step functions using
* auto differentiation
*
* Author: Gowtham Garimella ggarime1[at]jhu[dot]edu
*/
template <typename T = VectorXd, int _nx = Dynamic, int _nu = Dynamic,
          int _np = Dynamic>
class CasadiSystem : public System<T, _nx, _nu, _np> {
private:
  using Base = System<T, _nx, _nu, _np>;
  static int compile_count;

public:
  CasadiSystem(Manifold<T, _nx> &X, typename Base::Vectormd p, int nu = 0,
               int np = 0, bool use_code_generation = false)
      : Base(X, nu, np), default_parameters_(p),
        use_code_generation_(use_code_generation),
        step_function_instantiated_(false) {}

  /**
  * @brief Generate a casadi function that performs a single step and produces
  * gradients
  *
  * @return Casadi function
  */
  virtual cs::Function casadi_step() = 0;

  void instantiateStepFunction() {
    step_function_ = casadi_step();
    if (use_code_generation_) {
      string function_name = step_function_.name();
      step_function_.generate(function_name);
      string compile_command = ("gcc -fPIC -shared -O3 " + function_name +
                                ".c -o " + function_name + ".so");
      int flag = std::system(compile_command.c_str());
      casadi_assert(flag == 0, "Compilation failed");
      step_function_ = cs::external(function_name);
    }
    step_function_instantiated_ = true;
  }

  double Step(T &xb, double t, const T &xa, const typename Base::Vectorcd &u,
              double h, const typename Base::Vectormd *p = 0,
              typename Base::Matrixnd *A = 0, typename Base::Matrixncd *B = 0,
              typename Base::Matrixnmd *C = 0) {

    // Prepare inputs h, t,xa,u,p
    std::vector<cs::DM> args;
    args.push_back(cs::DM(t));
    args.push_back(cs::DM(h));
    args.push_back(convertEigenToDM(xa));
    args.push_back(convertEigenToDM(u));
    if (p == 0) {
      args.push_back(convertEigenToDM(default_parameters_));
    } else {
      args.push_back(convertEigenToDM(*p));
    }
    std::vector<cs::DM> result = step_function_(args);
    if (result.size() == 0 || result.size() > 4) {
      throw std::runtime_error(
          "The output of the casadi function should be between 1 and 4");
    }
    // Extract results, xb,
    xb = convertDMToEigen(result.at(0));

    if (result.size() > 1) {
      if (A != 0) {
        (*A) = convertDMToEigen(result.at(1));
      }
    }
    if (result.size() > 2) {
      if (B != 0) {
        (*B) = convertDMToEigen(result.at(2));
      }
    }
    if (result.size() > 3) {
      if (C != 0) {
        (*C) = convertDMToEigen(result.at(3));
      }
    }
  }

  static Eigen::MatrixXd convertDMToEigen(const cs::DM &in) {
    int rows = in.rows();
    int cols = in.columns();
    Eigen::MatrixXd out(rows, cols);
    std::vector<double> elems = in.get_elements();
    std::memcpy(out.data(), elems.data(), sizeof(double) * rows * cols);
    return out;
  }

  static cs::DM convertEigenToDM(const Eigen::MatrixXd &in) {
    int rows = in.rows();
    int cols = in.cols();
    cs::DM out = cs::DM::zeros(rows, cols);
    std::memcpy(out.ptr(), in.data(), sizeof(double) * rows * cols);
    return out;
  }

private:
  typename Base::Vectormd default_parameters_;
  bool step_function_instantiated_;
  cs::Function step_function_;
  bool use_code_generation_;
};
}

#endif
