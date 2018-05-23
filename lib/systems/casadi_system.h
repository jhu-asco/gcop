#ifndef GCOP_CASADI_SYSTEM_H
#define GCOP_CASADI_SYSTEM_H

#include "eigen_casadi_conversions.h"
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
  using Base = System<T, _nx, _nu, _np>; ///< Base class

public:
  /**
 * @brief CasadiSystem constructor
 * @param X  The state manifold
 * @param p  The default parameters for the system
 * @param nu The control size. Only relevant if _nu is Dynamic
 * @param np The parameter vector size. Only relevant if _np is Dynamic
 * @param use_code_generation  If true will compile the casadi function into a
 * shared library on the fly
 */
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

  /**
   * @brief instantiateStepFunction generate the shared library or create the
   * step function
   */
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

  /**
   * @brief Step Perform a single step of the dynamics
   * @param xb  The state at the next step
   * @param t The curren time
   * @param xa The current state of the system
   * @param u  The current control
   * @param h  The time step for integration
   * @param p  The parameters of the system. If not provided will use default
   * parameters
   * @param A The jacobian wrt current state is returned if provided
   * @param B The jacobian wrt control is returned if provided
   * @param C The jacobian wrt parameters p is returned if provided
   * @return
   */
  double Step(T &xb, double t, const T &xa, const typename Base::Vectorcd &u,
              double h, const typename Base::Vectormd *p = 0,
              typename Base::Matrixnd *A = 0, typename Base::Matrixncd *B = 0,
              typename Base::Matrixnmd *C = 0) {

    // Prepare inputs t, h, xa,u,p
    std::vector<cs::DM> args;
    args.push_back(cs::DM(t));
    args.push_back(cs::DM(h));
    args.push_back(eigen_casadi_conversions::convertEigenToDM(xa));
    args.push_back(eigen_casadi_conversions::convertEigenToDM(u));
    if (p == 0) {
      args.push_back(
          eigen_casadi_conversions::convertEigenToDM(default_parameters_));
    } else {
      args.push_back(eigen_casadi_conversions::convertEigenToDM(*p));
    }
    std::vector<cs::DM> result = step_function_(args);
    if (result.size() == 0 || result.size() > 4) {
      throw std::runtime_error(
          "The output of the casadi function should be between 1 and 4");
    }
    // Extract results, xb,
    xb = eigen_casadi_conversions::convertDMToEigen(result.at(0));

    if (result.size() > 1) {
      if (A != 0) {
        (*A) = eigen_casadi_conversions::convertDMToEigen(result.at(1));
      }
    }
    if (result.size() > 2) {
      if (B != 0) {
        (*B) = eigen_casadi_conversions::convertDMToEigen(result.at(2));
      }
    }
    if (result.size() > 3) {
      if (C != 0) {
        (*C) = eigen_casadi_conversions::convertDMToEigen(result.at(3));
      }
    }
  }

private:
  typename Base::Vectormd default_parameters_; ///< Default system parameters
  bool step_function_instantiated_; ///< Flag to check if casadi step function
                                    /// is created
  cs::Function step_function_;      ///< The instantiated step function
  bool use_code_generation_; ///< Flag to specify whether code generation should
                             /// be used
};
}

#endif
