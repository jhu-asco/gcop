#ifndef CASADI_EIGEN_CONVERSIONS_H
#define CASADI_EIGEN_CONVERSIONS_H
#include <Eigen/Dense>
#include <casadi/casadi.hpp>

namespace gcop {
namespace eigen_casadi_conversions {
namespace cs = casadi;
/**
 * @brief convertDMToEigen Helper function to convert casadi matrix to Eigen
 * matrix
 * @param in A casadi matrix can also be a vector if ncols = 1
 * @return  Eigen Matrix with same data as casadi matrix
 */
Eigen::MatrixXd convertDMToEigen(const cs::DM &in);

/**
 * @brief convertEigenToDM Helper function to convert Eigen matrix to casadi
 * matrix
 * @param in An eigen matrix
 * @return A casadi matrix
 */
cs::DM convertEigenToDM(const Eigen::MatrixXd &in);
}
}

#endif // CASADI_EIGEN_CONVERSIONS_H
