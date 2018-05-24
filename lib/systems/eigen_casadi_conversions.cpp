#include "eigen_casadi_conversions.h"

using namespace casadi;

/**
 * @brief convertDMToEigen Helper function to convert casadi matrix to Eigen
 * matrix
 * @param in A casadi matrix can also be a vector if ncols = 1
 * @return  Eigen Matrix with same data as casadi matrix
 */
Eigen::MatrixXd gcop::eigen_casadi_conversions::convertDMToEigen(const DM &in) {
  int rows = in.rows();
  int cols = in.columns();
  Eigen::MatrixXd out(rows, cols);
  std::vector<double> elems = in.get_elements();
  std::memcpy(out.data(), elems.data(), sizeof(double) * rows * cols);
  return out;
}

/**
 * @brief convertEigenToDM Helper function to convert Eigen matrix to casadi
 * matrix
 * @param in An eigen matrix
 * @return A casadi matrix
 */
DM gcop::eigen_casadi_conversions::convertEigenToDM(const Eigen::MatrixXd &in) {
  int rows = in.rows();
  int cols = in.cols();
  cs::DM out = cs::DM::zeros(rows, cols);
  std::memcpy(out.ptr(), in.data(), sizeof(double) * rows * cols);
  return out;
}
