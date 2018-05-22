#ifndef LOAD_EIGEN_MATRIX_H
#define LOAD_EIGEN_MATRIX_H
#include <Eigen/Dense>
#include <fstream>
#include <iostream>
#include <string>

namespace gcop {
Eigen::MatrixXd loadEigenMatrix(std::string in_file) {
  std::ifstream ifile(in_file);
  int rows, cols;
  ifile >> rows >> cols;
  if (rows <= 0 || cols <= 0) {
    throw std::runtime_error("Rows or cols <=0: " + std::to_string(rows) +
                             ", " + std::to_string(cols));
  }
  Eigen::MatrixXd result(rows, cols);
  double val = 0;
  char delim;
  int index = 0;
  while ((ifile >> val >> delim)) {
    int col_ind = (index % cols);
    int row_ind = (index - col_ind) / cols;
    result(row_ind, col_ind) = val;
    ++index;
  }
  return result;
}
}

#endif // LOAD_EIGEN_MATRIX_H
