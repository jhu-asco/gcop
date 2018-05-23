#ifndef LOAD_EIGEN_MATRIX_H
#define LOAD_EIGEN_MATRIX_H
#include <Eigen/Dense>
#include <fstream>
#include <iostream>
#include <string>

namespace gcop {
Eigen::MatrixXd loadEigenMatrix(std::string in_file_path) {
  std::ifstream ifile(in_file_path);
  if (!ifile.good())
    throw std::runtime_error("File not found: " + in_file_path);
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
    if (row_ind >= rows || col_ind >= cols) {
      throw std::runtime_error("Size of matrix wrong: " + std::to_string(rows) +
                               ", " + std::to_string(cols) + ", " +
                               std::to_string(index));
    }
    result(row_ind, col_ind) = val;
    ++index;
  }
  if (index != rows * cols) {
    throw std::runtime_error("Size of matrix wrong: " + std::to_string(rows) +
                             ", " + std::to_string(cols) + ", " +
                             std::to_string(index));
  }
  return result;
}
}

#endif // LOAD_EIGEN_MATRIX_H
