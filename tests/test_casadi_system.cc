#include "gcop_conversions.h"
#include <Eigen/Dense>
#include <casadi/casadi.hpp>
#include <gtest/gtest.h>

using namespace gcop;

TEST(EigenConversionToDM, Vector3d) {
  Eigen::Vector3d input = Eigen::Vector3d::Zero();
  input(1) = 1.0;
  casadi::DM output = conversions::convertEigenToDM(input);
  std::vector<double> elems = output.get_elements();
  // verify output is same as input
  for (int i = 0; i < 3; ++i) {
    ASSERT_EQ(elems.at(i), input(i));
  }
}

TEST(EigenConversionToDM, VectorXd) {
  Eigen::VectorXd input(20);
  input.setRandom();
  casadi::DM output = conversions::convertEigenToDM(input);
  std::vector<double> elems = output.get_elements();
  // verify output is same as input
  for (int i = 0; i < 3; ++i) {
    ASSERT_EQ(elems.at(i), input(i));
  }
}

TEST(EigenConversionToDM, Matrix3d) {
  Eigen::Matrix3d input;
  input.setRandom();
  casadi::DM output = conversions::convertEigenToDM(input);
  std::vector<double> elems = output.get_elements();
  auto row_vec = output.get_row();
  int ncols = 3;
  int nrows = 3;
  for (int col_ind = 0; col_ind < ncols; ++col_ind) {
    for (int row_ind = 0; row_ind < nrows; ++row_ind) {
      int index = row_ind + col_ind * nrows;
      ASSERT_EQ(row_vec[index], row_ind);
      ASSERT_EQ(elems.at(index), input(row_ind, col_ind));
    }
  }
}

TEST(EigenConversionToDM, MatrixXd) {
  int nrows = 12;
  int ncols = 10;
  Eigen::MatrixXd input(nrows, ncols);
  input.setRandom();
  casadi::DM output = conversions::convertEigenToDM(input);
  std::vector<double> elems = output.get_elements();
  auto row_vec = output.get_row();
  for (int col_ind = 0; col_ind < ncols; ++col_ind) {
    for (int row_ind = 0; row_ind < nrows; ++row_ind) {
      int index = row_ind + col_ind * nrows;
      ASSERT_EQ(row_vec[index], row_ind);
      ASSERT_EQ(elems.at(index), input(row_ind, col_ind));
    }
  }
}

TEST(DMConversionToEigen, Vector3d) {
  casadi::DM input = casadi::DM::rand(3, 1);
  Eigen::Vector3d output = conversions::convertDMToEigen(input);
  std::vector<double> elems = input.get_elements();
  // verify output is same as input
  for (int i = 0; i < 3; ++i) {
    ASSERT_EQ(elems.at(i), output(i));
  }
}

TEST(DMConversionToEigen, MatrixXd) {
  int nrows = 12;
  int ncols = 10;
  casadi::DM input = casadi::DM::rand(nrows, ncols);
  Eigen::MatrixXd output = conversions::convertDMToEigen(input);
  std::vector<double> elems = input.get_elements();
  auto row_vec = input.get_row();
  for (int col_ind = 0; col_ind < ncols; ++col_ind) {
    for (int row_ind = 0; row_ind < nrows; ++row_ind) {
      int index = row_ind + col_ind * nrows;
      ASSERT_EQ(row_vec[index], row_ind);
      ASSERT_EQ(elems.at(index), output(row_ind, col_ind));
    }
  }
}

TEST(DMConversionToEigen, SparseMatrixXd) {
  int nrows = 12;
  int ncols = 10;
  casadi::Sparsity sinput(nrows, ncols);
  sinput.add_nz(4, 4);
  sinput.add_nz(4, 9);
  sinput.add_nz(6, 9);
  sinput.add_nz(2, 3);
  casadi::DM input(sinput);
  for (int i = 0; i < 4; ++i) {
    input.ptr()[i] = 5.0;
  }
  Eigen::MatrixXd output = conversions::convertDMToEigen(input);
  std::vector<double> elems = input.get_nonzeros();
  auto row_vec = sinput.get_row();
  auto col_vec = sinput.get_col();
  for (int i = 0; i < 4; ++i) {
    ASSERT_EQ(output(row_vec.at(i), col_vec.at(i)), elems[i]);
  }
}

int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
