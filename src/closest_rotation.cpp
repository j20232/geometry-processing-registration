#include "closest_rotation.h"
#include <Eigen/Dense>
#include <Eigen/LU>

void closest_rotation(const Eigen::Matrix3d& M, Eigen::Matrix3d& R) {
  Eigen::JacobiSVD<Eigen::MatrixXd> SVD(
      M, Eigen::ComputeThinU | Eigen::ComputeThinV);
  Eigen::Matrix3d U = SVD.matrixU();
  Eigen::Matrix3d Vt = (SVD.matrixV()).transpose();
  Eigen::Matrix3d Omega = Eigen::Matrix3d::Identity();
  Omega(2, 2) = (U * Vt).determinant();
  R = U * Omega * Vt;
}
