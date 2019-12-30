#include "point_to_point_rigid_matching.h"
#include <igl/polar_svd.h>
#include "closest_rotation.h"

void point_to_point_rigid_matching(const Eigen::MatrixXd& X,
                                   const Eigen::MatrixXd& P, Eigen::Matrix3d& R,
                                   Eigen::RowVector3d& t) {
  // [Inputs] X: source points, P: target points
  // [Outputs] R: rotation matrix, t: translation
  // Closed-form point-to-point minimizer
  Eigen::MatrixXd x_cor = X.rowwise() - X.colwise().mean();
  Eigen::MatrixXd p_cor = P.rowwise() - P.colwise().mean();
  Eigen::MatrixXd M = (p_cor.transpose() * x_cor).transpose();
  closest_rotation(M, R);
  t = (P.colwise().mean()).transpose() - R * ((X.colwise().mean()).transpose());
}
