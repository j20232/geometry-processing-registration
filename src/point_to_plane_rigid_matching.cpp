#include "point_to_plane_rigid_matching.h"
#include <Eigen/Dense>
#include "closest_rotation.h"

void point_to_plane_rigid_matching(const Eigen::MatrixXd& X,
                                   const Eigen::MatrixXd& P,
                                   const Eigen::MatrixXd& N, Eigen::Matrix3d& R,
                                   Eigen::RowVector3d& t) {
  // [Inputs] X: source points, P: target points, N: target normals
  // [Outputs] R: rotation matrix, t: translation
  // (50) ~ (52)
  typedef Eigen::MatrixXd Mat;
  typedef Eigen::VectorXd Vec;
  int k = X.rows();

  // (51)
  Mat A(3 * k, 6);
  A << Vec::Zero(k), X.col(2), -X.col(1),        // (0, 0:3)
      Vec::Ones(k), Vec::Zero(k), Vec::Zero(k),  // (0, 3:6)
      -X.col(2), Vec::Zero(k), X.col(0),         // (1, 0:3)
      Vec::Zero(k), Vec::Ones(k), Vec::Zero(k),  // (1, 3:6)
      X.col(1), -X.col(0), Vec::Zero(k),         // (2, 0:3)
      Vec::Zero(k), Vec::Zero(k), Vec::Ones(k);  // (2, 3:6)

  Vec B(3 * k);
  B << X.col(0) - P.col(0), X.col(1) - P.col(1), X.col(2) - P.col(2);

  Mat Nd(k, 3 * k);
  Mat Nd0[] = {N.col(0).asDiagonal(), N.col(1).asDiagonal(),
               N.col(2).asDiagonal()};
  Nd << Nd0[0], Nd0[1], Nd0[2];

  A = Nd * A;
  B = Nd * B;

  // differentiate (52)
  Vec u = (A.transpose() * A).inverse() * (-A.transpose() * B);

  // (21)
  Eigen::Matrix3d M(3, 3);
  M << 1, u(2), -u(1), -u(2), 1, u(0), u(1), -u(0), 1;

  closest_rotation(M, R);
  t = Eigen::Vector3d(u(3), u(4), u(5));
}
