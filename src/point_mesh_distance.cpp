#include "point_mesh_distance.h"
#include <igl/per_face_normals.h>
#include "point_triangle_distance.h"

void point_mesh_distance(const Eigen::MatrixXd& X, const Eigen::MatrixXd& VY,
                         const Eigen::MatrixXi& FY, Eigen::VectorXd& D,
                         Eigen::MatrixXd& P, Eigen::MatrixXd& N) {
  // Compute distances D
  // between a set of given points X and their closest points P
  // on a given mesh with vertex positions VY and face indices FY.
  P.resizeLike(X);
  N = Eigen::MatrixXd::Zero(X.rows(), X.cols());
  D.resize(X.rows());
  Eigen::MatrixXd normals;
  igl::per_face_normals(VY, FY, Eigen::Vector3d(1, 1, 1).normalized(), normals);
  for (int i = 0; i < X.rows(); i++) {
    double min_distance = 10000000, d0;
    int tri_index;
    Eigen::RowVector3d min_p, p0;
    for (int j = 0; j < FY.rows(); j++) {
      point_triangle_distance(X.row(i), VY.row(FY(j, 0)), VY.row(FY(j, 1)),
                              VY.row(FY(j, 2)), d0, p0);
      if (d0 < min_distance) {
        min_distance = d0;
        min_p = p0;
        tri_index = j;
      }
    }
    D(i) = min_distance;
    P.row(i) = min_p;
    N.row(i) = normals.row(tri_index);
  }
}
