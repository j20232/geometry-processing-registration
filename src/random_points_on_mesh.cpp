#include "random_points_on_mesh.h"
#include <igl/cumsum.h>
#include <igl/doublearea.h>
#include <random>

void random_points_on_mesh(const int n, const Eigen::MatrixXd& V,
                           const Eigen::MatrixXi& F, Eigen::MatrixXd& X) {
  // Sample n points from (V, F) -> X
  X.resize(n, 3);

  Eigen::VectorXd Ax, cumarea;
  // https://github.com/libigl/libigl/blob/master/include/igl/doublearea.h
  igl::doublearea(V, F, Ax);          // Compute twice triangle areas
  igl::cumsum(Ax * 0.5, 1, cumarea);  // Compute cumulative sum

  std::random_device gen;
  std::uniform_real_distribution<> distribution(0.0, 1.0);

  for (int i = 0; i < n; i++) {
    double ra = distribution(gen) * cumarea(cumarea.rows() - 1);
    int tri = 0;
    while (cumarea(tri) < ra) tri++;
    Eigen::Vector3d v[] = {V.row(F(tri, 0)), V.row(F(tri, 1)),
                           V.row(F(tri, 2))};
    double alpha = distribution(gen);
    double beta = distribution(gen);

    // outside case
    if (alpha + beta > 1) {
      alpha = 1 - beta;
      beta = 1 - alpha;
    }

    X.row(i) = v[0] + alpha * (v[1] - v[0]) + beta * (v[2] - v[0]);
  }
}
