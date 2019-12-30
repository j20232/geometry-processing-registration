#include "point_triangle_distance.h"
#include <Eigen/Dense>

void point_triangle_distance(const Eigen::RowVector3d& x,
                             const Eigen::RowVector3d& a,
                             const Eigen::RowVector3d& b,
                             const Eigen::RowVector3d& c, double& d,
                             Eigen::RowVector3d& p) {
  Eigen::RowVector3d normal = (c - a).cross(c - b);
  normal.normalize();

  // project x onto the triangle abc
  Eigen::RowVector3d px = x + (a - x).dot(normal) * normal;

  // Ref: http://www.sousakuba.com/Programming/gs_hittest_point_triangle.html
  bool in_ab = (b - a).cross(px - a).dot(normal) >= 0;
  bool in_bc = (c - b).cross(px - b).dot(normal) >= 0;
  bool in_ca = (a - c).cross(px - c).dot(normal) >= 0;

  if (in_ab && in_bc && in_ca) {
    p = px;  // inside the triangle abc
  } else if (in_bc && in_ca) {
    p = a + (px - a).dot(b - a) * (b - a);
  } else if (in_ab && in_ca) {
    p = b + (px - b).dot(c - b) * (c - b);
  } else if (in_ab && in_bc) {
    p = c + (px - c).dot(a - c) * (a - c);
  } else if (!in_ab && !in_ca) {
    p = a;
  } else if (!in_ab && !in_bc) {
    p = b;
  } else {
    p - c;
  }

  d = (x - p).norm();
}
