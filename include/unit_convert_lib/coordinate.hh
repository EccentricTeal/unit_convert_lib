#ifndef UNIT_CONVERT_LIB__COORDINATE__HH
#define UNIT_CONVERT_LIB__COORDINATE__HH

#include <cmath>
#include <array>
#include <eigen3/Eigen/Dense>
#include <unit_convert_lib/angle.hh>


namespace unitcon::coordinate
{
  /* This calculation is based on next theory:
   * https://www.gsi.go.jp/common/000061216.pdf */

  //Eigen::Vector(0) = X[m] or LAT[deg],  Eigen::Vector(1) = Y[m] or LON[deg]
  Eigen::Vector2d global2xy( Eigen::Vector2d gnsspos, Eigen::Vector2d origin );
  Eigen::Vector2d xy2global( Eigen::Vector2d localpos, Eigen::Vector2d origin );

}

#endif