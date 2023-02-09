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

  class GlobalLocalConvert
  {
    /* Coefficients */
    private:
      const double a = 6378137.0; //[m] : long radius of an ellipsoid.
      const double F = 298.257222101; //rate of inverse oblateness.
      const double m_0 = 0.9999;//Map Correction factor on X-Axis
      const double n = ( 1.0 / ( 2.0 * F - 1.0 ) );
      std::array<double, 6> alpha;
      std::array<double, 6> beta;
      std::array<double, 7> theta;
      std::array<double, 6> A;
      double A_bar;

    /* Constructor */
    public:
      GlobalLocalConvert();

    /* Private Method */
    private:
      void init_alpha( void );
      void init_beta( void );
      void init_theta( void );
      void init_A( void );
      void init_Abar( void );
      double calc_Sbar( double origin_lat_rad );
      double calc_t( double lat_rad );
      double calc_epsilon_dash( double epsilon, double eta );
      double calc_eta_dash( double epsilon, double eta );
      
    /* Public Method */
    //Eigen::Vector(0) = X[m] or LAT[deg],  Eigen::Vector(1) = Y[m] or LON[deg]
    public:
      Eigen::Vector2d global2xy(
        double lat_deg, double lon_deg,
        double origin_lat_deg, double origin_lon_deg
      );
      Eigen::Vector2d xy2global(
        double local_x, double local_y,
        double origin_lat_deg, double origin_lon_deg
      );
  };
  
  

}

#endif