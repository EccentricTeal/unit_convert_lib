#include "unit_convert_lib/coordinate.hh"


namespace unitcon::coordinate
{
  /* This calculation is based on next theory:
   * https://www.gsi.go.jp/common/000061216.pdf */

  //Eigen::Vector(0) = X[m] or LAT[deg],  Eigen::Vector(1) = Y[m] or LON[deg]
  Eigen::Vector2d global2xy( Eigen::Vector2d gnsspos, Eigen::Vector2d origin )
  {
    /* Calculate Coefficients */
    const double a = 6378137.0; //[m] : long radius of an ellipsoid.
    const double F = 298.257222101; //rate of inverse oblateness.
    const double m_0 = 0.9999;//Map Correction factor on X-Axis
    const double n = ( 1.0 / ( 2.0 * F - 1.0 ) );

    std::array<double, 6> alpha;
    alpha[0] = 0.0;
    alpha[1] = 
        0.5 * n  -
      ( 2.0 / 3.0 ) * std::pow( n, 2.0 ) + 
      ( 5.0 / 16.0 ) * std::pow( n, 3.0 ) +
      ( 41.0 / 180.0 ) * std::pow( n, 4.0 ) -
      ( 127.0 / 288.0 ) * std::pow( n, 5.0 ) ;
    alpha[2] = 
      ( 13.0 / 48.0 ) * std::pow( n, 2.0 ) -
      ( 3.0 / 5.0 ) * std::pow( n, 3.0 ) +
      ( 557.0 / 1440.0 ) * std::pow( n, 4.0 ) +
      ( 281.0 / 630.0 ) * std::pow( n, 5.0 ) ;
    alpha[3] = 
      ( 61.0 / 240.0 ) * std::pow( n, 3.0 ) -
      ( 103.0 / 140.0 ) * std::pow( n, 4.0 ) +
      ( 15061.0 / 26880.0 ) * std::pow( n, 5.0 ) ;
    alpha[4] = 
      ( 49561.0 / 161280.0 ) * std::pow( n, 4.0 ) -
      ( 179.0 / 168.0 ) * std::pow( n, 5.0 ) ;
    alpha[5] = 
      ( 34729.0 / 80640.0 ) * std::pow( n, 5.0 );

    std::array<double, 6> A;
    A[0] =
      1.0 +
      ( 1.0 / 4.0 ) * std::pow( n, 2.0 ) +
      ( 1.0 / 64.0 ) * std::pow( n, 4.0 );
    A[1] =
      -( 3.0 / 2.0 ) * n +
      ( 3.0 / 16.0 ) * std::pow( n, 3.0 ) +
      ( 3.0 / 128.0 ) * std::pow( n, 5.0 );
    A[2] =
      ( 15.0 / 16.0 ) * std::pow( n, 2.0 ) -
      ( 15.0 / 64.0 ) * std::pow( n, 4.0 );
    A[3] =
      -( 35.0 / 48.0 ) * std::pow( n, 3.0 ) +
      ( 175.0 / 630.0 ) * std::pow( n, 5.0 );
    A[4] =
      ( 315.0 / 512.0 ) * std::pow( n, 4.0 );
    A[5] =
      -( 693.0 / 1280.0 ) * std::pow( n, 5.0 );
    double A_bar = ( m_0 * a / ( 1.0 + n ) ) * A[0];

    origin(0) = unitcon::angle::deg2rad( origin(0) );
    origin(1) = unitcon::angle::deg2rad( origin(1) );
    gnsspos(0) = unitcon::angle::deg2rad( gnsspos(0) );
    gnsspos(1) = unitcon::angle::deg2rad( gnsspos(1) );

    double t =
      std::sinh(
        std::atanh ( std::sin( gnsspos(0) ) ) -
        ( 2.0 * std::sqrt(n) / ( 1.0 + n ) ) * atanh( ( 2.0 * std::sqrt(n) / ( 1.0 + n ) ) * std::sin( gnsspos(0) ) )
      );
    double t_bar = std::sqrt( 1 + std::pow( t, 2.0 ) );
    double lambda_c = std::cos( gnsspos(0) - origin(0) );
    double lambda_s = std::sin( gnsspos(0) - origin(0) );
    double xi = std::atan( t / lambda_c );
    double eta = std::atanh( lambda_s / t_bar );

    double S_bar = A[0] * origin(0);
    for( int i = 1; i <= 5; i++ )
    {
      S_bar += A[i] * std::sin( 2.0 * static_cast<double>(i) * origin(0) );
    }
    S_bar *= m_0 * a / ( 1.0 + n );
    
    /* Calculate Point */
    double x = xi;
    for( int i = 1; i <= 5; i++ )
    {
      x += alpha[i] * std::sin( 2.0 * static_cast<double>(i) * xi ) * std::cosh( 2.0 * static_cast<double>(i) * eta );
    }
    x = A_bar * x - S_bar;

    double y = eta;
    for( int i = 1; i <= 5; i++ )
    {
      y += alpha[i] * std::cos( 2.0 * static_cast<double>(i) * xi ) * std::sinh( 2.0 * static_cast<double>(i) * eta );
    }
    y *= A_bar;

    Eigen::Vector2d localpos = { x, y };
    return localpos;
  }


  Eigen::Vector2d xy2global( Eigen::Vector2d localpos, Eigen::Vector2d origin )
  {
    ;//TODO
  }

}