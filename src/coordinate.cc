#include "unit_convert_lib/coordinate.hh"


namespace unitcon::coordinate
{
  GlobalLocalConvert::GlobalLocalConvert()
  {
    init_A();
    init_Abar();
    init_alpha();
    init_beta();
    init_theta();
  }

  void GlobalLocalConvert::init_alpha( void )
  {
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
  }

  void GlobalLocalConvert::init_beta( void )
  {
    beta[0] = 0.0;
    beta[1] = 
      0.5 * n  -
      ( 2.0 / 3.0 ) * std::pow( n, 2.0 ) + 
      ( 37.0 / 96.0 ) * std::pow( n, 3.0 ) -
      ( 1.0 / 360.0 ) * std::pow( n, 4.0 ) -
      ( 81.0 / 512.0 ) * std::pow( n, 5.0 ) ;
    beta[2] = 
      ( 1.0 / 48.0 ) * std::pow( n, 2.0 ) +
      ( 1.0 / 15.0 ) * std::pow( n, 3.0 ) -
      ( 437.0 / 1440.0 ) * std::pow( n, 4.0 ) +
      ( 46.0 / 105.0 ) * std::pow( n, 5.0 ) ;
    beta[3] = 
      ( 17.0 / 480.0 ) * std::pow( n, 3.0 ) -
      ( 37.0 / 840.0 ) * std::pow( n, 4.0 ) -
      ( 209.0 / 4480.0 ) * std::pow( n, 5.0 ) ;
    beta[4] = 
      ( 4397.0 / 161280.0 ) * std::pow( n, 4.0 ) -
      ( 11.0 / 504.0 ) * std::pow( n, 5.0 ) ;
    beta[5] = 
      ( 4583.0 / 161280.0 ) * std::pow( n, 5.0 );
  }

  void GlobalLocalConvert::init_theta( void )
  {
    theta[0] = 0.0;
    theta[1] = 
      2.0 * n  -
      ( 2.0 / 3.0 ) * std::pow( n, 2.0 ) -
      2.0 * std::pow( n, 3.0 ) +
      ( 116.0 / 45.0 ) * std::pow( n, 4.0 ) +
      ( 26.0 / 45.0 ) * std::pow( n, 5.0 ) -
      ( 2854.0 / 675.0 ) * std::pow( n, 6.0 );
    theta[2] = 
      ( 7.0 / 3.0 ) * std::pow( n, 2.0 ) -
      ( 8.0 / 5.0 ) * std::pow( n, 3.0 ) -
      ( 227.0 / 45.0 ) * std::pow( n, 4.0 ) +
      ( 2704.0 / 315.0 ) * std::pow( n, 5.0 ) +
      ( 2323.0 / 945.0 ) * std::pow( n, 6.0 );
    theta[3] = 
      ( 56.0 / 15.0 ) * std::pow( n, 3.0 ) -
      ( 136.0 / 35.0 ) * std::pow( n, 4.0 ) -
      ( 1262.0 / 105.0 ) * std::pow( n, 5.0 ) +
      ( 73814.0 / 2835.0 ) * std::pow( n, 6.0 );
    theta[4] = 
      ( 4279.0 / 630.0 ) * std::pow( n, 4.0 ) -
      ( 332.0 / 35.0 ) * std::pow( n, 5.0 ) -
      ( 399572.0 / 14175.0 ) * std::pow( n, 6.0 );
    theta[5] = 
      ( 4174.0 / 315.0 ) * std::pow( n, 5.0 ) -
      ( 144838.0 / 6237.0 ) * std::pow( n, 6.0 );
    theta[6] = 
    ( 601676.0 / 22275.0 ) * std::pow( n, 6.0 );
  }

  void GlobalLocalConvert::init_A( void )
  {
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
      ( 175.0 / 768.0 ) * std::pow( n, 5.0 );
    A[4] =
      ( 315.0 / 512.0 ) * std::pow( n, 4.0 );
    A[5] =
      -( 693.0 / 1280.0 ) * std::pow( n, 5.0 );
  }

  void GlobalLocalConvert::init_Abar( void )
  {
    A_bar = ( m_0 * a * A[0] ) / ( 1.0 + n );
  }

  double GlobalLocalConvert::calc_Sbar( double origin_lat_rad )
  {
    double S_bar = A[0] * origin_lat_rad;
    for( int i = 1; i <= 5; i++ )
    {
      S_bar += A[i] * std::sin( 2.0 * static_cast<double>(i) * origin_lat_rad );
    }
    S_bar *= ( ( m_0 * a ) / ( 1.0 + n ) );

    return S_bar;
  }

  double GlobalLocalConvert::calc_t( double lat_rad )
  {
    double t = std::sinh(
        std::atanh ( std::sin( lat_rad ) ) -
        ( 2.0 * std::sqrt(n) / ( 1.0 + n ) ) * atanh( ( 2.0 * std::sqrt(n) / ( 1.0 + n ) ) * std::sin( lat_rad ) )
      );

    return t;
  }

  double GlobalLocalConvert::calc_epsilon_dash( double epsilon, double eta )
  {
    double epsilon_dash = epsilon;
    for( int i = 1; i <= 5; i++ )
    {
      epsilon_dash -= beta[i] * std::sin( 2.0 * static_cast<double>(i) * epsilon ) * std::cosh( 2.0 * static_cast<double>(i) * eta );
    }

    return epsilon_dash;
  }


  double GlobalLocalConvert::calc_eta_dash( double epsilon, double eta )
  {
    double eta_dash = eta;
    for( int i = 1; i <= 5; i++ )
    {
      eta_dash -= beta[i] * std::cos( 2.0 * static_cast<double>(i) * epsilon ) * std::sinh( 2.0 * static_cast<double>(i) * eta );
    }
    
    return eta_dash;
  }


    
  Eigen::Vector2d GlobalLocalConvert::global2xy(
    double lat_deg, double lon_deg,
    double origin_lat_deg, double origin_lon_deg
  )
  {
    /* Coef Preparation */
    Eigen::Vector2d origin = {
      unitcon::angle::deg2rad( origin_lat_deg ),
      unitcon::angle::deg2rad( origin_lon_deg )
    };
    Eigen::Vector2d gnsspos = {
      unitcon::angle::deg2rad( lat_deg ),
      unitcon::angle::deg2rad( lon_deg )
    };
    double t = calc_t( gnsspos(0) );
    double t_bar = std::sqrt( 1 + std::pow( t, 2.0 ) );
    double lambda_c = std::cos( gnsspos(1) - origin(1) );
    double lambda_s = std::sin( gnsspos(1) - origin(1) );
    double xi = std::atan2( t, lambda_c );
    double eta = std::atanh( lambda_s / t_bar );
    double S_bar = calc_Sbar( origin(0) );

    
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


  Eigen::Vector2d GlobalLocalConvert::xy2global(
    double local_x, double local_y,
    double origin_lat_deg, double origin_lon_deg
  )
  {
    /* Coef Preparation */
    Eigen::Vector2d origin = {
      unitcon::angle::deg2rad( origin_lat_deg ),
      unitcon::angle::deg2rad( origin_lon_deg )
    };

    double S_bar = calc_Sbar( origin(0) );
    double epsilon = ( local_x + S_bar ) / A_bar;
    double eta = local_y / A_bar;
    double epsilon_dash = calc_epsilon_dash( epsilon, eta );
    double eta_dash = calc_eta_dash( epsilon, eta );
    double chi = std::asin( std::sin( epsilon_dash ) / std::cosh( eta_dash ) );


    /* calculate point */
    double lat = unitcon::angle::rad2deg( chi );
    for( int i = 1; i <= 6; i++ )
    {
      lat += unitcon::angle::rad2deg( theta[i] * std::sin( 2.0 * static_cast<double>(i) * chi ) );
    }

    double lon =
      origin_lon_deg +
      unitcon::angle::rad2deg( std::atan2( std::sinh( eta_dash ), std::cos( epsilon_dash ) ) );

    Eigen::Vector2d global_pos = { lat, lon };
    return global_pos;
  }

}