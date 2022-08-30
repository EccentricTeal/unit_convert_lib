#ifndef UNIT_CONVERT_LIB__ANGLE__HH
#define UNIT_CONVERT_LIB__ANGLE__HH

#include <cmath>

namespace unitcon::angle
{
  double rad2deg( double rad ){ return( ( rad * 180.0 ) / M_PI ); }
  double deg2rad( double deg ){ return( ( deg * M_PI ) / 180.0 ); }

}

#endif