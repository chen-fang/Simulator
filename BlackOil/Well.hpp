#pragma once

#include <cmath>

struct WellFunction
{
  static double Compute_re ( double _dX, double _dY,
			     double _Kx, double _Ky,
			     double _rw)
  {
    double r = _Kx / _Ky;
    double t = _dY / _dX;
    return 0.28 * _dX * sqrt( 1.0 + r * t * t ) / ( 1.0 + sqrt(r) );
  }

  static double Compute_WI ( double _Kx, double _Ky,
			     double _rw, double _re, double _s = 0.0 )
  {
    return sqrt( _Kx * _Ky ) / ( 141.2 * log( _re /_rw ) + _s );
  }
};

struct WellInfo
{
  WellInfo ( std::size_t _i, std::size_t _j, std::size_t _k, double _rw,
	     double _dX, double _dY, double _Kx, double _Ky )
    : i( _i ), j( _j ), k( _k ), rw( _rw ),
      re( WellFunction :: Compute_re( _dX, _dY, _Kx, _Ky, _rw ) ),
      WI( WellFunction :: Compute_WI( _Kx, _Ky, _rw, re ) )
  {}
  
  const std::size_t i, j, k;
  const double rw;
  const double re;
  const double WI;
};


class Well_ConstBHP
{
public:
  Well_ConstBHP ( double _bhp,
		  std::size_t _i, std::size_t _j, std::size_t _k, double _rw,
		  double _dX, double _dY, double _Kx, double _Ky )
    : pressure( _pressure ),
      info( _i, _j, _k, _rw, _dX, _dY, _Kx, _Ky )
  {}

  double BHP () const
  {
    return bhp;
  }

  std::size_t I () const  { return info.i; }
  std::size_t J () const  { return info.j; }
  std::size_t K () const  { return info.k; }
  double      WI () const { return info.WI; }
  
private:
  const double bhp;
  const WellInfo info;
};

class Well_ConstRate
{
public:
  Well_ConstRate ( double _rate,
		   std::size_t _i, std::size_t _j, std::size_t _k, double _rw,
		   double _dX, double _dY, double _Kx, double _Ky )
    : rate( _rate ),
      info( _i, _j, _k, _rw, _dX, _dY, _Kx, _Ky )
  {}

  const static int WellType = 1;
  
  double Rate () const
  {
    return rate;
  }
  
private:
  const double rate;
  const WellInfo info;
};
