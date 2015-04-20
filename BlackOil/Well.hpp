#pragma once

#include <cmath>
#include "adetl/scalars/ADscalar.hpp"
#include "adetl/systems/ADvector.hpp"

typedef adetl::ADscalar<> ADscalar;
typedef adetl::ADvector   ADvector;

// struct WellFunction
// {
//    static double Compute_re ( double _dX, double _dY,
// 			      double _Kx, double _Ky,
// 			      double _rw)
//    {
//       double r = _Kx / _Ky;
//       double t = _dY / _dX;
//       return 0.28 * _dX * sqrt( 1.0 + r * t * t ) / ( 1.0 + sqrt(r) );
//    }

//    static double Compute_WI ( double _Kx, double _Ky,
// 			      double _rw, double _re, double _s = 0.0 )
//    {
//       return sqrt( _Kx * _Ky ) / ( 141.2 * log( _re /_rw ) + _s );
//    }
// };

class WellInfo
{
public:
   WellInfo ( std::size_t _i, std::size_t _j, std::size_t _k, double _rw,
	      double _dX, double _dY, double _Kx, double _Ky )
      : m_i( _i ), m_j( _j ), m_k( _k ), m_rw( _rw ),
	m_re( Compute_re( _dX, _dY, _Kx, _Ky, _rw ) ),
	m_WI( Compute_WI( _Kx, _Ky, _rw, m_re ) )
   {}

   std::size_t I ()  const { return m_i; }
   std::size_t J ()  const { return m_j; }
   std::size_t K ()  const { return m_k; }
   double      WI () const { return m_WI; }
   double      Re () const { return m_re; }

protected:
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
      return sqrt( _Kx * _Ky ) / 141.2 / ( log( _re /_rw ) + _s );
   }

   const std::size_t m_i, m_j, m_k;
   const double m_rw;
   const double m_re;
   const double m_WI;
   //bool  active;
};


class Producer : public WellInfo // producer with constant bottom pressure
{
public:
   Producer ( double _bhp,
	      std::size_t _i, std::size_t _j, std::size_t _k, double _rw,
	      double _dX, double _dY, double _Kx, double _Ky )
      : WellInfo( _i, _j, _k, _rw, _dX, _dY, _Kx, _Ky ),
	m_bhp( _bhp )
   {}

   int& Number ()       { return m_number; }
   int  Number () const { return m_number; }
   double BHP  () const { return m_bhp;            }

private:
   int m_number;
   const double m_bhp;
};

class Injector : public WellInfo // injector with constant water rate
{
public:
   Injector ( double _water_rate,
	      std::size_t _i, std::size_t _j, std::size_t _k, double _rw,
	      double _dX, double _dY, double _Kx, double _Ky )
      : WellInfo( _i, _j, _k, _rw, _dX, _dY, _Kx, _Ky ),
	m_water_rate( _water_rate )
   {}

   int& Number ()        { return m_number; }
   int  Number () const  { return m_number; }
   double Rate () const  { return m_water_rate; }
   
private:
   int m_number;
   const double m_water_rate;
};
