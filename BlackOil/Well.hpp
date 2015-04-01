#pragma once

#include <cmath>

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

struct WellInfo
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
      return sqrt( _Kx * _Ky ) / ( 141.2 * log( _re /_rw ) + _s );
   }

   const std::size_t m_i, m_j, m_k;
   const double m_rw;
   const double m_re;
   const double m_WI;
};


class Well_ConstBHP : public WellInfo
{
public:
   Well_ConstBHP ( double _bhp,
		   std::size_t _i, std::size_t _j, std::size_t _k, double _rw,
		   double _dX, double _dY, double _Kx, double _Ky )
      : WellInfo( _i, _j, _k, _rw, _dX, _dY, _Kx, _Ky ),
	m_bhp( _bhp )
   {}

   double BHP () const
   {
      return m_bhp;
   }
  
private:
   const double m_bhp;
   //const WellInfo info;
};

class Well_ConstRate : public WellInfo
{
public:
   Well_ConstRate ( int _tag, double _rate,
		    std::size_t _i, std::size_t _j, std::size_t _k, double _rw,
		    double _dX, double _dY, double _Kx, double _Ky )
      : WellInfo( _i, _j, _k, _rw, _dX, _dY, _Kx, _Ky ),
	tag( _tag ), m_rate( _rate )
   {
      // tag == 0, const oil rate
      // tag == 1, const water rate
      // tag == 2, const total rate
   }

   int Tag ()     const  { return tag; }
   double Rate () const  { return m_rate; }
   
private:
   const int tag;
   const double m_rate;
   //const WellInfo info;
};
