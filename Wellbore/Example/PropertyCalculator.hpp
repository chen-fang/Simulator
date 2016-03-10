#pragma once
#include <fstream>
#include <vector>
#include "adetl/scalars/ADscalar.hpp"
#include "func_table.hpp"


class PropertyCalculator
{
private:
   typedef adetl::ADscalar< >         scalar;
   typedef fastl::sorted_key<double>  TableKey;
   typedef fastl::interpolatable_range<double,
				       fastl::forward_linear,
				       fastl::centered_linear,
				       fastl::backward_linear>  TableValue;

public:
   enum Phase { O=0, W=1 };
   enum StatusID { OW };
   typedef struct
   {
      // O, W
      StatusID   status;
      scalar S[2];
      scalar P[2];
      scalar Kr[2];
      scalar b[2];
      scalar bmu[2];
      scalar Rho[2];
      scalar phi;
   }                                            CellProps;

public:
   PropertyCalculator() : CR(3.4e-4), PREF_ROCK(2500.0), 
			  RHO_W_S(63.02), RHO_O_S(45.0),
			  PREF_W(3600.0), BW_REF(1.00341), CW(3.0e-4), MUW_REF(0.52341), CVW(1.2e-6)
   {
      read_PVDO_table( );
      read_SWOF_table( );
   };

   template< typename T >
   void calculate( double PHI, const T & _state,  CellProps & p )
   {
      p.status = _state.status;
      p.P[O]   = _state.Po;
      p.S[W]   = _state.Sw;
      p.S[O]   = 1.0 - p.S[W];
      calc_rock_fluid( p );
      calc_water_props( p );
      calc_oil_props( p );
      calc_phi( PHI, p );
   }

private:

   void calc_phi( double PHI, CellProps &p ) const
   {
      const scalar X = CR * ( p.P[O] - PREF_ROCK );
      p.phi = PHI * ( 1.0 + X + 0.5 * X * X );
      if ( p.phi.value() < 0.005 ) p.phi = 0.05;
   }

   void calc_rock_fluid( CellProps &p )
   {
      std::size_t i = SW_key( p.S[W] );
      KRW_val(  i, SW_key, p.S[W], p.Kr[W] );
      KRO_val(  i, SW_key, p.S[W], p.Kr[O] );
      PCOW_val(  i, SW_key,p.S[W], mTmp[0] );
      p.P[W] = p.P[O] - mTmp[0];
   }

   void calc_water_props( CellProps &p )
   { 
      mTmp[0]    = CW*(p.P[W]-PREF_W);
      p.b[W]     = ( 1.0 + mTmp[0] + 0.5 * mTmp[0] * mTmp[0] ) / BW_REF;
      mTmp[0]    = ( CW - CVW ) * ( p.P[W] - PREF_W );
      p.bmu[W]   = ( 1.0 + mTmp[0] + 0.5 * mTmp[0] * mTmp[0] ) / ( BW_REF * MUW_REF );
      p.Rho[W]   = RHO_W_S * p.b[W];
   }  

   void calc_oil_props( CellProps &p )
   { 
	 std::size_t i = PO_key( p.P[O] );
	 bo_val(   i, PO_key, p.P[O], p.b[O] );
	 bmuo_val( i, PO_key, p.P[O], p.bmu[O] );	 
         p.Rho[O] = RHO_O_S * p.b[O];
  }  

   void
   read_PVDO_table( )
   {
      // read PVDO dead oil table
      double tmp1, tmp2, tmp3;
      std::ifstream strm ( "./input/PVDO.dat" );
      while ( strm.good( ) )
      {
	 if (strm >> tmp1) PO_key.push_back( tmp1 );
	 if (strm >> tmp2) bo_val.push_back( 1.0/tmp2 );
	 if (strm >> tmp3) bmuo_val.push_back( 1.0/(tmp2*tmp3) );
      }
      strm.close(); 
   }

   void
   read_SWOF_table( )
   {
      // read SWOF table
      double tmp;
      std::ifstream strm ( "./input/SWOF.dat" );
      while ( strm.good( ) )
      {
	 if (strm >> tmp) SW_key.push_back( tmp );
	 if (strm >> tmp) KRW_val.push_back( tmp );
	 if (strm >> tmp) KRO_val.push_back( tmp );
	 if (strm >> tmp) PCOW_val.push_back( tmp );
      }
      strm.close(); 
   }

private:
   const double CR;
   const double PREF_ROCK;
   const double RHO_W_S;
   const double RHO_O_S;
   const double PREF_W;
   const double BW_REF;
   const double CW;
   const double MUW_REF;
   const double CVW;
   TableKey     PO_key;
   TableValue   bo_val;
   TableValue   bmuo_val;
   TableKey     SW_key;
   TableValue   KRW_val;
   TableValue   KRO_val;
   TableValue   PCOW_val;
   scalar       mTmp[5];
};
