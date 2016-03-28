#pragma once
#include <vector>
#include "adetl/scalars/ADscalar.hpp"

#include "Solve.hpp" // EPRI Drift-Flux

const double surface_tension = 0.072; // N/m
const double diameter = 0.15; // m
const double g = 9.81; // m/s^2
class PropertyCalculator
{
private:
   typedef adetl::ADscalar< >         scalar;

public:
   typedef struct
   {
      // This struct contains ALL properties needed in the model
      scalar P;
      scalar V[2];
      scalar alpha[2]; // [ air | water ]
      scalar rho[2];   // [ air | water ]
      scalar vis[2];   // [ air | water ]
      scalar Re[2];
      scalar CCFL;
      scalar C0j;
      scalar Vgj;
      scalar friction;
   }                                            CellProps;

public:
   PropertyCalculator() 
   {
     // constant parameters asigned as data members
   };

  // T state struct with P, V0, V1, alpha0
   template< typename T >
   void calculate( const T & _state,  CellProps & p )
   {
     p.P        = _state.P;
     p.V[0]     = _state.V0; // air
     p.V[1]     = _state.V1; // water
     p.alpha[0] = _state.alpha0;
     
     // calculate the rest for example
     p.alpha[1] = 1.0 - p.alpha[0];
     calc_rho( p );
     calc_vis( p );
     calc_Re(  p );
     calc_CCFL( p );
     calc_C0j( p );
     calc_Vgj( p );
     calc_wall_friction( p );
   }

private:
   void calc_rho( CellProps& p )
   {
      // air | water
      p.rho[0] = 1.184; // kg/m3
      p.rho[1] = 1000.0; // kg/m3
   }

   void calc_vis( CellProps& p )
   {
      // air | water
      p.vis[0] = 1.81E-05; // Pa.s
      p.vis[1] = 1.00E-03; // Pa.s
      
   }

   void calc_Re( CellProps& p )
   {
      // air | water
      EPRI_DF epri;
      p.Re[0] = epri.ReyNum( p.rho[0], p.V[0], p.vis[0], diameter );
      p.Re[1] = epri.ReyNum( p.rho[1], p.V[1], p.vis[1], diameter );
   }

   void calc_CCFL( CellProps& p )
   {
      mTmp[0] = p.alpha[0] * p.V[0]; // j_gas
      mTmp[1] = 0.0; // Hg_CCFL

      Solve::Solve_CCFL( p.rho[0], p.rho[1], p.vis[0], p.vis[1], diameter, g, surface_tension, mTmp[0],
			 p.CCFL, mTmp[1] );

      if( p.alpha[0] < mTmp[1] )
      {
	 // cross CCFL boundary
	 // correction
      }

   }
     
   void calc_C0j( CellProps& p )
   {
      EPRI_DF epri;
      mTmp[0] = epri.C0( p.rho[0], p.rho[1], p.Re[0], p.Re[1], p.alpha[0] ); // C0
      mTmp[1] = p.alpha[0] * p.V[0] + p.alpha[1] * p.V[1]; // j = jL + jG
      p.C0j = mTmp[0] * mTmp[1];
   }

   void calc_Vgj( CellProps& p )
   {
      EPRI_DF epri;
      mTmp[1] = p.alpha[1] * p.V[1]; // j_liquid
      
      p.Vgj = epri.Vgj( p.rho[0], p.rho[1], p.Re[0], p.Re[1], mTmp[1], p.CCFL, p.alpha[0], surface_tension, diameter, g );
   }

   void calc_wall_friction( CellProps& p )
   {
      mTmp[0] = p.alpha[0] * p.rho[0] + p.alpha[1] * p.rho[1]; // mixture density
      mTmp[1] = p.alpha[0] * p.vis[0] + p.alpha[1] * p.vis[1]; // mixture viscosity
      mTmp[2] = ( p.alpha[0] * p.rho[0] * p.V[0] + p.alpha[1] * p.rho[1] * p.V[1] ) / mTmp[0]; // mixture velocity

      EPRI_DF epri;
      mTmp[3] = epri.ReyNum( mTmp[0], mTmp[2], mTmp[1], diameter ); // mixuture Re

      // mTmp[4] is friction factor
      mTmp[4] = 64.0 / mTmp[3]; // !!! for test only. Remove when below is completed !!!
      
      if( mTmp[2].value() <= 2000.0 )
	 mTmp[4] = 64.0 / mTmp[3];
      else {}
	 // iteration required

      p.friction = mTmp[4] * mTmp[0] * mTmp[2]*mTmp[2] / ( 2.0 * diameter );
   }
/*   
   void calc_rho_oil( CellProps& p ) // oil in liquid phase
   {
      // constant
      const scalar Pb( 14.7 ); // bubble point [psi]
      const scalar Co( 1.0E-05 );  // compressibility
      const scalar Bob( 1.0 ); // volume factor
      const scalar rho_oil_sc( 53.0 ); // lbm/ft3

      mTmp[0] = Bob * ( 1.0 - Co * ( p.P - _Pb ) ); // Bo
      p.rho[0] = rho_oil_sc / mTmp[0];
   }

   
   void cal_rho_water( CellProp &p ) // water in liquid phase
   {
      // constant
      const scalar Pb( 14.7 ); // bubble point [psi]
      const scalar Cw( 3.0E-06 );  // compressibility
      const scalar Bwb( 1.0 ); // volume factor
      const scalar rho_water_sc( 65.0 ); // lbm/ft3
      
      mTmp[0] = _Bwb * ( 1.0 - _Cw * ( _Pw - _Pb ) ); // Bw
      p.rho[1] = rho_water_sc / mTmp[0];
   }
*/

private:
  scalar       mTmp[5]; // pool of intermediate adscalars
};
