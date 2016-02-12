#pragma once
#include <vector>
#include <iostream>

#include "adetl/scalars/ADscalar.hpp"
#include "adetl/systems/ADvector.hpp"

typedef adetl::ADscalar<> ADs;
typedef adetl::ADvector   ADv;

using std::vector;


struct ConstFluidProperty
{
   /** reference pressure **///------------------------------------------------------------
   const ADs Pb = 14.7; // psi

   /** compressibility **///---------------------------------------------------------------
   const ADs Co = 1.0E-05;
   const ADs Cw = 3.0E-06;

   /** viscosoty **///---------------------------------------------------------------------
   const ADs VisO = 1.1;
   const ADs VisW = 0.8;

   /** volume factor **///-----------------------------------------------------------------
   const ADs Bob = 1.0;
   const ADs Bwb = 1.0;

   /** density **///----------------------------------------------------------------------
   ADs DensOil_sc = 53;  // [ lbm/ft3 ]
   ADs DensWat_sc = 65;  // [ lbm/ft3 ]
};

struct FluidProperty
{
   static ADs Bo( const ADs& Po,
		  const ADs& Bob, const ADs& Co, const ADs& Pb )
   {
      ADs ret( 0.0 );
      ret = Bob * ( 1.0 - Co * ( Po - Pb ) );
      return ret;
   }

   static ADs Bw( const ADs& Pw,
		  const ADs& Bwb, const ADs& Cw, const ADs& Pb )
   {
      ADs ret( 0.0 );
      ret = Bwb * ( 1.0 - Cw * ( Pw - Pb ) );
      return ret;
   }

   static ADs Density_Oil( const ADs& Bo,
			   const ADs& DensOil_sc )
   {
      ADs ret( 0.0 );
      ret = DensOil_sc / Bo;
      return ret;
   }

   static ADs Density_Wat ( const ADs& Bw,
			    const ADs& DensWat_sc )
   {
      ADs ret( 0.0 );
      ret = DensWat_sc / Bw;
      return ret;
   }

   /** Two-Phase Mixture Density **///-------------------------------------------------------
   static ADs Density_Mix( const ADs& dens_oil, const ADs& dens_wat, const ADs& oil_frac )
   {
      ADs ret( 0.0 );
      ret = dens_oil * oil_frac + dens_wat * ( 1.0 - oil_frac );
      return ret;
   }

   /** Specific Weight **///-----------------------------------------------------------------
   static ADs SpecificWeight( const ADs& density )
   {
      ADs ret( 0.0 );
      ret = density / 144.0;
      return ret;
   }

   // /** Potential **///-----------------------------------------------------------------------
   // static ADs Potential( const ADs& _pressure_1, const ADs& _pressure_2,
   // 			       const ADs& _spec_weight, double _del_depth )
   // {
   //    ADs ret( 0.0 );
   //    ret = _pressure_1 - _pressure_2 - _spec_weight * _del_depth;
   //    return ret;
   // }

   // static double Potential( double _pressure_1, double _pressure_2,
   // 			     double _spec_weight, double _del_depth )
   // {
   //    return ( _pressure_1 - _pressure_2 - _spec_weight * _del_depth );
   // }

};
