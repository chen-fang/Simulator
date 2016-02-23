#pragma once
#include <vector>
#include <iostream>

#include "adetl/scalars/ADscalar.hpp"
#include "adetl/systems/ADvector.hpp"

typedef adetl::ADscalar<> ADs;
typedef adetl::ADvector   ADv;

using std::vector;

struct FluidProperty
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


   static ADs Bo( const ADs& Po )
   {
      ADs ret( 0.0 );
      ret = Bob * ( 1.0 - Co * ( Po - Pb ) );
      return ret;
   }

   static ADs Bw( const ADs& Pw )
   {
      ADs ret( 0.0 );
      ret = Bwb * ( 1.0 - Cw * ( Pw - Pb ) );
      return ret;
   }

   static ADs Density_Oil( const ADs& Po )
   {
      ADs Bo( 0.0 );
      Bo = Bo( Po );
      ADs ret( 0.0 );
      ret = DensOil_sc / Bo;
      return ret;
   }

   static ADs Density_Oil( const ADs& Bo )
   {
      ADs ret( 0.0 );
      ret = DensOil_sc / Bo;
      return ret;
   }

   static ADs Density_Wat ( const ADs& Pw )
   {
      ADs Bw( 0.0 );
      Bw = Bw( Pw );
      ADs ret( 0.0 );
      ret = DensWat_sc / Bw;
      return ret;
   }

   static ADs Density_Wat ( const ADs& Bw )
   {
      ADs ret( 0.0 );
      ret = DensWat_sc / Bw;
      return ret;
   }

   /** Two-Phase Mixture Density **///-------------------------------------------------------
   static ADs Density_Mix( const ADs& denL, const ADs& denG, const ADs& Hg )
   {
      ADs ret( 0.0 );
      ret = denL * ( 1.0 - Hg ) + denG * Hg;
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
