#pragma once
#include <vector>
#include <iostream>

#include "adetl/scalars/ADscalar.hpp"
#include "adetl/systems/ADvector.hpp"

typedef adetl::ADscalar<> ADscalar;
typedef adetl::ADvector   ADvector;

using std::vector;


struct ConstProperty
{
   /** compressibility **///---------------------------------------------------------------
   const ADscalar Co = 1.0E-05;
   const ADscalar Cw = 3.0E-06;

   /** viscosoty **///---------------------------------------------------------------------
   const ADscalar VisO = 1.1;
   const ADscalar VisW = 0.8;

   /** volume factor **///-----------------------------------------------------------------
   const ADscalar Bob = 1.0;
   const ADscalar Bwb = 1.0;

   /** density **///----------------------------------------------------------------------
   ADscalar DensOil_sc = 53;  // [ lbm/ft3 ]
   ADscalar DensWat_sc = 65;  // [ lbm/ft3 ]
};

struct PROPERTY
{
   static ADscalar Pw ( const ADscalar& _Po )
   {
      ADscalar ret( 0.0 );
      ret = _Po;
      return ret;
   }

   static ADscalar Bo ( const ADscalar& _Po,
			const ADscalar& _Bob, const ADscalar& _Co, const ADscalar& _Pb )
   {
      ADscalar ret( 0.0 );
      ret = _Bob * ( 1.0 - _Co * ( _Po - _Pb ) );
      return ret;
   }
   static ADscalar Bw ( const ADscalar& _Pw,
			const ADscalar& _Bwb, const ADscalar& _Cw, const ADscalar& _Pb )
   {
      ADscalar ret( 0.0 );
      ret = _Bwb * ( 1.0 - _Cw * ( _Pw - _Pb ) );
      return ret;
   }

   static ADscalar Density_Oil ( const ADscalar& _Bo,
				 const ADscalar& _DensOil_sc )
   {
      ADscalar ret( 0.0 );
      ret = _DensOil_sc / _Bo;
      return ret;
   }

   static ADscalar Density_Wat ( const ADscalar& _Bw,
				 const ADscalar& _DensWat_sc )
   {
      ADscalar ret( 0.0 );
      ret = _DensWat_sc / _Bw;
      return ret;
   }

   /** Specific Weight **///-----------------------------------------------------------------
   static ADscalar SpecificWeight ( const ADscalar& _density )
   {
      ADscalar ret( 0.0 );
      ret = _density / 144.0;
      return ret;
   }


   // /** Potential **///-----------------------------------------------------------------------
   // static ADscalar Potential ( const ADscalar& _pressure_1, const ADscalar& _pressure_2,
   // 			       const ADscalar& _spec_weight, double _del_depth )
   // {
   //    ADscalar ret( 0.0 );
   //    ret = _pressure_1 - _pressure_2 - _spec_weight * _del_depth;
   //    return ret;
   // }

   // static double Potential ( double _pressure_1, double _pressure_2,
   // 			     double _spec_weight, double _del_depth )
   // {
   //    return ( _pressure_1 - _pressure_2 - _spec_weight * _del_depth );
   // }

};
