#pragma once
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <iostream>

#include "adetl/scalars/ADscalar.hpp"
#include "adetl/systems/ADvector.hpp"

typedef adetl::ADscalar<> ADscalar;
typedef adetl::ADvector   ADvector;

using std::vector;

struct ConstProperty
{
   ConstProperty() : Siw( 0.5 ),
		     Poi( 2.5E+04 ),
		     Co( 3.0E-02 ),
		     VisO( 1.1E-03 ), VisW( 0.8E-03 ),
		     porosity( 0.3 ),
		     Kx( 5.0E-12 )
   {}
   
   /** saturation **///---------------------------------------------------------------------
   const ADscalar Siw;
   //const ADscalar Sor = 0.2;
   
   /** pressure **///-----------------------------------------------------------------------
   const ADscalar Poi; // [ Pa ]

   /** compressibility **///---------------------------------------------------------------
   const ADscalar Co; // [ Kg/m3/Pa ]

   /** viscosoty **///---------------------------------------------------------------------
   const ADscalar VisO; // [ Pa.s ]
   const ADscalar VisW; // [ Pa.s ]
   
   /** porosity **///---------------------------------------------------------------------
   double porosity;
   
   /** permeability **///-----------------------------------------------------------------
   // vector<double> vec_Kx;
   // vector<double> vec_Ky;
   // vector<double> vec_Kz;
   const double Kx;
   // const double Ky = 5.0E-12;
   // const double Kz = 5.0E-12;
};

struct PROPERTY
{
   static ADscalar Density_Oil ( const ADscalar& _Co,
				 const ADscalar& _Pressure )
   {
      ADscalar ret( 0.0 );
      ret = _Co * _Pressure;
      return ret;
   }

   static ADscalar Density_Wat ()
   {
      ADscalar ret( 1.0E+03 ); // [ Kg/m3 ]
      return ret;
   }

   static ADscalar Krw ( const ADscalar& _Sw )
   {
      ADscalar ret( 0.0 );
      ret = pow( _Sw, 2.0 );
      return ret;
   }

   static ADscalar Kro ( const ADscalar& _Sw )
   {
      ADscalar ret( 0.0 );
      ret = pow( 1.0-_Sw, 2.0 );
      return ret;
   }

   static double K_Interface ( double _K )
   {
      return _K;
   }

   /** Transmissibility **///--------------------------------------------------------------
   static double Trans_Constant ( double _DY, double _DZ, double _DX )
   {
      return _DY * _DZ / _DX;
   }

   static ADscalar Transmissibility ( double _K, double _trans_constant,
				      const ADscalar& _Kr,
				      const ADscalar& _density,
				      const ADscalar& _viscosity )
   {
      ADscalar ret( 0.0 );
      ret = _trans_constant * _K * _Kr * _density / _viscosity;
      return ret;
   }

   /** Potential **///-----------------------------------------------------------------------
   static ADscalar Potential ( const ADscalar& _pressure_1, const ADscalar& _pressure_2 )

   {
      ADscalar ret( 0.0 );
      ret = _pressure_1 - _pressure_2;
      return ret;
   }
};
