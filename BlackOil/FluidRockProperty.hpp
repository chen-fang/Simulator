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

struct ReadFunction
{
   static void Read_Layer_Property ( std::string _filename,
				     std::vector<double>& _porosity_b,
				     std::vector<double>& _Kx,
				     std::vector<double>& _Ky,
				     std::vector<double>& _Kz )
   {
      /* File Format
       * For each row: poro_b -- Kx -- Ky -- Kz
       */
      std::ifstream input( _filename );
      std::string tmp_line;
      
      getline( input, tmp_line ); // skip
      getline( input, tmp_line ); // skip
      
      while( getline( input, tmp_line ) )
      {
	 std::istringstream ss( tmp_line );
	 double value;
	 ss >> value;
	 _porosity_b.push_back( value );
	 ss >> value;
	 _Kx.push_back( value );
	 ss >> value;
	 _Ky.push_back( value );
	 ss >> value;
	 _Kz.push_back( value );
      }
   }
};

struct ConstProperty
{
   /** saturation **///---------------------------------------------------------------------
   const ADscalar Siw = 0.3;
   const ADscalar Sor = 0.2;
   
   /** pressure **///-----------------------------------------------------------------------
   const ADscalar Pb = 14.7;  // [ psi ]
   const ADscalar Poi = 3500; // [ psi ], initial oil phase pressure @ bottom of reservoir
   const ADscalar Pc = 0.0;   // [ psi ]

   /** compressibility **///---------------------------------------------------------------
   const ADscalar Co = 1.0E-05;
   const ADscalar Cw = 3.0E-06;
   const ADscalar Cr = 0.0E-00;

   /** viscosoty **///---------------------------------------------------------------------
   const ADscalar VisO = 1.1;
   const ADscalar VisW = 0.8;

   /** volume factor **///-----------------------------------------------------------------
   const ADscalar Bob = 1.0;
   const ADscalar Bwb = 1.0;

   /** density **///----------------------------------------------------------------------
   ADscalar DensOil_sc = 53;  // [ lbm/ft3 ]
   ADscalar DensWat_sc = 65;  // [ lbm/ft3 ]

   // The followings will be read from files
   void Read_From_File( std::string _filename )
   {
      ReadFunction :: Read_Layer_Property( _filename, vec_porosity_b, vec_Kx, vec_Ky, vec_Kz );
   }
   
   /** porosity **///---------------------------------------------------------------------
   vector<double> vec_porosity_b; // read from file
   
   /** permeability **///-----------------------------------------------------------------
   vector<double> vec_Kx;
   vector<double> vec_Ky;
   vector<double> vec_Kz;
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

   static ADscalar Swd ( const ADscalar& _Sw,
			 const ADscalar& _Siw, const ADscalar& _Sor )
   {
      ADscalar ret( 0.0 );
      ret =  ( _Sw - _Siw ) / ( 1.0 - _Sor - _Siw );
      return ret;
   }

   static ADscalar Krw ( const ADscalar& _Swd )
   {
      double a1 = 0.2;
      double m1 = 2.5;
      ADscalar ret( 0.0 );
      ret = a1 * pow( _Swd, m1 );
      return ret;
   }

   static ADscalar Kro ( const ADscalar& _Swd )
   {
      double a2 = 1.0;
      double m2 = 2.5;
      ADscalar ret( 0.0 );
      ret = a2 * pow( 1.0-_Swd, m2 );
      return ret;
   }

   static double K_Interface ( double _K_1, double _K_2 )
   {
      return 2.0 * _K_1 * _K_2 / ( _K_1 + _K_2 );
   }

   static ADscalar Porosity ( const ADscalar& _Po, const double _Porosity_b,
			      const ADscalar& _Cr, const ADscalar& _Pb )
   {
      ADscalar ret( 0.0 );
      ret = _Porosity_b * ( 1.0 + _Cr * ( _Po - _Pb ) );
      return ret;
   }

   /** Transmissibility **///--------------------------------------------------------------
   static double Trans_Constant ( double _del_1, double _del_2, double _del_3 )
   {
      return 1.127E-03 * _del_1 * _del_2 / _del_3;
   }

   static ADscalar Transmissibility ( double _K, double _trans_constant, const ADscalar& _Kr,
				      const ADscalar& _viscosity, const ADscalar& _volume_factor )
   {
      ADscalar ret( 0.0 );
      ret = _trans_constant * _K * _Kr / _viscosity / _volume_factor;
      return ret;
   }

   /** Specific Weight **///-----------------------------------------------------------------
   static ADscalar SpecificWeight ( const ADscalar& _density )
   {
      ADscalar ret( 0.0 );
      ret = _density / 144.0;
      return ret;
   }

   /** Potential **///-----------------------------------------------------------------------
   static ADscalar Potential ( const ADscalar& _pressure_1, const ADscalar& _pressure_2,
			       const ADscalar& _spec_weight, double _del_depth )
   {
      ADscalar ret( 0.0 );
      ret = _pressure_1 - _pressure_2 - _spec_weight * _del_depth;
      return ret;
   }
};
