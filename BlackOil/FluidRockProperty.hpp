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

namespace PROPERTY
{
   /** saturation **///---------------------------------------------------------------------
   const ADscalar Siw = 0.3;
   const ADscalar Sor = 0.2;
   
   /** pressure **///-----------------------------------------------------------------------
   const ADscalar Pb = 14.7;  // [ psi ]
   const ADscalar Poi = 3500; // [ psi ], initial oil phase pressure @ bottom of reservoir
   const ADscalar Pc = 0.0;   // [ psi ]
   ADscalar Pw ( const ADscalar& _Po )
   {
      return _Po;
   }

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
   ADscalar Bo ( const ADscalar& _Po,
		 const ADscalar& _Bob = Bob, const ADscalar& _Co = Co, const ADscalar& _Pb = Pb )
   {
      return _Bob * ( 1.0 - _Co * ( _Po - _Pb ) );
   }
   ADscalar Bw ( const ADscalar& _Pw,
		 const ADscalar& _Bwb = Bwb, const ADscalar& _Cw = Cw, const ADscalar& _Pb = Pb )
   {
      return _Bwb * ( 1.0 - _Cw * ( _Pw - _Pb ) );
   }

   /** density **///----------------------------------------------------------------------
   ADscalar DensOil_sc = 53;  // [ lbm/ft3 ]
   ADscalar DensWat_sc = 65;  // [ lbm/ft3 ]
   ADscalar Density_Oil ( const ADscalar& _Bo,
			  const ADscalar& _DensOil_sc = DensOil_sc )
   {
      return _DensOil_sc / _Bo;
   }

   ADscalar Density_Wat ( const ADscalar& _Bw,
			  const ADscalar& _DensWat_sc = DensWat_sc )
   {
      return _DensWat_sc / _Bw;
   }

   /** permeability **///-----------------------------------------------------------------
   struct AbsK
   {
      static vector<double> vec_Kx;
      static vector<double> vec_Ky;
      static vector<double> vec_Kz;
   };
   vector<double> AbsK:: vec_Kx;
   vector<double> AbsK:: vec_Ky;
   vector<double> AbsK:: vec_Kz;
   
   ADscalar Swd ( const ADscalar& _Sw,
		  const ADscalar& _Siw = Siw, const ADscalar& _Sor = Sor )
   {
      return ( _Sw - _Siw ) / ( 1.0 - _Sor - _Siw );
   }

   ADscalar Krw ( const ADscalar& _Swd )
   {
      double a1 = 0.2;
      double m1 = 2.5;
      return a1 * pow( _Swd, m1 );
   }

   ADscalar Kro ( const ADscalar& _Swd )
   {
      double a2 = 1.0;
      double m2 = 2.5;
      return a2 * pow( 1.0-_Swd, m2 );
   }

   double K_Interface ( double _K_1, double _K_2 )
   {
      return 2.0 * _K_1 * _K_2 / ( _K_1 + _K_2 );
   }

   /** porosity **///---------------------------------------------------------------------
   static std::vector<double> vec_porosity_b; // read from file
   
   ADscalar Porosity ( const ADscalar& _Po, double _Porosity_b,
		       const ADscalar& _Cr = Cr, const ADscalar& _Pb = Pb )
   {
      return _Porosity_b * ( 1.0 + _Cr * ( _Po - _Pb ) );
   }

   /** read from file **/
   void Read_Layer_Property ( std::string _filename )
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
	 vec_porosity_b.push_back( value );
	 ss >> value;
	 AbsK :: vec_Kx.push_back( value );
	 ss >> value;
	 AbsK :: vec_Ky.push_back( value );
	 ss >> value;
	 AbsK :: vec_Kz.push_back( value );
      }
   }


   /** Transmissibility **///--------------------------------------------------------------
   double Trans_Constant ( double _del_1, double _del_2, double _del_3 )
   {
      return 1.127E-03 * _del_1 * _del_2 / _del_3;
   }

   ADscalar Transmissibility ( double _K, double _trans_constant, const ADscalar& _Kr,
			       const ADscalar& _viscosity, const ADscalar& _volume_factor )
   {
      return _trans_constant * _K * _Kr / _viscosity / _volume_factor;
   }

   /** Specific Weight **///-----------------------------------------------------------------
   ADscalar SpecificWeight ( const ADscalar& _density )
   {
      return _density / 144.0;
   }

   /** Potential **///-----------------------------------------------------------------------
   ADscalar Potential ( const ADscalar& _pressure_1, const ADscalar& _pressure_2,
			const ADscalar& _spec_weight, double _del_depth )
   {
      return _pressure_1 - _pressure_2 - _spec_weight * _del_depth;
   }




   double Find_Po ( double _Po, double _spec_weight, double _del_depth )
   {
      ADscalar tmp_Po_1 =  _Po - _spec_weight * _del_depth;
      // iterate
      ADscalar tmp_Bo = Bo( tmp_Po_1 );
      ADscalar tmp_density = Density_Oil( tmp_Bo );
      ADscalar tmp_spec_weight = SpecificWeight( tmp_density );
      double tmp_Po_2 = tmp_Po_1.value() - (_spec_weight + tmp_spec_weight.value() )/2.0 * _del_depth;
      return tmp_Po_2;
   }

   double Find_Pw ( double _Pw, double _spec_weight, double _del_depth )
   {
      ADscalar tmp_Pw_1 =  _Pw - _spec_weight * _del_depth;

      //std::cout << tmp_Pw_1 << std::endl;
      
      // iterate
      ADscalar tmp_Bw = Bw( tmp_Pw_1 );
      ADscalar tmp_density = Density_Wat( tmp_Bw );
      ADscalar tmp_spec_weight = SpecificWeight( tmp_density );
      double tmp_Pw_2 = tmp_Pw_1.value() - (_spec_weight + tmp_spec_weight.value() )/2.0 * _del_depth;

      //std::cout << tmp_Pw_2 << std::endl;
	    
      return tmp_Pw_2;
   }
}
