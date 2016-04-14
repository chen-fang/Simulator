#pragma once
#include <vector>
#include <iostream>

#include "adetl/scalars/ADscalar.hpp"
#include "adetl/systems/ADvector.hpp"

typedef adetl::ADscalar<> ADs;
typedef adetl::ADvector   ADv;

using std::vector;

// Note:
// Without special notice (i.e. XXX_FD), the standard unit (SI) is applied to variables.

// Density
// SI Unit: [ Kg/m3 ]
// FD Unit: [ Lbm/ft3 ]

// Pressure
// SI Unit: [ Pa ]
// FD Unit: [ psi ]


struct FluidProperty
{
private:
   /** reference pressure [ Pa ] **///--------------------------------------------------
   const double Pb = 1.01325E+05;

   /** compressibility factor [ 1/Pa ] **///--------------------------------------------
   const double Cw = 4.35E-08;

   /** SC density [ Kg/m3 ] **///------------------------------------------------------
   // Standard condition: 15['C] + atmosphere pressure
   const double DenA_SC = 1.225;
   const double DenW_SC = 999.1026;

   /** viscosoty [ Pa.s ] at [ 50C ] **///--------------------------------------------------
   //const ADs VisO = 1.1; // cp
   const ADs VisW = 0.547E-03;
   const ADs VisA = 1.962E-05; 

   /** reference formation volume factor **///--------------------------------------------
   const double Bwb = 1.0;

   /** specific gas constant **///--------------------------------------------------------
   const double Rsp_FD = 53.35;      // [ ft.Lbf/(Lbm.R) ]
   const double Rsp = 287.058; // [ J/(Kg.K) ]

   /** standard gravititional accecleration [ m/s2 ] **///--------------------------------
   const double g = 9.80665;

   double Temperature_C_to_K( double Tc )    { return Tc + 274.15; }
   double Temperature_F_to_R( double Tf )    { return Tf + 459.67; }

   // double Density_SI_to_FD( double Den_SI )  { return Den_SI * 0.062428; }
   // double Density_FD_to_SI( double Den_FD )  { return Den_FD / 0.062428; }

   // double Viscosity_cp_to_SI( double vis )   { return vis * 0.001; }

   // static ADs Bo( const ADs& Po )
   // {
   //    ADs ret( 0.0 );
   //    ret = Bob * ( 1.0 - Co * ( Po - Pb ) );
   //    return ret;
   // }

public:
   ADs Bw( const ADs& Pw )
   {
      ADs ret( 0.0 );
      ret = FluidProperty::Bwb * ( 1.0 - Cw * ( Pw - Pb ) );
      return ret;
   }

   /** viscosity **///----------------------------------------------------------------------
   ADs Viscosity_Wat()   { return VisW; }
   ADs Viscosity_Air()   { return VisA; }

   ADs Viscosity_Mix( const ADs& HL, const ADs& visL, 
		      const ADs& HG, const ADs& visG )
   {
      ADs ret( 0.0 );
      ret = HL * visL + HG * visG;
      return ret;
   }

   /** density **///----------------------------------------------------------------------
   // static ADs Density_Oil( const ADs& Po )
   // {
   //    ADs Bo( 0.0 );
   //    Bo = Bo( Po );
   //    ADs ret( 0.0 );
   //    ret = DensOil_sc / Bo;
   //    return ret;
   // }
   // static ADs Density_Oil( const ADs& Bo )
   // {
   //    ADs ret( 0.0 );
   //    ret = DensOil_sc / Bo;
   //    return ret;
   // }

   ADs Density_Wat( const ADs& Pw )
   {
      ADs ret( 0.0 );
      ret = DenW_SC / Bw(Pw);
      return ret;
   }

   ADs Density_Air( const ADs& P, double Tc )
   {
      // SI Unit
      // pressure:              [ Pa ]
      // specific gas constant: [ J/(Kg.K) ]
      // temperature:           [ C ]
      // density:               [ Kg/m3 ]
      ADs ret( 0.0 );
      ret = P / ( Rsp * Temperature_C_to_K(Tc) );
      return ret;
   }
   ADs Density_Air_FD( const ADs& P, double Tf )
   {
      // Field Unit
      // pressure:              [ psi ]
      // specific gas constant: [ ft.Lbf/(Lbm.R) ]
      // temperature:           [ F ]
      // density:               [ Lbm/ft3 ]
      ADs ret( 0.0 );
      ret = 144.0 * P / ( Rsp_FD * Temperature_F_to_R(Tf) );
      return ret;
   }

   /** Two-Phase Mixture Density **///-------------------------------------------------------
   ADs Density_Mix( const ADs& HL_denL, const ADs& HG_denG )
   {
      ADs ret( 0.0 );
      ret = HL_denL + HG_denG;
      return ret;
   }

   /** Specific Weight **///-----------------------------------------------------------------
   ADs SpecificWeight( const ADs& density )
   {
      // density:  [ Kg/m3 ]
      // SpWeight: [ N/m3 ] or [ Kg/(m2.s2) ]
      ADs ret( 0.0 );
      ret = density * g;
      return ret;
   }
   ADs SpecificWeight_FD( const ADs& density )
   {
      // density:  [ Lbm/ft3 ]
      // SpWeight: [ Lbf/ft3 ]
      ADs ret( 0.0 );
      ret = density;
      return ret;
   }

   /** Reynolds Number **///-----------------------------------------------------------------
   ADs ReynoldsNumber( const ADs& den, const ADs& vel, const ADs& vis, double D )
   {
      // SI Unit
      // density:   [ Kg/m3 ]
      // velocity:  [ m/s ]
      // diameter:  [ m ]
      // viscosity: [ Kg/(m.s) ]
      ADs ret( 0.0 );
      ret = den * vel * D / vis;
      return ret;
   }
   ADs ReynoldsNumber_FD( const ADs& den, const ADs& vel, const ADs vis, double D )
   {
      // Field Unit
      // density:   [ Lbm/ft3 ]
      // velocity:  [ ft/s ]
      // diameter:  [ ft ]
      // viscosity: [ cp ]
      ADs ret( 0.0 );
      ret = 1489.344 * den * vel * D / vis;
      return ret;
   }

   /** Frictional Factor **///--------------------------------------------------------
   ADs Fanning_Friction_Factor( const ADs& Re )
   {
      if( 0.0 == Re.value() )
	return 0.0;

      ADs ret( 0.0 );
      if( Re.value() <= 2100.0 )
      {
	 ret = 16.0 / Re;
      }
      else
      {
	 ret = 0.046 / pow( Re, 0.2 );
      }
      return ret;
   }
   ADs Moody_Friction_Factor( const ADs& Re )
   {
      if( 0.0 == Re.value() )
	return 0.0;

      ADs ret( 0.0 );
      ret = 4.0 * Fanning_Friction_Factor( Re );
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
