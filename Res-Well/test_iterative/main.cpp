#include "EPRI.hpp"
#include "Jacobian_CSR.h"
#include <iostream>

int main()
{
   ADs density_liquid( 1000.0 ); // kg/m^3
   ADs density_gas( 1.184 );

   ADs viscosity_liquid( 1.0E-03 ); // Pa.s
   ADs viscosity_gas( 1.81E-05 );

   double diameter = 0.15; // m
   double g = 9.81; // m/s^2
   double surface_tension = 0.072; // N.m

   double m = 1.0;
   double C = 1.0;

   ADs density_diff( 0.0 );
   density_diff = density_liquid - density_gas;

   ADs j_gas( 0.0221 ); // m/s
   ADs j_liquid( 0.0 );
   // ADs j_liquid_1( 0.0 );
   // ADs j_liquid_2( 0.0 );
   // j_liquid_1 = pow( g*diameter*density_diff, 0.25 );
   // j_liquid_2 = pow( density_gas, 0.25 ) * pow( j_gas, 0.5 );
   // j_liquid = pow( ( C*j_liquid_1-j_liquid_2)/m, 2.0 ) / pow( density_liquid, 0.5 );
   // j_liquid *= -1.0;
   j_liquid = -0.05;
   j_liquid.make_independent( 1 );
   
   ADs Hg( 0.5 );
   Hg.make_independent( 0 );

   std::cout << "Hg = " << Hg << std::endl;
   std::cout << "jL = " << j_liquid << std::endl;

   
   EPRI_DF epri;
   double residual_norm;

   int count = 0;
   
   do
   {
      std::cout << "--------------------------------------------------- Iteration: "<< count++ << std::endl;
      
      // std::cout << "Hg = " << Hg.value() << std::endl;
      // std::cout << "jL = " << j_liquid.value() << std::endl;
      std::cout << "jG = " << j_gas.value() << std::endl;
      
      ADs u_gas( 0.0 ), u_gas_dHg( 0.0 );
      u_gas =     epri.u_Gas(     j_gas, Hg );
      u_gas_dHg = epri.u_Gas_dHg( j_gas, Hg );

      
      // isEqual( u_gas, u_gas_dHg );
      
      ADs u_liquid( 0.0 ), u_liquid_dHg( 0.0 );
      u_liquid =     epri.u_Liquid(     j_liquid, Hg );
      u_liquid_dHg = epri.u_Liquid_dHg( j_liquid, Hg );

      ADs re_gas( 0.0 ), re_gas_dHg( 0.0 );
      re_gas =     epri.ReyNum(     density_gas, j_gas, viscosity_gas, diameter );
      re_gas_dHg = epri.ReyNum_dHg( density_gas, j_gas, viscosity_gas, diameter );

      // isEqual( re_gas, re_gas_dHg );

      
      ADs re_liquid( 0.0 ), re_liquid_dHg( 0.0 );
      re_liquid =     epri.ReyNum(     density_liquid, j_liquid, viscosity_liquid, diameter );
      re_liquid_dHg = epri.ReyNum_dHg( density_liquid, j_liquid, viscosity_liquid, diameter );

      // isEqual( re_liquid, re_liquid_dHg );

      // std::cout << "density_gas = " <<density_gas << std::endl;
      // std::cout << "density_liquid = " <<density_liquid << std::endl;
      // std::cout << "re_gas = " <<re_gas << std::endl;
      // std::cout << "re_liquid = " <<re_liquid << std::endl;
      // std::cout << "Hg = " << Hg << std::endl;

      ADs c0( 0.0 ), c0_dHg( 0.0 );
      c0 =     epri.C0(     density_gas, density_liquid, re_gas, re_liquid, Hg );
      c0_dHg = epri.C0_dHg( density_gas, density_liquid, re_gas, re_liquid,
			    re_gas_dHg,  re_liquid_dHg,  Hg );

      // std::cout << "[ C0 ] ";
      // isEqual( c0, c0_dHg );
      // std::cout << "C0 = " << c0 << std::endl;

      ADs vgj( 0.0 ), vgj_dHg( 0.0 );
      vgj =     epri.Vgj(     density_gas, density_liquid,  re_gas,   re_liquid,
			      j_liquid,    j_liquid,
			      Hg,          surface_tension, diameter, g );

      // std::cout << "In main: re_liquid" << std::endl;
      // isEqual( re_liquid, re_liquid_dHg );
      
      vgj_dHg = epri.Vgj_dHg( density_gas, density_liquid,  re_gas,   re_liquid,
			      re_gas_dHg,  re_liquid_dHg,   j_liquid, j_liquid,
			      Hg,          surface_tension, diameter, g );

      // std::cout << "[ Vgj ] ";
      // isEqual( vgj, vgj_dHg );
      // std::cout << "Vgj = " << vgj << std::endl;

      ADs j_liquid_max( 0.0 );
      j_liquid_max = -vgj / c0;
      // std::cout << "Max j_liquid = " << j_liquid_max.value() << std::endl;

      ADs HgVgj( 0.0 ), HgVgj_dHg( 0.0 );
      HgVgj = Hg * vgj;
      HgVgj_dHg = vgj + Hg * vgj_dHg;
      // isEqual( HgVgj, HgVgj_dHg );
      
      ADs HgC0( 0.0 ), HgC0_dHg( 0.0 );
      HgC0 = Hg * c0;
      HgC0_dHg = c0 + Hg * c0_dHg;
      // isEqual( HgC0, HgC0_dHg );

      ADs residual;

      residual = j_liquid + HgVgj_dHg/HgC0_dHg * ( 1.0-HgC0 ) + HgVgj;

      residual_norm = std::abs( residual.value() );
      std::cout << "residual_norm =  " << residual_norm << std::endl;

      double J = residual.derivative(0);
      double update = residual.value() / J;

      // std::cout << "Hg_update: " << update[0] << std::endl;
      // std::cout << "jL_update: " << update[1] << std::endl;

      j_gas -= update;

      if( count == 20 ) break;
   }
   while( residual_norm > 1.0E-08 );
   // while( false );
   // while( count < 200 );
   // while( true );

   if( count < 200 )
   {
      // std::cout << "Converged..." << std::endl;
      // std::cout << "Hg = " << Hg.value() << "  |  j_liquid = " << j_liquid.value() << "  |  j_gas = " << j_gas.value() << std::endl;
   }
   else
      // std::cout << "Not converged" << std::endl;

   return 0;
}
