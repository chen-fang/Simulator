#include "EPRI.hpp"
#include "Jacobian_CSR.h"
#include <iostream>
#include <fstream>

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

   ADs j_liquid( -0.05 ); // m/s
   ADs j_gas_inti( 0.021 );
   // ADs j_gas_1( 0.0 );
   // ADs j_gas_2( 0.0 );
   // j_gas_1 = pow( g*diameter*density_diff, 0.25 );
   // j_gas_2 = pow( density_liquid, 0.25 ) * pow( j_liquid, 0.5 );
   // j_gas = pow( ( C*j_gas_1-j_gas_2)/m, 2.0 ) / pow( density_gas, 0.5 );
   // j_liquid *= -1.0;
   // j_gas.make_independent( 0 );

   ADs Hg( 0.45 );
   Hg.make_independent( 0 );

   double incre_Hg = 1.0E-02;
   double incre_jG = 2.0E-04;
   //j_gas = 0.0221;

   EPRI_DF epri;
   double residual_norm;

   FILE* ofp1;
   FILE* ofp2;
   ofp1 = fopen( "R1.txt", "w" );
   ofp2 = fopen( "R2.txt", "w" );

   for( int count_Hg = 0; count_Hg < 10; ++count_Hg )
   {
      Hg += incre_Hg;

      ADs j_gas( 0.0 );
      j_gas = j_gas_inti;

      for( int count_jg = 0; count_jg < 10; ++count_jg )
      {
	 j_gas += incre_jG;
	 
	 printf( "Hg = %6.3f\njG = %8.5f\n", Hg.value(), j_gas.value() );

	 ADs u_gas( 0.0 );
	 u_gas = epri.u_Gas( j_gas, Hg );
	 // std::cout << "u_Gas = " << u_gas << std::endl;
      
	 ADs u_liquid( 0.0 );
	 u_liquid = epri.u_Liquid( j_liquid, Hg );
	 // std::cout << "u_liquid = " << u_liquid << std::endl;

	 ADs re_gas( 0.0 ), re_liquid( 0.0 );
	 re_gas = epri.ReyNum(    density_gas,    j_gas,    viscosity_gas,    diameter);
	 re_liquid = epri.ReyNum( density_liquid, j_liquid, viscosity_liquid, diameter);

	 // std::cout << "density_gas = " <<density_gas << std::endl;
	 // std::cout << "density_liquid = " <<density_liquid << std::endl;
	 // std::cout << "re_gas = " <<re_gas << std::endl;
	 // std::cout << "re_liquid = " <<re_liquid << std::endl;
	 // std::cout << "Hg = " << Hg << std::endl;
      
	 ADs c0( 0.0 );
	 c0 = epri.C0( density_gas, density_liquid, re_gas, re_liquid, Hg );
	 // std::cout << "C0 =  " << c0 << std::endl;

	 ADs vgj( 0.0 );
	 vgj = epri.Vgj( density_gas, density_liquid,  re_gas,   re_liquid,
			 j_liquid,    j_liquid,
			 Hg,          surface_tension, diameter, g );
	 // std::cout << "Vgj =  " << vgj << std::endl;
	 // double j_liquid_max( 0.0 );
	 // j_liquid_max = -vgj.value() / c0.value();
	 // std::cout << "j_liquid_max =  " << j_liquid_max << std::endl;

	 ADs HgVgj( 0.0 );
	 HgVgj = Hg * vgj;
      
	 ADs HgC0( 0.0 );
	 HgC0 = Hg * c0;

	 ADs residual( 0.0 );
	 residual = c0 * j_liquid + vgj - j_gas * (1.0-Hg*c0) / Hg;
	 // printf( "jL = %7.4f  |  Hg = %3.1f  |  jG = %7.4f  |  residual = %7.4f\n", j_liquid.value(), Hg.value(), j_gas.value(), residual.value() );

      
	 ADs R1( 0.0 ), R2( 0.0 );
	 R1 = Hg * ( c0*j_liquid + vgj ) / ( 1.0 - Hg*c0 ) - j_gas;
	 ADs tmp1( 0.0 ), tmp2( 0.0 );
	 tmp1 = Hg * vgj;
	 tmp2 = Hg * c0;
	 R2 = j_liquid + tmp1.derivative(0)/tmp2.derivative(0) * ( 1.0-Hg*c0 ) + Hg*vgj;

	 printf( "R1 = %14.10f  |  R2 = %14.10f\n\n", R1.value(), R2.value() );

	 fprintf( ofp1, "%14.10f  ", R1.value() );
	 fprintf( ofp2, "%14.10f  ", R2.value() );

	 // residual_norm = std::abs( residual.value() );

	 // double J = residual.derivative(0);
	 // std::cout << "J =  " << J << std::endl;

	 // double update = residual.value() / J;
	 // // std::cout << "update =  " << update << std::endl;

	 // if( update > 0.5 ) update = 0.5;
	 // if( update < -0.5 ) update = -0.5;

	 // std::cout << "jL =  " << j_liquid.value() << std::endl;
	 // std::cout << "update =  " << update << std::endl;

	 // j_gas = j_gas - update;

	 // // if( std::abs(j_liquid.value()) > std::abs(j_liquid_max) )
	 // // {
	 // // 	 j_liquid.value() = j_liquid_max;
	 // // }
	 // if( j_gas.value() <= 0.0 ) j_gas.value() = 0.00001;


	 // if( residual_norm <= 1.0E-08 ) break;
      }
      fprintf( ofp1, "\n" );
      fprintf( ofp2, "\n" );
   }

   fclose( ofp1 );
   fclose( ofp2 );
   
   return 0;
}
