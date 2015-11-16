
#include "Solve.hpp"
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

   /* Note: --------------------------------------
    * Give & fix j_gas, solve for Hg & j_liquid
    * --------------------------------------------
    */

   /*
   ADs j_gas( 0.022107 ); // m/s
   ADs j_liquid( -0.05 ); // m/s
   ADs Hg( 0.5 );
   */
   
   ADs j_gas( 0.12 ); // m/s
   ADs j_liquid_CCFL( -0.03 ); // m/s
   j_liquid_CCFL.make_independent( 1 );

   ADs Hg_CCFL( 0.79 );
   Hg_CCFL.make_independent( 0 );

   Solve::Solve_CCFL( density_gas,   density_liquid,
		      viscosity_gas, viscosity_liquid,
		      diameter, g,        surface_tension,
		      j_gas,    j_liquid_CCFL, Hg_CCFL );

   printf( "-------------------------------\n" );
   printf( "-------------------------------\n" );   
   printf( "Hg_CCFL       = %f\n", Hg_CCFL.value() );
   printf( "j_liquid_CCFL = %f\n", j_liquid_CCFL.value() );

   Hg_CCFL.make_constant();
   j_liquid_CCFL.make_constant();
   
   ADs Hg( Hg_CCFL.value()/2.0 );
   ADs j_liquid( j_liquid_CCFL.value()/2.0 );
   j_liquid.make_independent(0);
   
   Solve::Solve_j_liquid (density_gas,   density_liquid,
			  viscosity_gas, viscosity_liquid,
			  diameter, g,        surface_tension,
			  j_gas,    j_liquid_CCFL, Hg, j_liquid );
   
   
   return 0;
}
