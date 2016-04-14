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

   ADs j_gas( 0.044 ); // m/s
   ADs j_liquid_CCFL( -0.04 ); // m/s
   j_liquid_CCFL.make_independent( 1 );

   ADs Hg_CCFL( 0.68 );
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

/*   
   ADs Hg( 0.9 );
   double step = j_liquid_CCFL.value() / 100.0;
 
   for( int i = 0; i < 100; ++i )
   {
      ADs j_liquid( i * step );
      //ADs j_liquid( j_liquid_CCFL.value() );
      j_liquid.make_independent(0);
   
      ADs residual( 0.0 );
      residual = Solve::Get_Residual_forTest( density_gas,   density_liquid,
					      viscosity_gas, viscosity_liquid,
					      diameter, g,        surface_tension,
					      j_gas,    j_liquid_CCFL, Hg, j_liquid );

      printf( "Hg = %6.4f (%6.4f)  | j_liquid = %9.6f (%9.6f)  |  residual = %e\n",
	      Hg.value(), Hg_CCFL.value(), j_liquid.value(), j_liquid_CCFL.value(), residual.value() );
   }
*/

   ADs Hg( 0.9 );
   ADs j_liquid( -0.01925 );
   j_liquid.make_independent(0);

   Solve::Solve_j_liquid( density_gas,   density_liquid,
			  viscosity_gas, viscosity_liquid,
			  diameter, g,        surface_tension,
			  j_gas,    j_liquid_CCFL, Hg, j_liquid );
   
   return 0;
}
