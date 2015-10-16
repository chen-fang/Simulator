#include "EPRI.hpp"
#include "Jacobian_CSR.h"
#include <iostream>

void extract( const ADv& residual,
	      std::vector<double>& R, std::vector< std::vector<double> >& J )
{
   for( int i = 0; i < 2; ++i )
   {
      R[i] = residual[i].value();
      //
      J[i].resize( 2, 0.0 );
      for( int col = 0; col <2; ++col )
      {
	 J[i][col] = residual[i].derivative(col);
      }
   }
}

void Solve22( std::vector<double>&                R,
	      std::vector< std::vector<double> >& J,
	      std::vector<double>&                update )
{
   update.resize( 2, 0.0 );
   if( J[1][0] != 0.0 )
   {
      double tmp = J[0][0] / J[1][0];
      for( int col = 0; col < 2; ++col )
      {
	 J[1][col] *= tmp;
	 J[1][col] -= J[0][col];
      }
      R[1] *= tmp;
      R[1] -= R[0];
   }

   update[1] = R[1] / J[1][1];
   update[0] = ( R[0] - update[1]*J[0][1] ) / J[0][0];
}

int main()
{
   ADs density_liquid( 1000.0 ); // kg/m^3
   ADs density_gas( 1.2 );

   ADs viscosity_liquid( 0.001 ); // Pa.s
   ADs viscosity_gas( 1.81E-05 );

   double diameter = 0.1; // m
   double g = 9.81; // m/s^2
   double surface_tension = 0.072; // N.m

   double m = 1.0;
   double C = 1.0;

   ADs density_diff( 0.0 );
   density_diff = density_liquid - density_gas;

   ADs j_gas( 0.5 );
   ADs j_liquid( 0.0 );
   ADs j_liquid_1( 0.0 );
   ADs j_liquid_2( 0.0 );
   j_liquid_1 = pow( g*diameter*density_diff, 0.25 );
   j_liquid_2 = pow( density_gas, 0.25 ) * pow( j_gas, 0.5 );
   j_liquid = pow( ( C*j_liquid_1-j_liquid_2)/m, 2.0 ) / pow( density_liquid, 0.5 );
   j_liquid *= -1.0;
   j_liquid.make_independent( 1 );
   
   ADs Hg( 0.35 );
   Hg.make_independent( 0 );
   
   EPRI_DF epri;
   double residual_norm;

   int count = 0;
   do
   {
      std::cout << "--------------------------------------------------- Iteration: "<< count++ << std::endl;
      std::cout << "Hg = " << Hg.value() << std::endl;
      std::cout << "jL = " << j_liquid.value() << std::endl;
      
      ADs u_gas( 0.0 );
      u_gas = epri.u_Gas( j_gas, Hg );
      ADs u_liquid( 0.0 );
      u_liquid = epri.u_Liquid( j_liquid, Hg );

      std::cout << "u_liquid" << std::endl << u_liquid << std::endl;

      ADs re_gas( 0.0 );
      re_gas = epri.ReyNum( density_gas, u_gas, viscosity_gas, diameter );
      ADs re_liquid( 0.0 );
      re_liquid = epri.ReyNum( density_liquid, u_liquid, viscosity_liquid, diameter );

      // std::cout << re_gas << std::endl;
      // std::cout << re_liquid << std::endl;

      ADs c0( 0.0 );
      c0 = epri.C0( density_gas, density_liquid, re_gas, re_liquid, Hg );

      ADs vgj( 0.0 );
      vgj = epri.Vgj( density_gas, density_liquid, re_gas, re_liquid,
		      j_liquid, j_liquid, Hg, surface_tension, diameter, g );

      // std::cout << c0 << std::endl;
      // std::cout << vgj << std::endl;

      ADv residual;
      residual.resize( 2, 0.0 );

      residual[0] = Hg * ( c0*j_liquid + vgj )/( 1.0 - Hg*c0 ) - j_gas;

      ADs tmp1( 0.0 );
      tmp1 = Hg * vgj;
      ADs tmp2( 0.0 );
      tmp2 = Hg * c0;

      residual[1] = j_liquid + tmp1.derivative(0)/tmp2.derivative(0) * ( 1.0-tmp2 ) + tmp1;

      ADs c0_dHg = epri.C0_dHg( density_gas, density_liquid, re_gas, re_liquid, Hg );
      
      ADs HgC0_dHg( 0.0 );
      HgC0_dHg = c0 + Hg * c0_dHg;

      // std::cout << "Hg * C0 --- " << std::endl << tmp2 << std::endl;
      // std::cout << "HgC0_dHg ---" << std::endl << HgC0_dHg << std::endl;

      ADs vgj_dHg( 0.0 );
      vgj_dHg = epri.Vgj_dHg( density_gas, density_liquid, re_gas, re_liquid,
			      j_liquid, j_liquid, Hg, surface_tension, diameter, g );

      std::cout << "Vgj --- " << std::endl << vgj << std::endl;
      std::cout << "Vgj_dHg ---" << std::endl << vgj_dHg << std::endl;
      
      ADs HgVgj_dHg( 0.0 );
      HgVgj_dHg = vgj + Hg * vgj_dHg;

      // std::cout << "Hg * Vgj --- " << std::endl << tmp1 << std::endl;
      // std::cout << "HgVgj_dHg ---" << std::endl << HgVgj_dHg << std::endl;

      // std::cout << "R0 | " << residual[0] << std::endl;
      // std::cout << "R1 | " << residual[1] << std::endl;

      std::vector<double> R;
      R.resize( 2 );

      // CSR<> J( 2, 2, 4 );
      // residual.extract_CSR( R, J.Row(), J.Col(), J.NZV() );

      std::vector< std::vector<double> > J;
      J.resize( 2 );
      extract( residual, R, J );

      for( int i = 0; i < 2; ++i )
      {
	 for( int j = 0; j < 2 ; ++j )
	 {
	    // std::cout << J[i][j] << std::endl;
	 }
      }
      

      std::vector<double> update;
      Solve22( R, J, update );

      for( int i = 0; i < 2; ++i )
      {
	 update[i] *= 0.1;
      }

      // std::cout << "Hg_update: " << update[0] << std::endl;
      // std::cout << "jL_update: " << update[1] << std::endl;

      Hg = Hg - update[0];
      j_liquid = j_liquid - update[1];

      if( Hg < 0.0 ) Hg.value() = 0.1;
      if( Hg > 1.0 ) Hg.value() = 0.999;

      residual_norm = std::sqrt( std::pow(residual[0].value(),2.0) + std::pow(residual[1].value(),2.0) );
      // std::cout << "norm: " << residual_norm << std::endl;
      
   }
   // while( residual_norm > 1.0E-03 );
   while( false );
   
   return 0;
}
