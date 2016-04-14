#include "EPRI.hpp"
#include "Jacobian_CSR.h"

using std::vector;
typedef vector<double> Vector;
struct Solve
{
   static void extract22( const ADv& residual, Vector& R, vector<Vector>& J );
   
   static void Solve22( Vector& R, vector<Vector>& J, Vector& update );

   static void Solve_CCFL( const ADs& dens_gas, const ADs& dens_liquid,
			   const ADs& vis_gas,  const ADs& vis_liquid,
			   double diameter, double g, double surf_tension,
			   const ADs& j_gas,    ADs& j_liquid, ADs& Hg );
   
   static void Solve_j_liquid( const ADs& dens_gas, const ADs& dens_liquid,
			       const ADs& vis_gas,  const ADs& vis_liquid,
			       double diameter,  double g, double surface_tension,
			       const ADs& j_gas, const ADs& j_liquid_CCFL,
			       const ADs& Hg,          ADs& j_liquid       );
};

void Solve::extract22( const ADv& residual, Vector& R, vector<Vector>& J )
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

void Solve::Solve22( Vector& R, vector<Vector>& J, Vector& update )
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

void Solve::Solve_CCFL( const ADs& dens_gas, const ADs& dens_liquid,
			const ADs& vis_gas,  const ADs& vis_liquid,
			double diameter,  double g, double surface_tension,
			const ADs& j_gas, ADs& j_liquid_CCFL, ADs& Hg )
{
   EPRI_DF epri;
   double residual_norm;

   for( int count = 0; count < 20; ++count )
   {
      std::cout << "--------------------------------------------------- Iteration: "<< count << std::endl;

      printf( "jG      = %10.6f\n", j_gas.value() );
      printf( "Hg      = %10.6f\n", Hg.value() );
      printf( "jL_CCFL = %10.6f\n", j_liquid_CCFL.value() );
      
      ADs u_gas( 0.0 ), u_gas_dHg( 0.0 );
      u_gas =     epri.u_Gas(     j_gas, Hg );
      u_gas_dHg = epri.u_Gas_dHg( j_gas, Hg );
      
      ADs u_liquid( 0.0 ), u_liquid_dHg( 0.0 );
      u_liquid =     epri.u_Liquid(     j_liquid_CCFL, Hg );
      u_liquid_dHg = epri.u_Liquid_dHg( j_liquid_CCFL, Hg );

      ADs re_gas( 0.0 ), re_gas_dHg( 0.0 );
      // Reynold's Number is real velocity rather than supificial
      //re_gas =     epri.ReyNum(     dens_gas, j_gas, vis_gas, diameter );
      //re_gas_dHg = epri.ReyNum_dHg( dens_gas, j_gas, vis_gas, diameter );
      
      ADs re_liquid( 0.0 ), re_liquid_dHg( 0.0 );
      re_liquid =     epri.ReyNum(     dens_liquid, j_liquid_CCFL, vis_liquid, diameter );
      re_liquid_dHg = epri.ReyNum_dHg( dens_liquid, j_liquid_CCFL, vis_liquid, diameter );
      // isEqual( re_liquid, re_liquid_dHg );

      ADs c0( 0.0 ), c0_dHg( 0.0 );
      c0 =     epri.C0(     dens_gas, dens_liquid, re_gas, re_liquid, Hg );
      c0_dHg = epri.C0_dHg( dens_gas, dens_liquid, re_gas, re_liquid,
			    re_gas_dHg,  re_liquid_dHg,  Hg );

      ADs vgj( 0.0 ), vgj_dHg( 0.0 );
      vgj =     epri.Vgj(     dens_gas, dens_liquid,  re_gas,   re_liquid,
			      j_liquid_CCFL,    j_liquid_CCFL,
			      Hg,          surface_tension, diameter, g );
      vgj_dHg = epri.Vgj_dHg( dens_gas, dens_liquid,  re_gas,   re_liquid,
			      re_gas_dHg,  re_liquid_dHg,   j_liquid_CCFL, j_liquid_CCFL,
			      Hg,          surface_tension, diameter, g );

      ADs HgVgj( 0.0 ), HgVgj_dHg( 0.0 );
      HgVgj = Hg * vgj;
      HgVgj_dHg = vgj + Hg * vgj_dHg;
      
      ADs HgC0( 0.0 ), HgC0_dHg( 0.0 );
      HgC0 = Hg * c0;
      HgC0_dHg = c0 + Hg * c0_dHg;

      ADv residual;
      residual.resize( 2, 0.0 );

      residual[0] = Hg * ( c0*j_liquid_CCFL + vgj )/( 1.0 - Hg*c0 ) - j_gas;
      residual[1] = j_liquid_CCFL + HgVgj_dHg/HgC0_dHg * ( 1.0-HgC0 ) + HgVgj;

      residual_norm = std::sqrt( std::pow(residual[0].value(),2.0) + std::pow(residual[1].value(),2.0) );

      printf( "residual-1    = %3.6e\n", residual[0].value() );
      printf( "residual-2    = %3.6e\n", residual[1].value() );
      printf( "residual_norm = %3.6e\n", residual_norm );
      
      if( residual_norm <= 1.0E-08 )
      {
	 std::cout << "Converged..." << std::endl;
	 break;
      }

      std::vector<double> R;
      R.resize( 2 );

      // CSR<> J( 2, 2, 4 );
      // residual.extract_CSR( R, J.Row(), J.Col(), J.NZV() );

      std::vector< std::vector<double> > J;
      J.resize( 2 );
      extract22( residual, R, J );

      for( int i = 0; i < 2; ++i )
      {
	 for( int j = 0; j < 2 ; ++j )
	 {
	    // printf( "%10.4f  ", J[i][j] );
	 }
	 // printf( "\n" );
      }
      
      std::vector<double> update;
      Solve22( R, J, update );

      for( int i = 0; i < 2; ++i )
      {
	 // update[i] *= 0.1;
      }
      std::cout << "Hg_update      = " << update[0] << std::endl;
      std::cout << "jL_CCFL_update = " << update[1] << std::endl;

      if( std::abs( update[0] ) > 0.01 )
	 update[0] = update[0] / std::abs( update[0] ) * 0.01;

      if( std::abs( update[1] ) > 0.1 )
	 update[1] = update[1] / std::abs( update[1] ) * 0.1;

      Hg = Hg - update[0];
      j_liquid_CCFL = j_liquid_CCFL - update[1];

   }
}


void Solve::Solve_j_liquid( const ADs& dens_gas, const ADs& dens_liquid,
			    const ADs& vis_gas,  const ADs& vis_liquid,
			    double diameter,  double g, double surface_tension,
			    const ADs& j_gas, const ADs& j_liquid_CCFL,
			    const ADs& Hg,          ADs& j_liquid       )
{
   EPRI_DF epri;
   double residual_norm;

   for( int count = 0; count < 20; ++count )
   {
      std::cout << "--------------------------------------------------- Iteration: "<< count << std::endl;

      printf( "jG      = %10.6f\n", j_gas.value() );
      printf( "Hg      = %10.6f\n", Hg.value() );
      printf( "jL      = %10.6f (%10.6f)\n", j_liquid.value(), j_liquid_CCFL.value() );
      
      ADs u_gas( 0.0 ), u_gas_dHg( 0.0 );
      u_gas =     epri.u_Gas(     j_gas, Hg );
      u_gas_dHg = epri.u_Gas_dHg( j_gas, Hg );
      
      ADs u_liquid( 0.0 ), u_liquid_dHg( 0.0 );
      u_liquid =     epri.u_Liquid(     j_liquid, Hg );
      u_liquid_dHg = epri.u_Liquid_dHg( j_liquid, Hg );

      ADs re_gas( 0.0 ), re_gas_dHg( 0.0 );
      re_gas =     epri.ReyNum(     dens_gas, j_gas, vis_gas, diameter );
      re_gas_dHg = epri.ReyNum_dHg( dens_gas, j_gas, vis_gas, diameter );
      
      ADs re_liquid( 0.0 ), re_liquid_dHg( 0.0 );
      re_liquid =     epri.ReyNum(     dens_liquid, j_liquid, vis_liquid, diameter );
      re_liquid_dHg = epri.ReyNum_dHg( dens_liquid, j_liquid, vis_liquid, diameter );
      // isEqual( re_liquid, re_liquid_dHg );

      ADs c0( 0.0 );
      c0 = epri.C0( dens_gas, dens_liquid, re_gas, re_liquid, Hg );

      ADs vgj( 0.0 );
      vgj = epri.Vgj( dens_gas, dens_liquid,  re_gas,   re_liquid,
		      j_liquid,    j_liquid_CCFL,
		      Hg,          surface_tension, diameter, g );

      ADs HgVgj( 0.0 );
      HgVgj = Hg * vgj;
      
      ADs HgC0( 0.0 );
      HgC0 = Hg * c0;

      ADs residual( 0.0 );
      residual = Hg * ( c0*j_liquid + vgj )/( 1.0 - Hg*c0 ) - j_gas;

      residual_norm = std::abs( residual.value() );

      printf( "residual      = %3.6e\n", residual.value() );
      printf( "residual_norm = %3.6e\n", residual_norm );
      
      if( residual_norm <= 1.0E-08 )
      {
	 std::cout << "Converged..." << std::endl;
	 break;
      }

      double R = residual.value();
      double J = residual.derivative(0);

      double update = R / J;
      std::cout << "jL_update = " << update << std::endl;

      // if( std::abs( update ) > 0.1 )
      //  	 update = update / std::abs( update ) * 0.1;

      j_liquid = j_liquid - update;

   }
}
