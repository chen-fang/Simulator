// Independent Variables
// Hg, jg, jl, P
#include "adetl/scalars/ADscalar.hpp"
#include "adetl/systems/ADvector.hpp"

typedef adetl::ADscalar<> ADs;
typedef adetl::ADvector   ADv;

struct EPRI_DF
{
   ADs u_Liquid( const ADs& j_liquid, const ADs& Hg );
   ADs j_Liquid( const ADs& u_liquid, const ADs& Hg );

   ADs u_Gas( const ADs& j_gas, const ADs& Hg );
   ADs j_Gas( const ADs& u_gas, const ADs& Hg );
   
   // ---------------------------------------------------------
   // ---------------------------------------------------------
   
   ADs L( const ADs& Hg );

   ADs ReyNum( const ADs& density,   const ADs& velocity,
	       const ADs& viscosity, double     diameter );
   
   ADs Re( const ADs& ReyNum_Gas, const ADs& ReyNum_Liquid );

   ADs A1( const ADs& ReyNum );
  
   ADs B1( const ADs& A1 );

   ADs K0( const ADs& Density_Gas, const ADs& Density_Liquid,
	   const ADs& B1 );

   ADs R( const ADs& Density_Gas, const ADs& Density_Liquid,
	  const ADs& B1 );

   ADs C0( const ADs& Density_Gas, const ADs& Density_Liquid,
	   const ADs& ReyNum_Gas,  const ADs& ReyNum_Liquid,
	   const ADs& Hg );

   // ---------------------------------------------------------
   // ---------------------------------------------------------

   ADs C1x( const ADs& B1 );

   ADs C1( const ADs& Hg, const ADs& C1x );

   ADs C5( const ADs& Density_Gas, const ADs& Density_Liquid );

   ADs C2( const ADs& Density_Gas, const ADs& Density_Liquid,
	   const ADs& C5 );

   ADs B2( const ADs& ReyNum_Liquid );

   ADs C3( const ADs& C10, const ADs& B2 );

   double C7( double Dh );

   double C4( double C7 );

   ADs Vgj( const ADs& Density_Gas, const ADs& Density_Liquid,
	    const ADs& ReyNum_Gas,  const ADs& ReyNum_Liquid,
	    const ADs& j_Liquid,    const ADs& j_Liquid_CCFL,
	    const ADs& Hg,          double     surface_tension,
	    double     Dh,          double     g );

   // ---------------------------------------------------------
   // ---------------------------------------------------------

   ADs Z( const ADs& j_Liquid, const ADs& j_Liquid_CCFL );
   
   ADs Jfrx( const ADs& j_Liquid, const ADs& j_Liquid_CCFL );

   ADs Gamma( const ADs& ReyNum_Gas, const ADs& ReyNum_Liquid,
	      const ADs& Jfrx,       double     Dh );

   ADs C10_1( const ADs& ReyNum_Liquid, const ADs& gamma,
	      double     Dh);

   ADs C10_2( const ADs& ReyNum_Liquid, const ADs& j_Liquid,
	      const ADs& j_Liquid_CCFL, double     Dh );

   ADs C10_3( const ADs& ReyNum_Liquid, const ADs& j_Liquid,
	      const ADs& j_Liquid_CCFL, double     Dh );

   ADs C10( const ADs& ReyNum_Gas, const ADs& ReyNum_Liquid,
	    const ADs& j_Liquid,   const ADs& j_Liquid_CCFL,
	    double Dh);
};

ADs EPRI_DF::u_Liquid( const ADs& j_liquid, const ADs& Hg )
{
   ADs ret( 0.0 );
   ret = j_liquid / Hg;
   return ret;
}

ADs EPRI_DF::j_Liquid( const ADs& u_liquid, const ADs& Hg )
{
   ADs ret( 0.0 );
   ret = u_liquid * Hg;
   return ret;
}

ADs EPRI_DF::u_Gas( const ADs& j_gas, const ADs& Hg )
{
   ADs ret( 0.0 );
   ret = j_gas / ( 1.0 - Hg );
   return ret;
}

ADs EPRI_DF::j_Gas( const ADs& u_gas, const ADs& Hg )
{
   ADs ret( 0.0 );
   ret = u_gas * ( 1.0 - Hg );
   return ret;
}

ADs EPRI_DF::L( const ADs& Hg )
{
   ADs tmp( 0.0 );
   tmp = 1.15 * pow( Hg, 0.45 );
   
   if( tmp.value() > 1.0 )
   {
      tmp = 1.0;
   }
   return tmp;
}

ADs EPRI_DF::ReyNum( const ADs& density,   const ADs& velocity,
		     const ADs& viscosity, double     diameter )
{
   ADs ret( 0.0 );
   ret = density*velocity*diameter/viscosity;
   return ret;
}

ADs EPRI_DF::Re( const ADs& ReyNum_Gas, const ADs& ReyNum_Liquid )
{
   return ReyNum_Gas;
}

ADs EPRI_DF::A1( const ADs& ReyNum )
{
   ADs tmp( 0.0 );
   tmp = 1.0 / ( 1.0 + exp(-ReyNum/60000) );
   return tmp;
}

ADs EPRI_DF::B1( const ADs& A1 )
{
   ADs tmp( 0.0 );
   if( A1.value() > 0.8 )
   {
       tmp = 0.8;
       return tmp;
   }
   return A1;
}

ADs EPRI_DF::K0( const ADs& Density_Gas, const ADs& Density_Liquid,
		 const ADs& B1 )
{
   ADs tmp( 0.0 );
   tmp = B1 + ( 1.0 - B1 ) * pow( Density_Gas/Density_Liquid, 0.25 );
   return tmp;
}

ADs EPRI_DF::R(  const ADs& Density_Gas, const ADs& Density_Liquid,
		 const ADs& B1 )
{
   ADs tmp( 0.0 );
   tmp = ( 1.0 + 1.57*Density_Gas/Density_Liquid ) / ( 1.0 - B1 );
   return tmp;
}

ADs EPRI_DF::C0( const ADs& Density_Gas, const ADs& Density_Liquid,
		 const ADs& ReyNum_Gas,  const ADs& ReyNum_Liquid,
		 const ADs& Hg )
{
   ADs l( 0.0 );
   l = L( Hg );
   
   ADs re( 0.0 );
   re = Re( ReyNum_Gas, ReyNum_Liquid );
      
   ADs a1( 0.0 );
   a1 = A1( re );
   
   ADs b1( 0.0 );
   b1 = B1( a1 );
   
   ADs k0( 0.0 );
   k0 = K0( Density_Gas, Density_Liquid, b1 );

   ADs r( 0.0 );
   r = R( Density_Gas, Density_Liquid, b1 );   
   
   ADs ret( 0.0 );
   ret = l / ( k0 + (1.0-k0) * pow( Hg, r ) );
   return ret;
}

// ------------------------------------------------------------
// ------------------------------------------------------------
ADs EPRI_DF::C1x( const ADs& B1 )
{
   return B1;
}

ADs EPRI_DF::C1( const ADs& Hg, const ADs& C1x )
{
   ADs tmp( 0.0 );
   tmp = pow( 1.0 - Hg, C1x );
   return tmp;
}

ADs EPRI_DF::C5( const ADs& Density_Gas, const ADs& Density_Liquid )
{
   ADs tmp( 0.0 );
   tmp = sqrt( 150.0 * Density_Gas/Density_Liquid );
   return tmp;
}

ADs EPRI_DF::C2( const ADs& Density_Gas, const ADs& Density_Liquid, const ADs& C5 )
{
   ADs tmp( 0.0 );
   double ratio = Density_Liquid.value() / Density_Gas.value();
   if( ratio <= 18.0 )
   {
      ADs ratio_ads( 0.0 );
      ratio_ads = Density_Liquid/Density_Gas;
      tmp = 0.4757 * pow( log( ratio_ads), 0.7 );
   }
   else if( ratio > 18.0 && C5.value() >= 1.0 )
   {
      tmp = 1.0;
   }
   else
   {
      tmp = 1.0 - exp( -C5 / (1.0-C5 ) );
   }
   return tmp;
}

ADs EPRI_DF::B2( const ADs& ReyNum_Liquid )
{
   ADs tmp( 0.0 );
   tmp = 1.0 + 0.05 * fabs( ReyNum_Liquid ) / 350000.0;
   ADs ret( 0.0 );
   ret = pow( 1.0/tmp, 0.4 );
   return ret;
}

ADs EPRI_DF::C3( const ADs& C10, const ADs& B2 )
{
   ADs ret( 0.0 );
   ret = 2.0 * pow( C10/2.0, B2 );
   return ret;
}

double EPRI_DF::C7( double Dh )
{
   return std::pow( 0.09144/Dh, 0.6 );
}

double EPRI_DF::C4( double C7 )
{
   double ret;
   if( C7 >= 1.0 )
      ret = 1.0;
   else
      ret = 1.0 - std::exp( -C7/(1.0-C7) );
   return ret;
}

ADs EPRI_DF::Vgj( const ADs& Density_Gas, const ADs& Density_Liquid,
		  const ADs& ReyNum_Gas,  const ADs& ReyNum_Liquid,
		  const ADs& j_Liquid,    const ADs& j_Liquid_CCFL,
		  const ADs& Hg,          double     surface_tension,
		  double     Dh,          double     g )
{
   ADs re( 0.0 );
   re = Re( ReyNum_Gas, ReyNum_Liquid );

   ADs a1( 0.0 );
   a1 = A1( re );

   ADs b1( 0.0 );
   b1 = B1( a1 );

   ADs c1x( 0.0 );
   c1x = C1x( b1 );

   ADs c1( 0.0 );
   c1 = C1( Hg, c1x );

   // std::cout << "--- c1 --- " << c1 << std::endl;

   ADs c5( 0.0 );
   c5 = C5( Density_Gas, Density_Liquid );

   ADs c2( 0.0 );
   c2 = C2( Density_Gas, Density_Liquid, c5 );

   // std::cout << "--- C2 --- " << c2 << std::endl;

   ADs b2( 0.0 );
   b2 = B2( ReyNum_Liquid );

   ADs c10( 0.0 );
   c10 = C10( ReyNum_Gas, ReyNum_Liquid, j_Liquid, j_Liquid_CCFL, Dh );

   // std::cout << "--- C10 --- " << c10 << std::endl;
   
   ADs c3( 0.0 );
   c3 = C3( c10, b2 );

   double c7 = C7( Dh );
   double c4 = C4( c7 );
   
   ADs tmp( 0.0 );
   tmp = pow( (Density_Liquid-Density_Gas)*g*surface_tension/(Density_Liquid*Density_Liquid), 0.25 );
   
   ADs ret( 0.0 );
   ret = 1.41 * tmp * c1 * c2 * c3 * c4;
   
   return ret;
}

ADs EPRI_DF::Z( const ADs& j_Liquid, const ADs& j_Liquid_CCFL )
{
   ADs ret( 0.0 );
   ADs ratio = j_Liquid.value() / j_Liquid_CCFL.value();
   if( ratio.value() < 0.3 )
      ret = 0.8;
   else
      ret = 0.8 - ( ratio - 0.3 );
   return ret;
}

ADs EPRI_DF::Jfrx( const ADs& j_Liquid, const ADs& j_Liquid_CCFL )
{
   ADs z = Z( j_Liquid, j_Liquid_CCFL );
   ADs ratio = j_Liquid / j_Liquid_CCFL;

   // std::cout << "--- Jfrx - ratio --- " << ratio << std::endl;
   
   ADs ret( 0.0 );
   if( ratio.value() != 1.0 )
   {
      ret = pow( 1.0-ratio, z );
   }
   return ret;
}

ADs EPRI_DF::Gamma( const ADs& ReyNum_Gas, const ADs&ReyNum_Liquid,
		    const ADs& Jfrx,     double   Dh )
{
   ADs ret( 0.0 );
   ret = std::pow( 8.0, 0.01905/Dh ) * ReyNum_Gas * Jfrx * exp( -10.0 / fabs( ReyNum_Liquid ) );
   return ret;
}

ADs EPRI_DF::C10_1( const ADs& ReyNum_Liquid, const ADs& gamma,
		    double     Dh )
{
   ADs ret( 0.0 );
   ret = 2.0 * exp( pow( (-ReyNum_Liquid + gamma)/300000, 0.4 ) );
   return ret;
}

ADs EPRI_DF::C10_2( const ADs& ReyNum_Liquid, const ADs& j_Liquid,
		    const ADs& j_Liquid_CCFL, double     Dh )
{
   double c = 0.0381 / Dh;
   c = c*c;

   ADs jfrx( 0.0 );
   jfrx = Jfrx( j_Liquid, j_Liquid_CCFL );
      
   ADs ret( 0.0 );
   ret = -1.72 * pow( fabs(ReyNum_Liquid), 0.035 ) * exp( -fabs(ReyNum_Liquid)/( 35000.0*jfrx+25000.0 ) * c );

   return ret;
}

ADs EPRI_DF::C10_3( const ADs& ReyNum_Liquid, const ADs& j_Liquid,
		    const ADs& j_Liquid_CCFL, double     Dh )
{
   double c = 0.0381 / Dh;
   c = std::pow( c, 0.1 );

   ADs jfrx( 0.0 );
   jfrx = Jfrx( j_Liquid, j_Liquid_CCFL );
   
   ADs ret( 0.0 );
   ret = ( 0.26*jfrx + 0.85*(1.0-jfrx) ) * c * pow( fabs(ReyNum_Liquid), 0.001 );
   return ret;
}

ADs EPRI_DF::C10( const ADs& ReyNum_Gas, const ADs& ReyNum_Liquid,
		  const ADs& j_Liquid,   const ADs& j_Liquid_CCFL,
		  double Dh)
{
   ADs jfrx( 0.0 );
   jfrx = Jfrx( j_Liquid, j_Liquid_CCFL );

   // std::cout << "--- Jfrx --- " << jfrx << std::endl;
   
   ADs gamma( 0.0 );
   gamma = Gamma( ReyNum_Gas, ReyNum_Liquid, jfrx, Dh );
   
   ADs c10_1( 0.0 );
   c10_1 = C10_1( ReyNum_Liquid, gamma, Dh );
   
   ADs c10_2( 0.0 );
   c10_2 = C10_2( ReyNum_Liquid, j_Liquid, j_Liquid_CCFL, Dh );
   
   ADs c10_3( 0,0 );
   c10_3 = C10_3( ReyNum_Liquid, j_Liquid, j_Liquid_CCFL, Dh );

   // std::cout << "--- C10_1 --- " << c10_1 << std::endl;
   
   ADs ret( 0.0 );
   ret = c10_1 + c10_2 + c10_3;
   
   return ret;
}




