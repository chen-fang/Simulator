// Independent Variables
// Hg, jg, jl, P
#include "adetl/scalars/ADscalar.hpp"
#include "adetl/systems/ADvector.hpp"

typedef adetl::ADscalar<> ADs;
typedef adetl::ADvector   ADv;

struct EPRI_DF
{
   ADs L( const ADs& Hg );

   ADs Re( const ADs& Re_Gas, const ADs& Re_Liquid );

   ADs A1( const ADs& ReyNum );
  
   ADs B1( const ADs& A1 );

   ADs K0( const ADs& Density_Gas, const ADs& Density_Liquid,
	   const ADs& B1 );

   ADs r( const ADs& Density_Gas, const ADs& Density_Liquid,
	  const ADs& B1 );

   ADs C0( const ADs& Hg, const ADs& L, const ADs& K0,
	   const ADs& r );

   // ---------------------------------------------------------
   // ---------------------------------------------------------

   ADs C1x( const ADs& B1 );

   ADs C1( const ADs& Hg, const ADs& C1x );

   ADs C5( const ADs& Density_Gas, const ADs& Density_Liquid );

   ADs C2( const ADs& Density_Gas, const ADs& Density_Liquid,
	   const ADs& C5 );
   
};

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

ADs EPRI_DF::Re( const ADs& Re_Gas, const ADs& Re_Liquid )
{
   if( Re_Gas.value() > Re_Liquid.value() )
      return Re_Gas;
   return Re_Liquid;
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

ADs EPRI_DF::r(  const ADs& Density_Gas, const ADs& Density_Liquid,
		 const ADs& B1 )
{
   ADs tmp( 0.0 );
   tmp = ( 1.0 + 1.57*Density_Gas/Density_Liquid ) / ( 1.0 - B1 );
   return tmp;
}

ADs EPRI_DF::C0( const ADs& Hg, const ADs& L, const ADs& K0, const ADs& r )
{
   ADs tmp( 0.0 );
   tmp = L / ( K0 + (1.0-K0) * pow( Hg, r ) );
   return tmp;
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
      tmp = 0.4757 * pow(log( Density_Liquid/Density_Gas, 0.7 ));
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








