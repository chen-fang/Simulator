// Independent Variables
// Hg, jg, jl, P
#include "adetl/scalars/ADscalar.hpp"
#include "adetl/systems/ADvector.hpp"

typedef adetl::ADscalar<> ADs;
typedef adetl::ADvector   ADv;

struct EPRI
{
   // ---------------------------------------------------------
   // --------------------- Function --------------------------
   //
   ADs Min( double a, const ADs& b );
   ADs Max( double a, const ADs& b );

   // ---------------------------------------------------------
   // ------------------------- C0 ----------------------------
   //
   ADs L_CU( const ADs& Hg );
   ADs L_CD( const ADs& Hg );

   ADs ReyNum( const ADs& density,   const ADs& vel,
	       const ADs& viscosity, double     diameter );

   ADs Re_CU( const ADs& ReG, const ADs& ReL );
   ADs Re_CD( const ADs& ReG );

   ADs A1( const ADs& Re );
   ADs B1( const ADs& A1 );

   ADs K0( const ADs& denG, const ADs& denL, const ADs& B1 );
   ADs R(  const ADs& denG, const ADs& denL, const ADs& B1 );

   ADs C0_CU( const ADs& L, const ADs& K0, const ADs& r, const ADs& Hg );

   ADs C0_CD( const ADs& L,   const ADs& K0, const ADs& r,   const ADs& Hg,
	      const ADs& Vgj, const ADs& C1, const ADs& VsG, const ADs& VsL );

   // ---------------------------------------------------------
   // ------------------------ Vgj ----------------------------


   ADs C1x_CU( const ADs& B1 );
   ADs C1x_CD();

   ADs C1( const ADs& Hg, const ADs& C1x );

   ADs C5( const ADs& denG, const ADs& denL );

   ADs C2( const ADs& denG, const ADs& denL, const ADs& C5 );

   ADs B2( const ADs& ReL );

   ADs C3_CU( const ADs& ReL );
   ADs C3_CD( const ADs& C10, const ADs& B2 );

   double C7( double Dh );
   double C4( double C7 );

   ADs Vgj( const ADs& denG, const ADs& denL,
	    const ADs& C1,   const ADs& C2,
	    const ADs& C3,   double     C4,
	    double     Dh,   double     g,  
	    double surface_tension );

   ADs Z( const ADs& VsL, const ADs& VsL_CCFL );

   ADs Jfrx_CD();

   ADs Gamma_CD( const ADs& ReG, const ADs& ReL, const ADs& Jfrx, double Dh );

   ADs C10_1( const ADs& ReL, const ADs& gamma, double Dh );

   ADs C10_2( const ADs& ReL, const ADs& Jfrx, double Dh );

   ADs C10_3( const ADs& ReL, const ADs& Jfrx, double Dh );

   ADs C10( const ADs& C10_1, const ADs& C10_2, const ADs& C10_3 );
};


// Section: Function
ADs EPRI::Min( double a, const ADs& b )
{
   if( b.value() >= a )
   {
      return a;
   }
   return b;
}

// Section: Function
ADs EPRI::Max( double a, const ADs& b )
{
   if( b.value() <= a )
   {
      return a;
   }
   return b;
}


// Section: C0/
ADs EPRI::L_CU( const ADs& Hg )
{
   ADs tmp( 0.0 );
   tmp = 1.15 * pow( Hg, 0.45 );

   return Min( 1.0, tmp );
}

// Section: C0/
ADs EPRI::L_CD( const ADs& Hg )
{
   ADs tmp( 0.0 );
   tmp = 1.05 * pow( Hg, 0.25 );
   
   return Min( 1.0, tmp );
}

// Section: C0/
ADs EPRI::ReyNum( const ADs& density,   const ADs& vel,
		  const ADs& viscosity, double     diameter )
{
   ADs ret( 0.0 );
   ret = density * vel * diameter / viscosity;
   return ret;
}

// Section: C0/
ADs EPRI::Re_CU( const ADs& ReG, const ADs& ReL )
{
   if( ReG.value() > ReL.value() )
   {
      return ReG;
   }
   else
   {
      return ReL;
   }
}

// Section: C0/
ADs EPRI::Re_CD( const ADs& ReG )
{
   return ReG;
}

// Section: C0/
ADs EPRI::A1( const ADs& Re )
{
   ADs ret( 0.0 );
   ret = 1.0 / ( 1.0 + exp(-Re/60000.0) );
   return ret;
}

// Section: C0/
ADs EPRI::B1( const ADs& A1 )
{
   return Min( 0.8, A1 );
}

// Section: C0/
ADs EPRI::K0( const ADs& denG, const ADs& denL, const ADs& B1 )
{
   ADs ret( 0.0 );
   ret = B1 + ( 1.0 - B1 ) * pow( denG/denL, 0.25 );
   return ret;
}

// Section: C0/
ADs EPRI::R( const ADs& denG, const ADs& denL, const ADs& B1 )
{
   ADs ret( 0.0 );
   ret = ( 1.0 + 1.57*denG/denL ) / ( 1.0 - B1 );
   return ret;
}

// Section: C0/
ADs EPRI::C0_CU( const ADs& L, const ADs& K0, const ADs& r, const ADs& Hg )
{
   ADs ret( 0.0 );
   ret = L / ( K0 + (1.0-K0) * pow( Hg, r ) );
   return ret;
}


ADs EPRI::C0_CD( const ADs& L,   const ADs& K0, const ADs& r,   const ADs& Hg,
		 const ADs& Vgj, const ADs& C1, const ADs& VsG, const ADs& VsL )
{
   ADs tmp1( 0.0 );
   tmp1 = C0_CU( L, K0, r, Hg );

   ADs tmp2( 0.0 );
   tmp2 = Vgj/C1 * pow( 1.0-Hg, 0.2 ) / ( fabs(VsG) + fabs(VsL) );

   if( tmp1.value() >= tmp2.value() )
   {
      return tmp1;
   }
   else
   {
      return tmp2;
   }
}

// ------------------------------------------------------------
// ------------------------------------------------------------
//
// Section: Vgj/
ADs EPRI::C1x_CU( const ADs& B1 )
{
   return B1;
}

// Section: Vgj/
ADs EPRI::C1x_CD()
{
   return 0.5;
}

// Section: Vgj/
ADs EPRI::C1( const ADs& Hg, const ADs& C1x )
{
   ADs tmp( 0.0 );
   tmp = pow( 1.0 - Hg, C1x );
   return tmp;
}

// Section: Vgj/
ADs EPRI::C5( const ADs& denG, const ADs& denL )
{
   ADs tmp( 0.0 );
   tmp = sqrt( 150.0 * denG/denL );
   return tmp;
}

// Section: Vgj/
ADs EPRI::C2( const ADs& denG, const ADs& denL, const ADs& C5 )
{
   ADs tmp( 0.0 );
   double ratio = denL.value() / denG.value();
   if( ratio <= 18.0 )
   {
      ADs ratio_ADs( 0.0 );
      ratio_ADs = denL/denG;
      tmp = 0.4757 * pow( log(ratio_ADs), 0.7 );
   }
   else
   {
      if( C5.value() >= 1.0 )
      {
	 tmp = 1.0;
      }
      else
      {
	 tmp = 1.0 - exp( -C5 / (1.0-C5 ) );
      }
   }
   return tmp;
}

// Section: Vgj/
ADs EPRI::B2( const ADs& ReL )
{
   ADs tmp( 0.0 );
   tmp = 1.0 + 0.05 * fabs( ReL ) / 350000.0;
   ADs ret( 0.0 );
   ret = pow( 1.0/tmp, 0.4 );
   return ret;
}

// Section: Vgj/
ADs EPRI::C3_CU( const ADs& ReL )
{
   ADs tmp( 0.0 );
   tmp = 2.0 * exp( -ReL / 300000.0 );

   return Max( 0.5, tmp );
}

// Section: Vgj/
ADs EPRI::C3_CD( const ADs& C10, const ADs& B2 )
{
   ADs tmp( 0.0 );
   tmp = 2.0 * pow( C10/2.0, B2 );

   return Min( 10.0, tmp );
}

// Section: Vgj/
double EPRI::C7( double Dh )
{
   return std::pow( 0.09144/Dh, 0.6 );
}

double EPRI::C4( double C7 )
{
   double ret;
   if( C7 >= 1.0 )
   {
      ret = 1.0;
   }
   else
   {
      ret = 1.0 / ( 1.0 - std::exp(-C7/(1.0-C7)) );
   }
   return ret;
}

ADs EPRI::Vgj( const ADs& denG, const ADs& denL,
	       const ADs& c1,          const ADs& c2,
	       const ADs& c3,          double     c4,
	       double     Dh,          double     g,  double surface_tension )
{
   ADs tmp( 0.0 );
   tmp = pow( (denL-denG) *g *surface_tension /(denL*denL), 0.25 );
   
   ADs ret( 0.0 );
   ret = 1.41 * tmp * c1 * c2 * c3 * c4;
   return ret;
}


ADs EPRI::Z( const ADs& VsL, const ADs& VsL_CCFL )
{
   ADs ret( 0.0 );
   ADs ratio = VsL.value() / VsL_CCFL.value();
   if( ratio.value() < 0.3 )
   {
      ret = 0.8;
   }
   else
      ret = 0.8 - ( VsL/VsL_CCFL - 0.3 );
   return ret;
}

ADs EPRI::Jfrx_CD()
{
   return 1.0;
}

ADs EPRI::Gamma_CD( const ADs& ReG, const ADs& ReL, const ADs& Jfrx, double Dh )
{
   return 0.0;
}

ADs EPRI::C10_1( const ADs& ReL, const ADs& gamma, double Dh )
{
   ADs tmp( 0.0 );
   tmp = pow( (-ReL + gamma)/300000.0, 0.4 );
   
   ADs ret( 0.0 );
   ret = 2.0 * exp( tmp );
   
   return ret;
}

ADs EPRI::C10_2( const ADs& ReL, const ADs& Jfrx, double Dh )
{
   double c = std::pow( 0.0381/Dh, 2.0 );

   ADs ret( 0.0 );
   ret = -1.72 * pow( fabs(ReL), 0.035 ) * exp( -fabs(ReL)/( 35000.0*Jfrx+25000.0 ) * c );
   return ret;
}

ADs EPRI::C10_3( const ADs& ReL, const ADs& Jfrx, double Dh )
{
   double c = std::pow( 0.0381/Dh, 0.1 );

   ADs ret( 0.0 );
   ret = ( 0.26*Jfrx + 0.85*(1.0-Jfrx) ) * c * pow( fabs(ReL), 0.001 );
   return ret;
}

ADs EPRI::C10( const ADs& C10_1, const ADs& C10_2, const ADs& C10_3 )
{
   ADs ret( 0.0 );
   ret = C10_1 + C10_2 + C10_3;
   return ret;
}
