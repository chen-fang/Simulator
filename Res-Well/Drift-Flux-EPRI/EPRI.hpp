// Independent Variables
// Hg, jg, jl, P
#include "adetl/scalars/ADscalar.hpp"
#include "adetl/systems/ADvector.hpp"

typedef adetl::ADscalar<> ADs;
typedef adetl::ADvector   ADv;

void isEqual( const ADs& a, const ADs& b )
{
   if( a.derivative(0) == b.value() )
   {
      std::cout << "Equal: " << b.value() << std::endl;
   }
   else
   {
      std::cout << "NOT Equal: " << a.derivative(0) << "  |  " << b.value() << std::endl;
   }
}

struct EPRI_DF
{
   ADs u_Liquid(     const ADs& j_liquid, const ADs& Hg );
   ADs u_Liquid_dHg( const ADs& j_liquid, const ADs& Hg );

   ADs u_Gas(     const ADs& j_gas, const ADs& Hg );
   ADs u_Gas_dHg( const ADs& j_gas, const ADs& Hg );
   
   // ---------------------------------------------------------
   // ---------------------------------------------------------
   //
   ADs L(     const ADs& Hg );
   ADs L_dHg( const ADs& Hg );

   //
   ADs ReyNum(         const ADs& density,   const ADs& velocity,
		       const ADs& viscosity, double     diameter );
   ADs ReyNum_dHg(     const ADs& density,   const ADs& velocity_dHg,
		       const ADs& viscosity, double     diameter );
   ADs ReyNum_Abs_dHg( const ADs& ReyNum,    const ADs& ReyNum_dHg );

   //
   ADs Re(     const ADs& ReyNum_Gas,     const ADs& ReyNum_Liquid );
   ADs Re_dHg( const ADs& ReyNum_Gas_dHg, const ADs& ReyNum_Liquid_dHg );

   ADs A1(     const ADs& ReyNum );
   ADs A1_dHg( const ADs& ReyNum, const ADs& ReyNum_dHg );
  
   ADs B1(     const ADs& A1 );
   ADs B1_dHg( const ADs& A1, const ADs& A1_dHg );

   ADs K0(     const ADs& Density_Gas, const ADs& Density_Liquid,
	       const ADs& B1 );
   ADs K0_dHg( const ADs& Density_Gas, const ADs& Density_Liquid,
	       const ADs& B1_dHg );

   ADs R(     const ADs& Density_Gas, const ADs& Density_Liquid,
	      const ADs& B1 );
   ADs R_dHg( const ADs& Density_Gas, const ADs& Density_Liquid,
	      const ADs& B1,          const ADs& B1_dHg );

   ADs C0( const ADs& Density_Gas, const ADs& Density_Liquid,
	   const ADs& ReyNum_Gas,  const ADs& ReyNum_Liquid,
	   const ADs& Hg );

   ADs C0_dHg( const ADs& Density_Gas,     const ADs& Density_Liquid,
	       const ADs& ReyNum_Gas,      const ADs& ReyNum_Liquid,
	       const ADs& ReyNum_Gas_dHg,  const ADs& ReyNum_Liquid_dHg,    	       
	       const ADs& Hg );
   
   // ---------------------------------------------------------
   // ---------------------------------------------------------

   ADs C1x(     const ADs& B1 );
   ADs C1x_dHg( const ADs& B1_dHg );

   ADs C1(     const ADs& Hg, const ADs& C1x );
   ADs C1_dHg( const ADs& Hg, const ADs& C1x, const ADs& C1, const ADs& C1x_dHg );

   ADs C5( const ADs& Density_Gas, const ADs& Density_Liquid );

   ADs C2( const ADs& Density_Gas, const ADs& Density_Liquid,
	   const ADs& C5 );

   ADs B2(     const ADs& ReyNum_Liquid );
   ADs B2_dHg( const ADs& ReyNum_Liquid, const ADs& ReyNum_Liquid_dHg );

   ADs C3(     const ADs& C10, const ADs& B2 );
   ADs C3_dHg( const ADs& C10, const ADs& B2, const ADs& C10_dHg, const ADs& B2_dHg );

   double C7( double Dh );

   double C4( double C7 );

   ADs Vgj( const ADs& Density_Gas, const ADs& Density_Liquid,
	    const ADs& ReyNum_Gas,  const ADs& ReyNum_Liquid,
	    const ADs& j_Liquid,    const ADs& j_Liquid_CCFL,
	    const ADs& Hg,          double     surface_tension,
	    double     Dh,          double     g );

   ADs Vgj_dHg( const ADs& Density_Gas,     const ADs& Density_Liquid,
		const ADs& ReyNum_Gas,      const ADs& ReyNum_Liquid,
		const ADs& ReyNum_Gas_dHg,  const ADs& ReyNum_Liquid_dHg,
		const ADs& j_Liquid,        const ADs& j_Liquid_CCFL,
		const ADs& Hg,              double     surface_tension,
		double     Dh,              double     g );

   // ---------------------------------------------------------
   // ---------------------------------------------------------

   ADs Z( const ADs& j_Liquid, const ADs& j_Liquid_CCFL );
   
   ADs Jfrx( const ADs& j_Liquid, const ADs& j_Liquid_CCFL );

   ADs Gamma(     const ADs& ReyNum_Gas,     const ADs& ReyNum_Liquid,
		  const ADs& Jfrx,           double     Dh );
   ADs Gamma_dHg( const ADs& ReyNum_Gas,     const ADs& ReyNum_Liquid,
		  const ADs& ReyNum_Gas_dHg, const ADs& ReyNum_Liquid_dHg,
		  const ADs& Jfrx,           double     Dh );

   ADs C10_1(     const ADs& ReyNum_Liquid, const ADs& gamma,
	          double     Dh );
   ADs C10_1_dHg( const ADs& ReyNum_Liquid, const ADs& gamma,
		  double     Dh,            const ADs& ReyNum_Liquid_dHg );

   ADs C10_2(     const ADs& ReyNum_Liquid, const ADs& j_Liquid,
		  const ADs& j_Liquid_CCFL, double     Dh );
   ADs C10_2_dHg( const ADs& ReyNum_Liquid, const ADs& j_Liquid,
		  const ADs& j_Liquid_CCFL, double     Dh,
		  const ADs& ReyNum_Liquid_dHg );

   ADs C10_3(     const ADs& ReyNum_Liquid, const ADs& j_Liquid,
		  const ADs& j_Liquid_CCFL, double     Dh );
   ADs C10_3_dHg( const ADs& ReyNum_Liquid, const ADs& j_Liquid,
		  const ADs& j_Liquid_CCFL, double     Dh,
		  const ADs& ReyNum_Liquid_dHg );

   
   ADs C10(     const ADs& ReyNum_Gas,     const ADs& ReyNum_Liquid,
		const ADs& j_Liquid,       const ADs& j_Liquid_CCFL,
		double     Dh);
   ADs C10_dHg( const ADs& ReyNum_Gas,     const ADs& ReyNum_Liquid,
		const ADs& ReyNum_Gas_dHg, const ADs& ReyNum_Liquid_dHg,
		const ADs& j_Liquid,       const ADs& j_Liquid_CCFL,
		double     Dh);   
};

ADs EPRI_DF::u_Liquid( const ADs& j_liquid, const ADs& Hg )
{
   ADs ret( 0.0 );
   ret = j_liquid / ( 1.0 - Hg );
   return ret;
}

ADs EPRI_DF::u_Liquid_dHg( const ADs& j_liquid, const ADs& Hg )
{
   ADs tmp( 0.0 );
   tmp = 1.0 - Hg;
   ADs ret( 0.0 );
   ret = j_liquid / ( tmp * tmp );
   return ret;
}

ADs EPRI_DF::u_Gas( const ADs& j_gas, const ADs& Hg )
{
   ADs ret( 0.0 );
   ret = j_gas / Hg;
   return ret;
}

ADs EPRI_DF::u_Gas_dHg( const ADs& j_gas, const ADs& Hg )
{
   ADs ret( 0.0 );
   ret = -j_gas / ( Hg * Hg );
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

ADs EPRI_DF::L_dHg( const ADs& Hg )
{
   ADs tmp( 0.0 );
   tmp = 1.15 * pow( Hg, 0.45 );

   ADs ret;
   if( tmp.value() > 1.0 )
   {
      ret = 0.0;
   }
   else
   {
      ret = 1.15*0.45* pow( Hg, 0.45-1.0 );
   }
   return ret;   
}

ADs EPRI_DF::ReyNum( const ADs& density,   const ADs& velocity,
		     const ADs& viscosity, double     diameter )
{
   ADs ret( 0.0 );
   ret = density*velocity*diameter/viscosity;
   return ret;
}
ADs EPRI_DF::ReyNum_dHg( const ADs& density,   const ADs& velocity_dHg,
			 const ADs& viscosity, double     diameter )
{
   ADs ret( 0.0 );
   ret = density*velocity_dHg*diameter/viscosity;
   return ret;
}
ADs EPRI_DF::ReyNum_Abs_dHg( const ADs& ReyNum, const ADs& ReyNum_dHg )
{
   ADs ret( 0.0 );
   if( ReyNum.value() >= 0.0 )
      ret = ReyNum_dHg;
   else
      ret = -ReyNum_dHg;

   return ret;
}


ADs EPRI_DF::Re( const ADs& ReyNum_Gas, const ADs& ReyNum_Liquid )
{
   return ReyNum_Gas;
}
ADs EPRI_DF::Re_dHg( const ADs& ReyNum_Gas_dHg, const ADs& ReyNum_Liquid_dHg )
{
   return ReyNum_Gas_dHg;
}


ADs EPRI_DF::A1( const ADs& ReyNum )
{
   ADs ret( 0.0 );
   ret = 1.0 / ( 1.0 + exp(-ReyNum/60000) );
   return ret;
}
ADs EPRI_DF::A1_dHg( const ADs& ReyNum, const ADs& ReyNum_dHg )
{
   ADs tmp1( 0.0 );
   tmp1 = exp(-ReyNum/60000);

   ADs tmp2( 0.0 );
   tmp2 = 1.0 + tmp1;

   ADs ret( 0.0 );
   ret = tmp1 / 60000.0 * ReyNum_dHg / ( tmp2 * tmp2 );
      
   return ret;
}


ADs EPRI_DF::B1( const ADs& A1 )
{
   ADs ret( 0.0 );
   if( A1.value() > 0.8 )
   {
       ret = 0.8;
       return ret;
   }
   return A1;
}
ADs EPRI_DF::B1_dHg( const ADs& A1, const ADs& A1_dHg )
{
   ADs ret( 0.0 );
   if( A1.value() > 0.8 )
   {
      ret = 0.0;
   }
   else
   {
      ret = A1_dHg;
   }
   return ret;
}


ADs EPRI_DF::K0( const ADs& Density_Gas, const ADs& Density_Liquid,
		 const ADs& B1 )
{
   ADs ret( 0.0 );
   ret = B1 + ( 1.0 - B1 ) * pow( Density_Gas/Density_Liquid, 0.25 );
   return ret;
}
ADs EPRI_DF::K0_dHg( const ADs& Density_Gas, const ADs& Density_Liquid,
		     const ADs& B1_dHg )
{
   ADs ret( 0.0 );
   ret = B1_dHg * ( 1.0 - pow( Density_Gas/Density_Liquid, 0.25 ) );
   return ret;
}

ADs EPRI_DF::R(  const ADs& Density_Gas, const ADs& Density_Liquid,
		 const ADs& B1 )
{
   ADs ret( 0.0 );
   ret = ( 1.0 + 1.57*Density_Gas/Density_Liquid ) / ( 1.0 - B1 );
   return ret;
}
ADs EPRI_DF::R_dHg( const ADs& Density_Gas, const ADs& Density_Liquid,
		    const ADs& B1,          const ADs& B1_dHg )
{
   ADs tmp( 0.0 );
   tmp = 1.0 - B1;

   ADs ret( 0.0 );
   ret = ( 1.0 + 1.57*Density_Gas/Density_Liquid ) * B1_dHg / ( tmp*tmp );
   return ret;
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


ADs EPRI_DF::C0_dHg( const ADs& Density_Gas,     const ADs& Density_Liquid,
		     const ADs& ReyNum_Gas,      const ADs& ReyNum_Liquid,
		     const ADs& ReyNum_Gas_dHg,  const ADs& ReyNum_Liquid_dHg,		     
		     const ADs& Hg )
{
   ADs l( 0.0 ), l_dHg( 0.0 );
   l = L( Hg );
   l_dHg = L_dHg( Hg );
   
   ADs re( 0.0 ), re_dHg( 0.0 );
   re = Re( ReyNum_Gas, ReyNum_Liquid );
   re_dHg = Re_dHg( ReyNum_Gas_dHg, ReyNum_Liquid_dHg );
      
   ADs a1( 0.0 ), a1_dHg( 0.0 );
   a1 = A1( re );
   a1_dHg = A1_dHg( re, re_dHg );
   
   ADs b1( 0.0 ), b1_dHg( 0.0 );
   b1 = B1( a1 );
   b1_dHg = B1_dHg( a1, a1_dHg );
   
   ADs k0( 0.0 ), k0_dHg( 0.0 );
   k0 = K0( Density_Gas, Density_Liquid, b1 );
   k0_dHg = K0_dHg( Density_Gas, Density_Liquid, b1_dHg );

   ADs r( 0.0 ), r_dHg( 0.0 );
   r = R( Density_Gas, Density_Liquid, b1 );
   r_dHg = R_dHg( Density_Gas, Density_Liquid, b1, b1_dHg );
   
   ADs X( 0.0 ), X_dHg( 0.0 );
   X = k0 + (1.0-k0) * pow( Hg, r );
   X_dHg = k0_dHg - k0_dHg * pow( Hg, r ) + (1.0-k0) * pow( Hg, r ) * ( r_dHg*log(Hg) + r/Hg );
   
   ADs ret( 0.0 );
   ret = ( l_dHg * X - X_dHg * l ) / ( X * X );

   return ret;
}

// ------------------------------------------------------------
// ------------------------------------------------------------
ADs EPRI_DF::C1x( const ADs& B1 )
{
   return B1;
}
ADs EPRI_DF::C1x_dHg( const ADs& B1_dHg )
{
   return B1_dHg;
}

ADs EPRI_DF::C1( const ADs& Hg, const ADs& C1x )
{
   ADs tmp( 0.0 );
   tmp = pow( 1.0 - Hg, C1x );
   return tmp;
}

ADs EPRI_DF::C1_dHg( const ADs& Hg, const ADs& C1x, const ADs& C1, const ADs& C1x_dHg )
{
   ADs ret( 0.0 );
   ret = C1 * ( C1x_dHg * log(1.0-Hg) - C1x / (1.0-Hg) );
   return ret;
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
ADs EPRI_DF::B2_dHg( const ADs& ReyNum_Liquid, const ADs& ReyNum_Liquid_dHg )
{
   ADs X( 0.0 );
   X = 1.0 + 0.05 * fabs( ReyNum_Liquid ) / 350000.0;
   ADs re_liquid_abs_dHg( 0.0 );
   re_liquid_abs_dHg = ReyNum_Abs_dHg( ReyNum_Liquid, ReyNum_Liquid_dHg );

   ADs X_dHg( 0.0 );
   X_dHg = 0.05 / 350000.0 * re_liquid_abs_dHg;

   ADs ret( 0.0 );
   ret = 0.4 * pow( 1.0/X, 0.4-1.0 ) * (-X_dHg) / ( X*X );

   return ret;
}


ADs EPRI_DF::C3( const ADs& C10, const ADs& B2 )
{
   ADs ret( 0.0 );
   ret = 2.0 * pow( C10/2.0, B2 );
   return ret;
}
ADs EPRI_DF::C3_dHg( const ADs& C10, const ADs& B2, const ADs& C10_dHg, const ADs& B2_dHg )
{
   ADs c3( 0.0 );
   c3 = C3( C10, B2 );

   ADs ret( 0.0 );
   ret = c3 * ( B2_dHg * log( C10/2.0 ) + B2 / C10 * C10_dHg );

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

   // std::cout << "--- Within Vgj: C10 = " << c10.value() << std::endl;
   
   ADs c3( 0.0 );
   c3 = C3( c10, b2 );

   double c7 = C7( Dh );
   double c4 = C4( c7 );
   
   ADs tmp( 0.0 );
   tmp = pow( (Density_Liquid-Density_Gas)*g*surface_tension/(Density_Liquid*Density_Liquid), 0.25 );
   
   ADs ret( 0.0 );
   ret = 1.41 * tmp * c1 * c2 * c3 * c4;

   // std::cout << "--- Within Vgj: C1 = " << c1.value() << std::endl;
   // std::cout << "--- Within Vgj: C2 = " << c2.value() << std::endl;
   // std::cout << "--- Within Vgj: C3 = " << c3.value() << std::endl;
   // std::cout << "--- Within Vgj: C4 = " << c4 << std::endl;
   
   
   return ret;
}

ADs EPRI_DF::Vgj_dHg( const ADs& Density_Gas,     const ADs& Density_Liquid,
		      const ADs& ReyNum_Gas,      const ADs& ReyNum_Liquid,
		      const ADs& ReyNum_Gas_dHg,  const ADs& ReyNum_Liquid_dHg,
		      const ADs& j_Liquid,        const ADs& j_Liquid_CCFL,
		      const ADs& Hg,              double     surface_tension,
		      double     Dh,              double     g )
{
   ADs re( 0.0 ), re_dHg( 0.0 );
   re = Re( ReyNum_Gas, ReyNum_Liquid );
   re_dHg = Re_dHg( ReyNum_Gas_dHg, ReyNum_Liquid_dHg );

   // isEqual( re, re_dHg );

   ADs a1( 0.0 ), a1_dHg( 0.0 );
   a1 = A1( re );
   a1_dHg = A1_dHg( re, re_dHg );
   // isEqual( a1, a1_dHg );

   ADs b1( 0.0 ), b1_dHg( 0.0 );
   b1 = B1( a1 );
   b1_dHg = B1_dHg( a1, a1_dHg );
   // isEqual( b1, b1_dHg );

   ADs c1x( 0.0 ), c1x_dHg( 0.0 );
   c1x = C1x( b1 );
   c1x_dHg = C1x_dHg( b1_dHg );
   // isEqual( c1x, c1x_dHg );

   ADs c1( 0.0 ), c1_dHg( 0.0 );
   c1 = C1( Hg, c1x );
   c1_dHg = C1_dHg( Hg, c1x, c1, c1x_dHg );
   // isEqual( c1, c1_dHg );

   ADs c5( 0.0 );
   c5 = C5( Density_Gas, Density_Liquid );

   ADs c2( 0.0 );
   c2 = C2( Density_Gas, Density_Liquid, c5 );

   ADs b2( 0.0 ), b2_dHg( 0.0 );
   b2 = B2( ReyNum_Liquid );
   b2_dHg = B2_dHg( ReyNum_Liquid, ReyNum_Liquid_dHg );
   // isEqual( b2, b2_dHg );

   /* test zone begins */
//   /*
   ADs jfrx( 0.0 ), gamma( 0.0 );
   jfrx = Jfrx( j_Liquid, j_Liquid_CCFL );
   gamma = Gamma( ReyNum_Gas, ReyNum_Liquid, jfrx, Dh );
   // std::cout << "--- Within C10_dHg: gamma " << std::endl << gamma << std::endl;   
   
   ADs c10_1( 0.0 ), c10_1_dHg( 0.0 );
   c10_1 = C10_1( ReyNum_Liquid, gamma, Dh );
   c10_1_dHg = C10_1_dHg( ReyNum_Liquid, gamma, Dh, ReyNum_Liquid_dHg );
   // std::cout << "--- Within C10_dHg: ";
   // isEqual( c10_1, c10_1_dHg );

   ADs c10_2( 0.0 ), c10_2_dHg( 0.0 );
   c10_2 = C10_2( ReyNum_Liquid, j_Liquid, j_Liquid_CCFL, Dh );
   c10_2_dHg = C10_2_dHg( ReyNum_Liquid, j_Liquid, j_Liquid_CCFL, Dh, ReyNum_Liquid_dHg );

   // isEqual( c10_2, c10_2_dHg );
   
   ADs c10_3( 0,0 ), c10_3_dHg( 0,0 );
   c10_3 = C10_3( ReyNum_Liquid, j_Liquid, j_Liquid_CCFL, Dh );
   c10_3_dHg = C10_3_dHg( ReyNum_Liquid, j_Liquid, j_Liquid_CCFL, Dh, ReyNum_Liquid_dHg );
   // isEqual( c10_3, c10_3_dHg );

//   */
   /* end of test zone */
   
   ADs c10( 0.0 ), c10_dHg( 0.0 );
   c10 = C10( ReyNum_Gas, ReyNum_Liquid, j_Liquid, j_Liquid_CCFL, Dh );
   c10_dHg = C10_dHg( ReyNum_Gas, ReyNum_Liquid, ReyNum_Gas_dHg, ReyNum_Liquid_dHg,
		      j_Liquid, j_Liquid_CCFL, Dh );

   // isEqual( c10, c10_dHg );
   
   ADs c3( 0.0 ), c3_dHg( 0.0 );
   c3 = C3( c10, b2 );
   c3_dHg = C3_dHg( c10, b2, c10_dHg, b2_dHg );

   double c7 = C7( Dh );
   double c4 = C4( c7 );

   ADs tmp( 0.0 );
   tmp = pow( (Density_Liquid-Density_Gas)*g*surface_tension/(Density_Liquid*Density_Liquid), 0.25 );
   
   ADs ret( 0.0 );
   ret = 1.41 * tmp * c2 * c4 * ( c1_dHg * c3 + c3_dHg * c1 );

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
ADs EPRI_DF::Gamma_dHg( const ADs& ReyNum_Gas,     const ADs& ReyNum_Liquid,
			const ADs& ReyNum_Gas_dHg, const ADs& ReyNum_Liquid_dHg,
			const ADs& Jfrx,           double     Dh )
{
   ADs re_liquid_abs_dHg = ReyNum_Abs_dHg( ReyNum_Liquid, ReyNum_Liquid_dHg );
   ADs ret( 0.0 );
   ret = std::pow( 8.0, 0.01905/Dh ) * exp( -10.0/fabs(ReyNum_Liquid) ) * ( ReyNum_Gas_dHg + 10.0 * ReyNum_Gas / ( ReyNum_Liquid * ReyNum_Liquid ) * re_liquid_abs_dHg );

   return ret;	    
}
      


ADs EPRI_DF::C10_1( const ADs& ReyNum_Liquid, const ADs& gamma,
		    double     Dh )
{
   ADs tmp( 0.0 );
   tmp = pow( (-ReyNum_Liquid + gamma)/350000.0, 0.4 );
   
   ADs ret( 0.0 );
   ret = 2.0 * exp( tmp );
   
   // std::cout << "--- --- --- Within C10_1: ReyNum = " << ReyNum_Liquid.value() << std::endl;
   // std::cout << "--- --- --- Within C10_1: Gamma = " << gamma.value() << std::endl;
   // std::cout << "--- --- --- Within C10_1: tmp0 = " << tmp0.value() << std::endl;
   // std::cout << "--- --- --- Within C10_1: tmp = " << tmp.value() << std::endl;
   
   return ret;
}
ADs EPRI_DF::C10_1_dHg( const ADs& ReyNum_Liquid, const ADs& gamma,
			double     Dh,            const ADs& ReyNum_Liquid_dHg )
{
   ADs c10_1( 0.0 );
   c10_1 = C10_1( ReyNum_Liquid, gamma, Dh );
   ADs tmp( 0.0 );
   tmp = ( -ReyNum_Liquid + gamma ) / 350000.0;

   // test below
   // std::cout << "Within C10_1_dHg: ReyNum_Liquid" << std::endl;
   // isEqual( ReyNum_Liquid, ReyNum_Liquid_dHg );

   
   ADs test1( 0.0 );
   test1 = pow( tmp, 0.4 );

   ADs test2( 0.0 );
   test2 = 0.4 * pow( tmp, 0.4-1.0 ) / (-350000.0) * ReyNum_Liquid_dHg;

   // isEqual( test1, test2 );

   // test above

   ADs ret( 0.0 );
   ret = c10_1 * 0.4 * pow( tmp, 0.4-1.0 ) / (-350000.0) * ReyNum_Liquid_dHg;

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
ADs EPRI_DF::C10_2_dHg( const ADs& ReyNum_Liquid, const ADs& j_Liquid,
			const ADs& j_Liquid_CCFL, double     Dh,
			const ADs& ReyNum_Liquid_dHg )
{
   double c = 0.0381 / Dh;
   c = c*c;

   ADs jfrx( 0.0 );
   jfrx = Jfrx( j_Liquid, j_Liquid_CCFL );

   ADs tmp( 0.0 );
   tmp = -c / ( 35000.0 * jfrx + 25000.0 );

   ADs re_liquid_abs_dHg( 0.0 );
   re_liquid_abs_dHg = ReyNum_Abs_dHg( ReyNum_Liquid, ReyNum_Liquid_dHg );

   ADs ret( 0.0 );
   ret = -1.72 * pow( fabs(ReyNum_Liquid), 0.035 ) * re_liquid_abs_dHg * exp( fabs(ReyNum_Liquid) * tmp ) * ( 0.035 / fabs(ReyNum_Liquid) + tmp );

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

ADs EPRI_DF::C10_3_dHg( const ADs& ReyNum_Liquid, const ADs& j_Liquid,
			const ADs& j_Liquid_CCFL, double     Dh,
			const ADs& ReyNum_Liquid_dHg )
{
   ADs c10_3( 0.0 );
   c10_3 = C10_3( ReyNum_Liquid, j_Liquid, j_Liquid_CCFL, Dh );

   ADs re_liquid_abs_dHg( 0.0 );
   re_liquid_abs_dHg = ReyNum_Abs_dHg( ReyNum_Liquid, ReyNum_Liquid_dHg );

   ADs ret( 0.0 );
   ret = 0.001 * c10_3 / fabs(ReyNum_Liquid) * re_liquid_abs_dHg;

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

   
   ADs ret( 0.0 );
   ret = c10_1 + c10_2 + c10_3;

   // std::cout << "--- --- Within C10: C10_1 = " << c10_1.value() << std::endl;
   // std::cout << "--- --- Within C10: C10_2 = " << c10_2.value() << std::endl;
   // std::cout << "--- --- Within C10: C10_3 = " << c10_3.value() << std::endl;
   
   
   return ret;
}

ADs EPRI_DF::C10_dHg( const ADs& ReyNum_Gas,     const ADs& ReyNum_Liquid,
		      const ADs& ReyNum_Gas_dHg, const ADs& ReyNum_Liquid_dHg,
		      const ADs& j_Liquid,       const ADs& j_Liquid_CCFL,
		      double     Dh)
{
   ADs jfrx( 0.0 );
   jfrx = Jfrx( j_Liquid, j_Liquid_CCFL );
   
   ADs gamma( 0.0 );
   gamma = Gamma( ReyNum_Gas, ReyNum_Liquid, jfrx, Dh );

   ADs c10_1_dHg( 0.0 );
   c10_1_dHg = C10_1_dHg( ReyNum_Liquid, gamma, Dh, ReyNum_Liquid_dHg );
   
   ADs c10_2_dHg( 0.0 );
   c10_2_dHg = C10_2_dHg( ReyNum_Liquid, j_Liquid, j_Liquid_CCFL, Dh, ReyNum_Liquid_dHg );
   
   ADs c10_3_dHg( 0.0 );
   c10_3_dHg = C10_3_dHg( ReyNum_Liquid, j_Liquid, j_Liquid_CCFL, Dh, ReyNum_Liquid_dHg );

   ADs ret( 0.0 );
   ret = c10_1_dHg + c10_2_dHg + c10_3_dHg;

   return ret;
}


