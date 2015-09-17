// Independent Variables
// Hg, jg, jl, P

template< typename T >
struct EPRI_DF
{
   ADs  L( const ADs& Hg );

   ADs  B1( const ADs& A1 );

   ADs  r(  const ADs& Density_Gas, const ADs& Density_Liquid,
	    const ADs& B1 );
};

template< typename T >
ADs EPRI_DF::L( const ADs& Hg )
{
   ADs tmp( 0.0 );
   tmp = (T)1.15 * pow( Hg, (T)0.45 );
   
   if( tmp.value() > (T)1.0 )
   {
      tmp = (T)1.0;
   }
   return tmp;
}

template< typename T >
ADs EPRI_DF::L( const ADs& A1 )
{
   if( A1.value() > (T)0.8 )
   {
      return (T)0.8;
   }
   return A1;
}

template< typename T >
ADs EPRI_DF::r(  const ADs& Density_Gas, const ADs& Density_Liquid,
		 const ADs& B1 )
{
   ADs tmp( 0.0 );
   return ( (T)1.0 + (T)1.57*Density_Gas/Density_Liquid ) / ( (T)1.0 - B1 );
}
