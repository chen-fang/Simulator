#include "PropertyCalculator.hpp"
#include "DiscreteProblem.hpp"

int main()
{
   ADs tmp( 0.0 ), tmp2( 0.0 );
   tmp.make_independent( 1 );

   tmp -= tmp2;

   std::cout << tmp << std::endl;

   return 0;
}
