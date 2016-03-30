#include "FluidProperty.hpp"
#include "Discretization.hpp"

int main()
{
   ADs tmp( 0.0 );
   tmp.make_independent( 1 );

   std::cout << tmp << std::endl;

   return 0;
}
