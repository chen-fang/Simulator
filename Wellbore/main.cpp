#include "FluidProperty.hpp"

int main()
{
   ADscalar tmp( 0.0 );
   tmp.make_independent( 1 );

   std::cout << tmp << std::endl;

   return 0;
}
