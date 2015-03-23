#include "DiscreteProblem.hpp"
#include <iostream>

template< typename T >
void print ( T a )
{
   std::cout << a << std::endl;
}

template< typename T1, typename T2 >
void print ( T1 a, T2 b )
{
   std::cout << a << "\t" << b << std::endl;
}

int main ()
{
   DiscreteProblem problem;
   problem.Initialize();

   
   return -1;
}
