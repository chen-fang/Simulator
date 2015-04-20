#pragma once
#include "Jacobian_CSR.hpp"

struct CSR_Example
{
   CSR_Example ()
   {
      // 1  2  0
      // 0  3  4
      
      matrix.nRow() = 2;
      matrix.nCol() = 3;
      matrix.nNZV() = 4;

      matrix.Row().push_back( 0 );
      matrix.Row().push_back( 2 );
      matrix.Row().push_back( 4 );

      matrix.Col().push_back( 0 );
      matrix.Col().push_back( 1 );
      matrix.Col().push_back( 1 );
      matrix.Col().push_back( 2 );

      matrix.NZV().push_back( 1 );
      matrix.NZV().push_back( 2 );
      matrix.NZV().push_back( 3 );
      matrix.NZV().push_back( 4 );
   }
   
   CSR<> matrix;
};
