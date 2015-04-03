#include "DiscreteProblem.hpp"
#include "Well.hpp"
#include "Jacobian_CSR.hpp"
#include <iostream>

#include "Intel_Pardiso.hpp"

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

template< typename V >
void print_vector ( const V& _vec )
{
   for( std::size_t i =  0; i < _vec.size(); ++i )
   {
      std::cout << i << "\t" << _vec[i] << std::endl;
   }
}

template< typename VectorType >
double TwoNorm ( const VectorType& _vector )
{
   double ret = 0.0;
   for( std::size_t i = 0; i < _vector.size(); ++i )
   {
      ret += ( _vector[i] * _vector[i] );
   }
   return sqrt( ret );
}

struct TestCSR
{
   int nRow;
   int nCol;
   int nNZ;

   std::vector<int> vecRow;
   std::vector<int> vecCol;
   std::vector<double> vecNZ;
};
   

int main ()
{
   TestCSR csr;
   csr.nRow = 2;
   csr.nCol = 2;
   csr.nNZ  = 4;
   csr.vecRow.push_back( 0 );
   csr.vecRow.push_back( 2 );
   csr.vecRow.push_back( 4 );
   
   csr.vecCol.push_back( 0 );
   csr.vecCol.push_back( 1 );
   csr.vecCol.push_back( 0 );
   csr.vecCol.push_back( 1 );
   
   csr.vecNZ. push_back( 1.07492292 );
   csr.vecNZ. push_back( -172909.77247459 );
   csr.vecNZ. push_back( 0.20451442 );
   csr.vecNZ. push_back( 168646.7012178 );
   
   std::vector<double> RHS;
   RHS.push_back( 103745.86348472 );
   RHS.push_back( 67458.68048712 );


   print( "==================================== Summray" );
   print( "nRow", csr.nRow );
   print( "nCol", csr.nCol );
   print( "nNNZ", csr.nNZ );
   print("-----------------row");
   print_vector( csr.vecRow );
   print("-----------------col");
   print_vector( csr.vecCol );
   print("-----------------nzv");
   print_vector( csr.vecNZ );
   print("-----------------RHS");
   print_vector( RHS );

   
   GENSOL::Intel_Pardiso solver( csr.nCol, csr.nNZ );
   std::vector<double> X;
   X.resize( RHS.size() );
   solver.solve( csr.nCol,       &csr.vecNZ[0],
   		 &csr.vecCol[0], &csr.vecRow[0],
   		 &X[0],          &RHS[0]        );


   print("-----------------Result");
   print_vector( X );


   
   DiscreteProblem problem;
   problem.Initialize( 0.2 );

   std::vector< double > residual;
   CSR<> Jacobian( 2, 2, 4 );
   
   problem.Evaluate( 0.2, Jacobian, residual );

   std::cout << std::endl;
   print( "==================================== Jacobian Summray" );
   print( "nRow", Jacobian.nRow() );
   print( "nCol", Jacobian.nCol() );
   print( "nNNZ", Jacobian.nNZV() );
   print("-----------------row");
   print_vector( Jacobian.Row() );
   print("-----------------col");
   print_vector( Jacobian.Col() );
   print("-----------------nzv");
   print_vector( Jacobian.NZV() );
   print("-----------------RHS");
   print_vector( residual );

   
   std::vector< double > update;
   update.resize( residual.size() );

   GENSOL :: Intel_Pardiso LSolver( Jacobian.nCol(), Jacobian.nNZV() );

   // // const int MaxIteration = 5;
   // // int count_iter = 0;
   // // while( count_iter <= MaxIteration )
   // // {
       LSolver.solve( Jacobian.nCol(),    &Jacobian.NZV()[0],
		      &Jacobian.Col()[0], &Jacobian.Row()[0],
		      &update[0],         &residual[0]           );

   //    //print_vector( residual );
       print("-----------------Result");
       print_vector( update );

   //    // double norm_update = TwoNorm( update );
   //    // double norm_residual = TwoNorm( residual );

   // //    std::cout << " ----------------------------------- Iteration: " << count_iter << std::endl;
   // //    std::cout << "norm (update) =\t" << norm_update << std::endl;
   // //    std::cout << "norm (residual) =\t" << norm_residual << std::endl;
      
   // //    if( norm_update <= 1.0E-06 && norm_residual <= 1.0E-06 )
   // //    {
   // // 	 std::cout << "converged..." << std::endl;
   // // 	 break;
   // //    }

   // //    problem.Update_Unknown( update );
   // //    problem.Update_AllProperty();
   // //    problem.Evaluate( 0.2, Jacobian, residual );
   // //    ++count_iter;
   // // }

      
   return -1;
}
