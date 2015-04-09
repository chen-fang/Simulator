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
template<>
void print_vector ( const ADvector& _vec )
{
   for( std::size_t i =  0; i < _vec.size(); ++i )
   {
      std::cout << i << "\t" << _vec[i].value() << std::endl;
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


int main ()
{
   DiscreteProblem problem;
   print( " --- initial unknown --- " );
   print_vector( problem.GetUnew() );

   std::vector< double > negative_newton_update;
   std::vector< double > residual;

   // int nRow = problem.GetVarNumber();
   // int nCol = problem.GetVarNumber();
   // int nNNZ = problem.GetVarNumber() * 8;
   int nRow = problem.GetVarNumber() + problem.GetInjectorNumber();
   int nCol = problem.GetVarNumber() + problem.GetInjectorNumber();
   int nNNZ = problem.GetVarNumber() * 8 + problem.GetInjectorNumber();
   
   CSR<> Jacobian( nRow, nCol, nNNZ );

   GENSOL :: Intel_Pardiso LSolver( Jacobian.nCol(), Jacobian.nNZV() );

   // for each time step
   double dT = 0.1;
   problem.Update_Accum_Old( dT );

   const int MaxIteration = 10;
   int count_iter = 0;

   // below is newton's iteration within one single time step
   while( count_iter <= MaxIteration )
   {
      std::cout << " ----------------------------------- Iteration: " << count_iter << std::endl;

      problem.Evaluate( dT, Jacobian, residual );
      // print( " *** *** *** *** *** *** *** *** **** " );      
      // print( " *** *** Jacobian Information *** *** " );
      // print( " --- Row --- " );
      // print_vector( Jacobian.Row() );
      // print( " --- Column --- " );
      // print_vector( Jacobian.Col() );
      // print( " --- Jacobian --- " );
      // print_vector( Jacobian.NZV() );
      // print( " *** *** *** *** *** *** *** *** **** " );
      // print( " *** *** *** *** *** *** *** *** **** " );
 
      negative_newton_update.resize( residual.size(), 0.0 );
	 
      LSolver.solve( Jacobian.nCol(),            &Jacobian.NZV()[0],
		     &Jacobian.Col()[0],         &Jacobian.Row()[0],
		     &negative_newton_update[0], &residual[0]           );
      // damping
      for( std::size_t i = 1; i < negative_newton_update.size(); i+=2 )
      {
	 if( negative_newton_update[i] > 0.2 ) negative_newton_update[i] = 0.2;
	 if( negative_newton_update[i] < -0.2) negative_newton_update[i] = -0.2;
      }

      // check convergence
      bool is_converge_update = false;
      bool is_converge_residual = false;

      double norm_update = TwoNorm( negative_newton_update );
      is_converge_update = ( norm_update <= 1.0E-06 );
      std::cout << "norm_update =\t" << norm_update << std::endl;
      //
      double norm_residual = TwoNorm( residual );
      is_converge_residual = ( norm_residual <= 1.0E-06 );
      std::cout << "norm_residual =\t" << norm_residual << std::endl;      

      if( is_converge_residual && is_converge_update )
      {
       	 std::cout << "converged..." << std::endl;
      	 break;
      }

      // update for next iterate
      // print( " --- residual ---" );
      // print_vector( residual );
      print( "--- --- --- negative_newton_update:" );
      print_vector( negative_newton_update ); 

      problem.Update_Unew( negative_newton_update );
      print( "--- --- --- updated unknown " );
      print_vector( problem.GetUnew() );

      problem.Update_Property();
 
      ++count_iter;
   } // end of newton iteration

   problem.Update_Uold();

      
   return -1;
}
