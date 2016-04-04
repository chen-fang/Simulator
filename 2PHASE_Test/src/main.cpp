#include "DiscreteProblem.hpp"
#include "Linear/Intel_ILU0.hpp"
#include "Linear/Intel_ILUT.hpp"
#include "Linear/Intel_Prec_GMRES.hpp"
#include "Linear/Intel_Pardiso.hpp"
#include "Linear/DRS_Precond_Solver.hpp"
#include "Linear/R_Precond_Solver.hpp"
#include "Linear/RC_Precond_Solver.hpp"
#include "Linear/TwoStageCombinative.hpp"
#include "Nonlinear/NewtonSolver.hpp"

//typedef GENSOL::Intel_Pardiso                         LINEARSOLVER;
typedef GENSOL::Intel_Prec_GMRES< GENSOL::Intel_ILU0 >  LINEARSOLVER;
typedef GENSOL::NewtonSolver< DiscreteProblem, LINEARSOLVER > STDN;

#include "fastl/containers/pod_vector_unbounded.hpp"
#include <fstream>

void 
dump_solution( const char * filename, const DiscreteProblem::StateVector & v )
{
   std::ofstream strm( filename );
   for ( std::size_t i = 0; i<v.size()-1; ++i )
      strm << v[i].HG.value() << "\t"
	   << v[i].P.value() << "\t"
	   << v[i].VL.value() << "\t"
	   << v[i].VG.value() << std::endl;
   strm << v[v.size()-1].HG.value() << "\t"
	<< v[v.size()-1].P.value() << std::endl;
   strm.close();
}

int main ( )
{
   const std::size_t MAX_NLNITER = 10;
   const double      DT_INIT     = 1.0;
   const double      DT_CUT      = 0.5;
   const double      DT_GROW     = 1.0;
   const double      DT_MAX      = 35.0;
   const double      T_FINAL     = 365.0;
   
   
   const double OPT_NLNITER = ( 3 < MAX_NLNITER ? MAX_NLNITER : 3 );

   const int PHASE_NUM = 2;
   const int TOTAL_CELL_NUM = 2;
   const double dX = 1.0;
   const double D  = 0.5;
   const double theta = 180;

   
   DiscreteProblem model( PHASE_NUM, TOTAL_CELL_NUM, dX, D, theta );

   LINEARSOLVER lnsolver( model.max_num_eqns(), 
			  model.max_num_nnz() );

   std::cout << "MAX = " << model.max_num_eqns() << std::endl;

   STDN newton( model, MAX_NLNITER, 1);

   DiscreteProblem::StateVector uOld, uNew;
   model.initialize_state( uOld );

   //   std::cout << uOld << std::endl;


   double DT   = DT_INIT;
   double time = 0.0;
   //   while ( time < T_FINAL )
   //   {
      uNew = uOld;
      // std::cout << uNew << std::endl;
      
      STDN::report_t stdsmry = newton.solve_timestep( uNew, uOld, DT, model, lnsolver );

      /*
      if ( stdsmry.is_converged )
      {
         uOld =  uNew;
         time += DT;
         if ( stdsmry.niter < OPT_NLNITER ) DT *= DT_GROW;
         DT = std::min( DT, DT_MAX );
         DT = std::min( DT, T_FINAL - time );
         std::cout << "CONVERGED t = " << time << " days" << std::endl;
      }
      else
      {
         DT *= DT_CUT;
         std::cout << "FAILED " << std::endl;
      }
      //   }
   dump_solution( "./output/results.out", uNew );

   */
   return -1;
}
