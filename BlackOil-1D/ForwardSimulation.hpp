#include "DiscreteProblem.hpp"
#include "Well.hpp"
#include "Jacobian_CSR.hpp"
#include <iostream>
#include <fstream>

#include "Intel_Pardiso.hpp"

void printR( const vector<double>& residual )
{
   for( int i = 0; i < residual.size(); ++i )
   {
      printf( "%10.6e\n", residual[i] );
   }
   n();
}

void printJ( const CSR<>& Jacobian )
{
   for( int row = 0; row < Jacobian.Row().size()-1; ++row )
   {
      int row_beg = Jacobian.Row()[row];
      int row_end = Jacobian.Row()[row+1];

      for( int j = row_beg; j < row_end; ++j )
      {
	 int col = Jacobian.Col()[ j ];
	 double val = Jacobian.NZV()[ j ];

	 printf( "%3d%3d   %10.6e\n", row, col, val );
      }
   } 
}


struct WellDataInfo
{
   int tag; // (0)producer,(1)injector
   int wellindex;
   std::size_t i, j, k;
};
   
struct WellData
{
   int                 mTimeStepIndex;
   double              mTimeStep;
   std::vector<double> mData;
};

class Simulator
{
public:
   Simulator ( );
      
   void Run ();

   const DiscreteProblem         & Get_Problem ()      const  { return problem;         }
   const std::vector<double>     & Get_TimeStep()      const  { return mTimeStep;       }
   const std::vector<double>     & Get_TotalTime()     const  { return mTotalTime;      }
   const std::vector<std::size_t>& Get_RecordTSIndex() const  { return mRecord_TSIndex; }
   
   const std::vector< std::vector<double> >  & Get_Record ()        const { return mRecord;         }
   const std::vector< WellDataInfo >         & Get_WellDataInfo()   const { return mWellDataInfo;   }
   const std::vector< WellData >             & Get_RecordWellData() const { return mRecordWellData; }

   const int Get_DataSetNumber() const { return mRecordWellData.size(); }
   const int GetNg() const { return mWellDataInfo.size(); }
   
private:
   void GenerateTimeStep ( double _init = 0.1, double _grow = 1.05,
			   double _max_step = 30.0, double _limit = 60.0 );

   // void RecordData ();
   // void RecordWellData ( int _time_step_index, double _time_step );

   template< typename VectorType >
   double TwoNorm ( const VectorType& _vector );

private:
   DiscreteProblem problem;
   std::vector< double >      mTimeStep;
   std::vector< double >      mTotalTime;
   std::vector< std::size_t > mRecord_TSIndex;
   
   std::vector< WellDataInfo >            mWellDataInfo;
   std::vector< std::vector<double> >     mRecord;
   std::vector< WellData >                mRecordWellData;
};

Simulator :: Simulator ()
{
   const std::vector< Producer >& vec_producer = problem.GetProducer();
   for( std::size_t i = 0; i < vec_producer.size(); ++i )
   {
      WellDataInfo info;
      info.tag = 0;
      info.wellindex = i;
      info.i = vec_producer[i].I();
      info.j = vec_producer[i].J();
      info.k = vec_producer[i].K();
      mWellDataInfo.push_back( info );
   }
   const std::vector< Injector >& vec_injector = problem.GetInjector();
   for( std::size_t i = 0; i < vec_injector.size(); ++i )
   {
      WellDataInfo info;
      info.tag = 1;
      info.wellindex = i;
      info.i = vec_injector[i].I();
      info.j = vec_injector[i].J();
      info.k = vec_injector[i].K();
      mWellDataInfo.push_back( info );      
   }   
}


void Simulator :: Run ()
{
   //DiscreteProblem problem;

   std::vector< double > negative_newton_update;
   std::vector< double > residual;

   int nRow = problem.GetTotalVarNumber();
   int nCol = problem.GetTotalVarNumber();
   int nNNZ = problem.GetBaseVarNumber() * 8 + problem.GetInjectorNumber();
   CSR<> Jacobian( nRow, nCol, nNNZ );

   // delete this, work on the loop ------------
   // double dT = 1;
   // problem.Update_Accum_Old( dT );
   // problem.Evaluate( dT, Jacobian, residual );
   // ------------------------------------------
/*
   std::ofstream ofs( "Jacobian.txt" );
   ofs << Jacobian.Row().size() << " "
       << Jacobian.Col().size() << " "
       << Jacobian.NZV().size() << std::endl;
   for( int i = 0; i < Jacobian.Row().size(); ++i )
      ofs << Jacobian.Row()[i] << " ";
   ofs << std::endl;
   for( int i = 0; i < Jacobian.Col().size(); ++i )
      ofs << Jacobian.Col()[i] << " ";
   ofs << std::endl;
   for( int i = 0; i < Jacobian.NZV().size(); ++i )
      ofs << Jacobian.NZV()[i] << " ";
   ofs << std::endl;
   
   ofs.close();

   //
   ofs.open( "Residual.txt" );
   for( int i = 0; i < residual.size(); ++i )
   {
      ofs << residual[i] << " ";
   }
   ofs.close();
*/
      
   

   // std::cout << Jacobian << std::endl;
   

   GENSOL :: Intel_Pardiso<> LSolver( Jacobian.nCol(), Jacobian.nNZV() );

   // GenerateTimeStep();
   mTimeStep.push_back( 100 );
   mTotalTime.push_back( 1 );


   int TSIndex_Record = 0;
   // for each time step
   for( std::size_t i = 0; i < mTimeStep.size(); ++i )
   {
      double dT = mTimeStep[i];

      std::cout << " -------------------------------- Time: "
		<< i << " --- " << dT << " --- " << mTotalTime[i] << std::endl;
      problem.Update_Accum_Old( dT );

      const int MaxIteration = 10;
      int count_iter = 0;
      
      // below is newton's iteration within one single time step
      while( count_iter <= MaxIteration )
      {
	 std::cout << " --- --- Iteration: " << count_iter << std::endl;

	 problem.Evaluate( dT, Jacobian, residual );

	 if( count_iter == 2 )
	 {
	    problem.printU();
	    // printR( residual );
	    // printJ( Jacobian );
	    n();
	 }
 
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
	 //
	 double norm_residual = TwoNorm( residual );
	 is_converge_residual = ( norm_residual <= 1.0E-06 );

	 // printf( "norm_update = %.6e\n", norm_update );
	 // printf( "norm_residual = %.6e\n", norm_residual );
	 // problem.printU();
	 if( count_iter == 1 )
	 {
	    for( int i = 0; i < negative_newton_update.size(); ++i )
	    {
	       // printf( "%.6e\n", negative_newton_update[i] );
	    }
	 }
	 n();

	 if( is_converge_residual && is_converge_update )
	 {
	    std::cout << "converged..." << std::endl;
	    problem.printU();
	    
	    break;
	 }

	 problem.Update_Unew( negative_newton_update );
	 problem.Update_Property();
 
	 ++count_iter;
      } // end of newton iteration

      problem.Update_Uold();
   }
}


void Simulator :: GenerateTimeStep ( double _init, double _grow,
				     double _max_step, double _limit )
{
   double step = _init;
   double sum  = _init;
   int count_index = 0;
   
   mTimeStep.push_back( step );
   mTotalTime.push_back( sum );
   count_index++;
      
   double N30  = 30.0;
   while( sum <= _limit )
   {
      step *= _grow;
      if( step > _max_step )  step = _max_step;

      sum += step;
      // if cross boundary of 30
      if( sum > N30 )
      {
	 double remain = sum - N30;
	 sum -=  remain;
	 double tmp_step = step - remain;
	 mTimeStep.push_back( tmp_step );
	 mTotalTime.push_back( sum );
	 mRecord_TSIndex.push_back( count_index );
	 N30 += 30.0;
      }
      else
      {
	 mTimeStep.push_back( step );
	 mTotalTime.push_back( sum );
      }
      count_index++;
   }

   for( std::size_t i = 0; i < mTimeStep.size(); ++i )
   {
   	 std::cout << i << " --- " << mTimeStep[i] << " --- " << mTotalTime[i] << std::endl;
   }

   for( std::size_t i = 0; i < mRecord_TSIndex.size(); ++i )
   {
   	 std::cout << i << " --- " << mRecord_TSIndex[i] << std::endl;
   }
}

/*
void Simulator :: RecordData ()
{
   std::vector< double > data;
   const ADvector& result = problem.GetUnew();
   // 1.
   for( std::size_t i = 0; i < result.size(); ++i )
   {
      // std::cout << i << " --- " << result[i].value() << std::endl;
      data.push_back( result[i].value() );
   }
   mRecord.push_back( data );
}

void Simulator :: RecordWellData( int _time_step_index, double _time_step )
{
   WellData welldata;
   welldata.mTimeStepIndex = _time_step_index;
   welldata.mTimeStep      = _time_step;
   
   // 2.1 producer
   const std::vector< Producer >& vec_producer = problem.GetProducer();
   for( std::size_t i = 0; i < vec_producer.size(); ++i )
   {
      const Producer& well = vec_producer[i];
      double Qo = problem.GetQ( well, 0 ).value();
      welldata.mData.push_back( Qo );
   }
   // 2.2 injector
   const std::vector< Injector >& vec_injector = problem.GetInjector();
   for( std::size_t i = 0; i < vec_injector.size(); ++i )
   {
      std::size_t WellIndex = i + problem.GetBaseVarNumber();
      double pwf = problem.GetUnew()[ WellIndex ].value();
      welldata.mData.push_back( pwf );
   }
   //
   mRecordWellData.push_back( welldata );

   // for( std::size_t i = 0; i < welldata.size(); ++i )
   // {
   //    const WellData& data = welldata[i];
   //    std::cout << "+++ Well#" << i << " +++ " << std::endl
   // 		<< data.tag << "  |  " << data.i << "  |  " << data.j << "  |  " << data.k<< std::endl
   // 		<< data.Qo << "  |  " << data.Qw << "  |  " << data.Ql << "  |  " << std::endl
   // 		<< data.WOR << std::endl << data.Pwf << std::endl
   // 		<< "+++++++++" << std::endl;
   // }
}
*/

template< typename VectorType >
double Simulator :: TwoNorm ( const VectorType& _vector )
{
   double ret = 0.0;
   for( std::size_t i = 0; i < _vector.size(); ++i )
   {
      ret += ( _vector[i] * _vector[i] );
   }
   return sqrt( ret );
}
