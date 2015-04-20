#pragma once
#include "FluidRockProperty.hpp"
#include "Jacobian_CSR.hpp"
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

std::vector< double > operator+ ( const std::vector<double>& a,
				  const std::vector<double>& b )
{
   std::vector<double> ret;
   std::size_t sz = a.size();
   ret.resize( sz, 0.0 );
   for( std::size_t i = 0; i < sz; ++i )
   {
      ret[i] = a[i] + b[i];
   }
   return ret;
}

std::vector< double > operator* ( double a,
				  const std::vector<double>& b )
{
   std::vector<double> ret;
   std::size_t sz = b.size();
   ret.resize( sz, 0.0 );
   for( std::size_t i = 0; i < sz; ++i )
   {
      ret[i] = a * b[i];
   }
   return ret;
}

template< typename CSR >
std::vector< double > operator* ( const CSR& matrix,
				  const std::vector<double>& vec )
{
   std::vector<double> ret;
   ret.resize( matrix.nRow(), 0.0 );
   
   for( int i = 0; i < matrix.nRow(); ++i )
   {
      for( int j = matrix.Row()[i]; j < matrix.Row()[i+1]; ++j )
      {
	 ret[i] += matrix.NZV()[j] * vec[ matrix.Col()[j] ];
      }
   }
   return ret;
}

template< typename CSR >
std::vector< double > operator* ( const std::vector<double>& vec,
				  const CSR& matrix )
{
   std::vector<double> ret;
   ret.resize( matrix.nCol(), 0.0 );

   for( int i = 0; i < matrix.nRow(); ++i )
   {
      for( int j = matrix.Row()[i]; j < matrix.Row()[i+1]; ++j )
      {
	 int c = matrix.Col()[j];
	 ret[c] += vec[i] * matrix.NZV()[j];
      }
   }
   return ret;
}


class ADJOINT
{
public:
   void Fnext_Dy ( const DiscreteProblem& _problem,
   		   const std::vector<double> & _data, CSR<>& _matrix )
   {
      const ConstProperty& constprop = _problem.GetConstProperty();
      int track_nzv = 0;

      for( int i = 0; i < _problem.GetGridNumber(); ++i )
      {
   	 int PoIndex = _problem.GetPoIndex( i );
   	 int SwIndex = _problem.GetSwIndex( i );
	 
   	 ADscalar Po( 0.0 ), Sw( 0.0 );
   	 Po.value() = _data[ PoIndex ];
   	 Sw.value() = _data[ SwIndex ];
	
   	 Po.make_independent( PoIndex );
   	 Sw.make_independent( SwIndex );

   	 ADscalar So( 0.0 ), pw( 0.0 );
   	 So = 1.0 - Sw;
   	 pw = PROPERTY :: Pw( Po );

   	 ADscalar porosity( 0.0 );
   	 porosity = PROPERTY :: Porosity( Po, constprop.vec_porosity_b[i],
   					  constprop.Cr, constprop.Pb );
   	 ADscalar bo( 0.0 ), bw( 0.0 );
   	 bo = PROPERTY :: Bo( Po, constprop.Bob, constprop.Co, constprop.Pb );
   	 bw = PROPERTY :: Bw( pw, constprop.Bwb, constprop.Cw, constprop.Pb );

   	 ADscalar Fo( 0.0 ), Fw( 0.0 );
   	 Fo = (-1) * porosity * So / bo;
   	 Fw = (-1) * porosity * Sw / bw;

   	 _matrix.Row().push_back( track_nzv );
   	 track_nzv++;
   	 _matrix.NZV().push_back( Fo.derivative( PoIndex ) );
   	 _matrix.Col().push_back( 2*i );
   	 track_nzv++;
   	 _matrix.NZV().push_back( Fw.derivative( PoIndex ) );
   	 _matrix.Col().push_back( 2*i+1 );

   	 _matrix.Row().push_back( track_nzv );
   	 track_nzv++;
   	 _matrix.NZV().push_back( Fo.derivative( SwIndex ) );
   	 _matrix.Col().push_back( 2*i );
   	 track_nzv++;
   	 _matrix.NZV().push_back( Fw.derivative( SwIndex ) );
   	 _matrix.Col().push_back( 2*i+1 );
      }
      
      // Pwf
      for( int i = _problem.GetBaseVarNumber(); i < _problem.GetTotalVarNumber(); ++i )
      {
   	 _matrix.Row().push_back( track_nzv );
   	 track_nzv++;
   	 _matrix.NZV().push_back( 0.0 );
   	 _matrix.Col().push_back( i );
      }
      _matrix.Row().push_back( track_nzv );

      //
      _matrix.nRow() = _problem.GetTotalVarNumber();
      _matrix.nCol() = _problem.GetTotalVarNumber();
      _matrix.nNZV() = track_nzv;
   }

   // Only Qo is recorded at a producer
   void Gi_Dy_Col_Producer ( const DiscreteProblem& _problem,
			     int _producer_index,
			     const std::vector<double>& _recorded_data,
			     std::vector<double>& _gi_dy_col )
   {
      int sz = _problem.GetTotalVarNumber();
      _gi_dy_col.clear();
      _gi_dy_col.resize( sz, 0.0 );
      
      const ConstProperty& constprop = _problem.GetConstProperty();
      const Producer& _producer = _problem.GetProducer()[ _producer_index ];

      std::size_t producer_gd_index = _problem.GetGBIndex( _producer.I(),
							   _producer.J(),
							   _producer.K() );
      int PoIndex = _problem.GetPoIndex( producer_gd_index );
      int SwIndex = _problem.GetSwIndex( producer_gd_index );

      ADscalar Po( 0.0 ), Sw( 0.0 ), pw( 0.0 );
      Po.value() = _recorded_data[ PoIndex ];
      Sw.value() = _recorded_data[ SwIndex ];
	
      Po.make_independent( PoIndex );
      Sw.make_independent( SwIndex );
	 
      pw = PROPERTY :: Pw( Po );

      ADscalar Swd( 0.0 ), Kro( 0.0 );
      Swd = PROPERTY :: Swd( Sw, constprop.Siw, constprop.Sor );
      Kro = PROPERTY :: Kro( Swd );

      ADscalar bo( 0.0 );
      bo = PROPERTY :: Bo( Po, constprop.Bob, constprop.Co, constprop.Pb );

      // Qo
      ADscalar Qo( 0.0 );
      Qo = _producer.WI() * Kro / ( constprop.VisO * bo ) * ( Po - _producer.BHP() );

      _gi_dy_col[ PoIndex ] = Qo.derivative( PoIndex );
      _gi_dy_col[ SwIndex ] = Qo.derivative( SwIndex );
   }      

   // Only Pwf is recorded at an injector
   void Gi_Dy_Col_Injector ( const DiscreteProblem& _problem,
			     int _injector_index,
			     const std::vector<double>& _recorded_data,
			     std::vector<double>& _gi_dy_col )
   {
      int sz = _problem.GetTotalVarNumber();
      
      _gi_dy_col.clear();
      _gi_dy_col.resize( sz, 0.0 );

      const Injector& _injector = _problem.GetInjector()[ _injector_index ];
      std::size_t bhp_index = _problem.GetBHPIndex( _injector.Number() );

      // Pwf
      _gi_dy_col[ bhp_index ] = 1.0;
   }


   // Jacobian
   template< typename CSR >
   void RecoverJacobian ( const std::vector<double>& _data,
			  double _time_step, CSR& _Jacobian )
   {
      DiscreteProblem _problem( _data );

      std::vector<double> tmp_residual;

      _Jacobian.nRow() = _problem.GetTotalVarNumber();
      _Jacobian.nCol() = _problem.GetTotalVarNumber();
      _Jacobian.nNZV() = _problem.GetBaseVarNumber() * 8 + _problem.GetInjectorNumber();

      _problem.Evaluate( _time_step, _Jacobian, tmp_residual );
   }

   template< typename CSR >
   void RecoverJacobian_Dm ( const std::vector<double>& _data, double _time_step, CSR& _jacobian_dm )
   {
      DiscreteProblem _problem( _data );
      _problem.Get_Jacobian_Dm( _time_step, _jacobian_dm );
   }

   template< typename CSR >
   void Lamda_Col ( CSR& _Jacobian, std::vector<double>& _RHS,
		    std::vector<double>& _lamda_col )
   {
      _lamda_col.resize( _RHS.size(), 0.0 );
      
      GENSOL :: Intel_Pardiso<1> LSolver( _Jacobian.nCol(), _Jacobian.nNZV() );

      LSolver.solve( _Jacobian.nCol(),       &_Jacobian.NZV()[0],
		     &_Jacobian.Col()[0],    &_Jacobian.Row()[0],
		     &_lamda_col[0],         &_RHS[0]           );
   }

   void Gi_Dm_i ( const Simulator& _sim,
		  const WellData & _welldata,
		  std::vector< std::vector<double> >& _gi_dm_i )
   {
      int init_timestep_index = _welldata.mTimeStepIndex;

      std::vector< std::vector<double> > lamda_old;

      for( int i = init_timestep_index; i >= 0; --i )
      {
	 print( "======== ========== =========== ========== ======== ==========--- ---", i );
	 
	 double timestep = _sim.Get_TimeStep()[ i ];
	 const std::vector<double>& record_data = _sim.Get_Record()[ i ];

	 print( "Time Step Index: ", i );
	 print( "Time Step:       ", timestep );      

	 CSR<> Jm;
	 RecoverJacobian_Dm( record_data, timestep, Jm );
	 // std::cout << "----------------- Jm" << std::endl;
	 // std::cout << Jm << std::endl;

	 CSR<> Jacobian;
	 RecoverJacobian( record_data, timestep, Jacobian );
	 // std::cout << "----------------- Jacobian" << std::endl;
	 // std::cout << Jacobian << std::endl;


	 if( i == init_timestep_index )
	 {
	    
	    for( int j = 0; j < _sim.GetNg(); ++j )
	    {
	       std::vector< double > gi_dy_col;
	       if( _sim.Get_WellDataInfo()[ j ].tag == 0 )
	       {
		  Gi_Dy_Col_Producer ( _sim.Get_Problem(),
				       _sim.Get_WellDataInfo()[ j ].wellindex,
				       record_data,
				       gi_dy_col );
	       }
	       else
	       {
		  Gi_Dy_Col_Injector ( _sim.Get_Problem(),
				       _sim.Get_WellDataInfo()[j].wellindex,
				       record_data,
				       gi_dy_col );
	       }

	       
	       
	       gi_dy_col = (-1.0) * gi_dy_col;

	       // print( "------------------- gi_dy_col" );
	       // print_vector( gi_dy_col );

	       std::vector<double> lamda_col;
	       Lamda_Col( Jacobian, gi_dy_col, lamda_col );

	       // print( "------------------- lamda_col" );
	       // print_vector( lamda_col );

	       lamda_old.push_back( lamda_col );

	       _gi_dm_i.push_back( lamda_col * Jm );
	       // print( "------------------- gi_dm_col" );
	       // print_vector( _gi_dm_i[j] );
	    }
	 }
	 else
	 {
	    for( int j = 0; j < _sim.GetNg(); ++j )
	    {
	       CSR<> fnext_dy;
	       Fnext_Dy (  _sim.Get_Problem(), record_data, fnext_dy );

	       std::vector<double> RHS = fnext_dy * lamda_old[ j ];
	       RHS = (-1.0) * RHS;
	       
	       std::vector<double> lamda_col;
	       Lamda_Col( Jacobian, RHS, lamda_col );

	       // print( "------------------- lamda_col" );
	       // print_vector( lamda_col );
	       
	       // print( "------------------- tmp" );
	       // print_vector( lamda_col * Jm );
	       
	       _gi_dm_i[ j ] = _gi_dm_i[ j ] + lamda_col * Jm;
	       // print( "------------------- gi_dm_col" );
	       // print_vector( _gi_dm_i[j] );
	       
	       lamda_old[ j ] = std::move( lamda_col );
	    }
	 }
      }
   }


};
