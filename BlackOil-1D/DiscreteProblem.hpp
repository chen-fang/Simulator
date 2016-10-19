#pragma once

#include "FluidRockProperty.hpp"
#include "Gridding.hpp"
#include "Well.hpp"

#include <iostream>

void n()
{
   printf( "\n" );
}

void p( double v )
{
   printf( "%10.6e\n", v );
}

void pp( const ADscalar& v )
{
   printf( "%10.6e\n", v.value() );
}


struct CentralProperty
{
   CentralProperty () : kro(0.0), krw(0.0),
			viscosity_oil(0.0), viscosity_wat(0.0),
			density_oil(0.0), density_wat(0.0),
			porosity(0.0)
   {
   }
   
   void Reset ()
   {
      kro = 0.0;
      krw = 0.0;
      viscosity_oil = 0.0;
      viscosity_wat = 0.0;
      density_oil = 0.0;
      density_wat = 0.0;
      porosity = 0.0;
   }
   
   ADscalar kro, krw;
   ADscalar viscosity_oil, viscosity_wat;
   ADscalar density_oil, density_wat;
   ADscalar porosity;
};

struct InterfcProperty
{			
   double trans_constant;
   double K_interface;
};

class DiscreteProblem
{
public:   
   DiscreteProblem();

private: // for initialization & preparation
   void Add_Wells ();
   void Activate_Uold ();
   void Add_Producer ( double _rate, std::size_t _i, std::size_t j, std::size_t k, double _rw = 0.5 );
   void Add_Injector ( double _bhp,  std::size_t _i, std::size_t j, std::size_t k, double _rw = 0.5 );
   

public:  // for update
   void Activate_U ( const std::vector<double>& _data );
   void Update_Uold ();                                                     // 1
   void Update_Unew ( const std::vector<double>& _negative_newton_update ); // 2
   void Update_Property    ();                                              // 3
   void Update_Accum_Old ( double _time_step );                             // 4

   
private: // for evaluation
   // Accumulation
   double Accum_Constant ( double _time_step, double _del_x, double _del_y, double _del_z ) const;

   ADscalar Accum_Term( double _accum_constant, const ADscalar& _porosity,
			const ADscalar& _saturation, const ADscalar& _density ) const;

   void Evaluate_Accumulation ( double _time_step );

   // Flux
   const CentralProperty& GetUpstream ( std::size_t _gdbk_left,
					const ADscalar& _pressure_left,
					std::size_t _gdbk_right,
					const ADscalar& _pressure_right ) const;

   ADscalar Flux_Term ( const ADscalar& _transmissibility, const ADscalar& _potential) const;
   void     Evaluate_Flux ();  


   // Source & Sink
   ADscalar BHP_Residual ( const Injector& _injector );
   ADscalar Sink_Term ( const Producer& _producer, int _tag ); // for constant BHP wells ONLY
   void     Evaluate_Source_Sink ();

   
public:  // for access & information
   const int GetNx               () const { return grid.Nx;          }
   const int GetNy               () const { return grid.Ny;          }
   const int GetNz               () const { return grid.Nz;          }
   const int GetGridNumber       () const { return GridNumber;       }
   const int GetBaseVarNumber    () const { return VarNumber_Base;   }
   const int GetInjectorNumber   () const { return VarNumber_Added;  }
   const int GetTotalVarNumber   () const { return VarNumber_Total;  }

   const ADvector& GetUold () const  { return mUold; }
   const ADvector& GetUnew () const  { return mUnew; }

   /*
   ADscalar GetQ ( const Producer& _producer, int _tag )
   {
      ADscalar tmp( 0.0 );
      tmp = Sink_Term( _producer, _tag );
      return tmp;
   }
   */
   
   const std::vector< Injector >& GetInjector () const   { return vec_injector; }
   const std::vector< Producer >& GetProducer () const   { return vec_producer; }

   const ConstProperty&   GetConstProperty ()    const   { return constprop;    }
   const InterfcProperty& GetInterfcProperty ()  const   { return interfcprop;  }
   const std::vector< CentralProperty >& GetCentralProperty () const { return vec_centralprop; }

private: // to help functionality
   void ClearCentralProperty ();

public:
   const std::size_t GetPoIndex( std::size_t _GridBlockIndex )                    const;
   const std::size_t GetSwIndex( std::size_t _GridBlockIndex )                    const;
   const std::size_t GetPoIndex( std::size_t _i, std::size_t _j, std::size_t _k ) const;
   const std::size_t GetSwIndex( std::size_t _i, std::size_t _j, std::size_t _k ) const;
   const std::size_t GetGBIndex( std::size_t _i, std::size_t _j, std::size_t _k ) const;

   const std::size_t GetBHPIndex( int _injector_number      )                     const;


public:  // for evaluation  
   template< typename VECTOR, typename CSR >
   void Evaluate ( double _time_step,
		   CSR& _Jacobian, VECTOR& _vec_residual );

private:
   int GetKxColIndex( std::size_t _Kx_GBIndex ) const;
   int GetKyColIndex( std::size_t _Ky_GBIndex ) const;
   int GetKzColIndex( std::size_t _Kz_GBIndex ) const;
   int GetPSColIndex( std::size_t _ps_GBIndex ) const; // porosity

private:
   Cartesian3D                    grid;
   const int                      GridNumber;
   const int                      VarNumber_Base;
   int                            VarNumber_Added;
   int                            VarNumber_Total;
   ConnectionList                 CList;

   ConstProperty                  constprop;
   std::vector< CentralProperty > vec_centralprop;
   InterfcProperty                interfcprop;
   
   std::vector< double >          vec_accum_old;

   std::vector< Injector >  vec_injector; // inject water @ constant rate
   std::vector< Producer >  vec_producer; // produce both @ constant bhp
		       
   ADvector mUnew; // water saturation & oil pressure
   ADvector mUold;
   ADvector residual;

public:
   void printU()
   {
      for( int i = 0; i < mUnew.size(); ++i )
      {
	 pp( mUnew[i] );
      }
      n();
   }

   void hack()
   {
      for( int i = 0; i < GridNumber; ++i )
      {
	 mUnew[ 2*i+0 ] += 1000 * i;
	 mUold[ 2*i+0 ] += 1000 * i;

	 mUnew[ 2*i+1 ] += 0.02 * i;
	 mUold[ 2*i+1 ] += 0.02 * i;
      }
   }

};


// ================================================================================
// ================================================================================
// ================================================================================
// ================================================================================
// ================================================================================
DiscreteProblem :: DiscreteProblem () : grid( 20, 20, 20, 3, 1, 1 ), //20-30-3
					GridNumber      ( grid.Nx * grid.Ny * grid.Nz ),
					VarNumber_Base  ( 2 * GridNumber ),
					VarNumber_Added ( 0 ),
					VarNumber_Total ( VarNumber_Base ),
					CList           ( grid ),
					vec_centralprop ( GridNumber, CentralProperty() ),
					interfcprop     (),
					vec_accum_old   ( VarNumber_Base, 0.0 ),
					mUnew           ( VarNumber_Base, 0.0 ),
					mUold           ( VarNumber_Base, 0.0 ),
					residual        ( VarNumber_Base, 0.0 )
{
   interfcprop.trans_constant = PROPERTY::Trans_Constant( grid.Dy, grid.Dz, grid.Dx );
   interfcprop.K_interface    = PROPERTY::K_Interface( constprop.Kx );

   Add_Wells();

   Activate_Uold();
   mUnew = mUold;

   hack();
   Update_Property();
}

// --- ---
void DiscreteProblem :: Add_Wells()
{
   // inject water with constant rate
   Add_Injector( -80, 0, 0, 0 );

   // producer with constant bottom pressure
   // Add_Producer( 2.0E+04, 2, 0, 0 );
}
// --- ---
void DiscreteProblem :: Activate_Uold()
{
   // 1. Initialize Po and Sw
   for( int i = 0; i < GridNumber; ++i )
   {
      const std::size_t PoIndex = GetPoIndex( i );
      const std::size_t SwIndex = GetSwIndex( i );
      mUold[ PoIndex ].value() = constprop.Poi.value();
      mUold[ PoIndex ].make_independent( PoIndex );
      
      mUold[ SwIndex ].value() = constprop.Siw.value();
      mUold[ SwIndex ].make_independent( SwIndex );
   }
   // 2. Initialize BHP as unknowns
   for( int i = VarNumber_Base; i < VarNumber_Total; ++i )
   {
      mUold[ i ].value() = constprop.Poi.value();
      mUold[ i ].make_independent( i );

      // std::cout << "Initializing BHP at " << i << std::endl;
      // std::cout << "Initializing BHP's value: " << mUold[ i ].value() << std::endl;
   }
}
// --- ---
void DiscreteProblem :: Activate_U ( const std::vector<double>& _data )
{
   mUold.clear();
   mUold.resize( _data.size(), 0.0 );
   
   // 1. Initialize Po and Sw
   for( int i = 0; i < VarNumber_Base; )
   {
      mUold[ i ].value() = _data[ i ];
      mUold[ i ].make_independent( i );
      i++;
      
      mUold[ i ].value() = _data[ i ];
      mUold[ i ].make_independent( i );
      i++;
   }
   // 2. Initialize BHP as unknowns
   for( int i = VarNumber_Base; i < VarNumber_Total; ++i )
   {
      mUold[ i ].value() = _data[ i ];
      mUold[ i ].make_independent( i );

      // std::cout << "Initial BHP at " << i << std::endl;
      // std::cout << "Initial BHP's value changes to: " << mUold[ i ].value() << std::endl;
   }
}


void DiscreteProblem :: Add_Injector ( double _rate,
				       std::size_t _i, std::size_t _j, std::size_t _k, double _rw )
{
   const std::size_t BlockIndex = grid.GetGridBlockIndex( _i, _j, _k );
   Injector well( _rate, _i, _j, _k, _rw,
		  grid.Dx, grid.Dy,
		  constprop.Kx, constprop.Kx );

   well.Number() = vec_injector.size();
   vec_injector.push_back( well );
   // expand
   VarNumber_Added++;
   VarNumber_Total++;
   
   ADscalar tmp( 0.0 );
   mUold.push_back( tmp );
   residual.push_back( tmp );

   // std::cout << "mUold expanded to size: " << mUold.size() << std::endl;
   // std::cout << "residual expanded to size: " << residual.size() << std::endl; 
   
}
//
void DiscreteProblem :: Add_Producer ( double _bhp,
				       std::size_t _i,std::size_t _j,std::size_t _k,
				       double _rw )
{
   std::size_t BlockIndex = grid.GetGridBlockIndex( _i, _j, _k );
   Producer well( _bhp, _i, _j, _k, _rw,
		  grid.Dx, grid.Dy,
		  constprop.Kx, constprop.Kx );

   well.Number() = vec_producer.size();
   vec_producer.push_back( well );
}


//1
void DiscreteProblem :: Update_Uold ()
{
   mUold = mUnew;
}
//2
void DiscreteProblem :: Update_Unew ( const std::vector<double>& _negative_newton_update )
{
  for( int i = 0; i < VarNumber_Total; ++i )
   {
      mUnew[i].value() -= _negative_newton_update[i];
   }
}
//3
void DiscreteProblem :: Update_Property ()
{
   for( int i = 0; i < GridNumber; ++i )
   {
      CentralProperty&  prop = vec_centralprop[ i ];
   
      const std::size_t PoIndex = GetPoIndex( i );
      const std::size_t SwIndex = GetSwIndex( i );
      const ADscalar&  _Po = mUnew[ PoIndex ];
      const ADscalar&  _Sw = mUnew[ SwIndex ];

      prop.kro = PROPERTY::Kro( _Sw );
      prop.krw = PROPERTY::Krw( _Sw );

      prop.density_oil = PROPERTY::Density_Oil( constprop.Co, _Po );
      prop.density_wat = PROPERTY::Density_Wat();
   
      prop.viscosity_oil = constprop.VisO;
      prop.viscosity_wat = constprop.VisW;
   
      prop.porosity = constprop.porosity;
   }
}
// 4
void DiscreteProblem :: Update_Accum_Old ( double _time_step )
{
   //std::cout << "Accum_Old Term Updated..." << std::endl;
   const double accum_const = Accum_Constant( _time_step, grid.Dx, grid.Dy, grid.Dz );
   for( int i = 0; i < GridNumber; ++i )
   {      
      const CentralProperty& prop = vec_centralprop[ i ];
      const std::size_t PoIndex = GetPoIndex( i );
      const std::size_t SwIndex = GetSwIndex( i );

      // oil phase
      ADscalar So( 0.0 );
      So = 1.0 - mUnew[ SwIndex ];
      vec_accum_old[ PoIndex ] = accum_const * constprop.porosity * So.value() * prop.density_oil.value();

      // water phase
      const ADscalar& Sw = mUnew[ SwIndex ];
      vec_accum_old[ SwIndex ] = accum_const * constprop.porosity * Sw.value() * prop.density_wat.value();
   }
}





template< typename VECTOR, typename CSR >
void DiscreteProblem :: Evaluate ( double _time_step,
				   CSR& _Jacobian, VECTOR& _vec_residual )
{
   residual.clear();
   residual.resize( VarNumber_Total, 0.0 );
   
   Evaluate_Accumulation( _time_step );
   Evaluate_Flux();
   Evaluate_Source_Sink();
   
   residual.extract_CSR( _vec_residual,
    			 _Jacobian.Row(), _Jacobian.Col(), _Jacobian.NZV() );
}




// ================================================================== Accumulation
// ================================================================== Accumulation
// ================================================================== Accumulation
// --- ---
double DiscreteProblem :: Accum_Constant ( double _time_step,
					   double _del_x, double _del_y, double _del_z ) const
{
   return _del_x * _del_y * _del_z / _time_step;
}
// --- ---
ADscalar DiscreteProblem :: Accum_Term ( double _accum_constant,
					 const ADscalar& _porosity,
					 const ADscalar& _saturation,
					 const ADscalar& _density ) const
{
   ADscalar ret( 0.0 );
   ret = _accum_constant * ( _porosity * _saturation * _density );
   return ret;
}

void DiscreteProblem :: Evaluate_Accumulation ( double _time_step )
{
   const double accum_constant = Accum_Constant( _time_step, grid.Dx, grid.Dy, grid.Dz );
   const double DV = grid.Dx * grid.Dy * grid.Dz;

   for( int i = 0; i < GridNumber; ++i )
   {
      // std::cout << " ------------------------------ " << i << std::endl;
      // accumulate
      const CentralProperty& prop = vec_centralprop[ i ];
      const std::size_t PoIndex = GetPoIndex( i );
      const std::size_t SwIndex = GetSwIndex( i );

      ADscalar accum_oil( 0.0 );
      ADscalar accum_wat( 0.0 );

      // oil phase
      // 1. use old accumulation value;
      ADscalar accum_old_oil( 0.0 );
      accum_old_oil = vec_accum_old[ PoIndex ] / DV;

      residual[ PoIndex ] -= accum_old_oil;


      // 2. update oil accumulation term, which will become the "old" one in the next iterate
      accum_oil = Accum_Term( accum_constant, prop.porosity, 1.0-mUnew[SwIndex], prop.density_oil ) / DV;
      residual[ PoIndex ] += accum_oil;

      
      printf( "%15s = %.6e\n", "accum_old_oil", accum_old_oil.value() );
      printf( "%15s = %.6e\n", "accum_oil",     accum_oil.value() );
      // printf( "%15s%d = %.6e\n", "residual_accum_oil at cell# ",
      //  	      i, residual[ PoIndex ].value() );
      n();
      

      // water phase
      // 1.
      residual[ SwIndex ] -= vec_accum_old[ SwIndex ];
      residual[ SwIndex ] /= DV;
      // 2.
      accum_wat = Accum_Term( accum_constant, prop.porosity, mUnew[SwIndex], prop.density_wat ) / DV;
      residual[ SwIndex ] += accum_wat;
   }
}



// ================================================================== Flux
// ================================================================== Flux
// ================================================================== Flux
// --- ---
const CentralProperty&
DiscreteProblem ::GetUpstream ( std::size_t _gdbk_left,
				const ADscalar& _pressure_left,
				std::size_t _gdbk_right,
				const ADscalar& _pressure_right ) const
{
   // take left
   /* Note
    * I'm lazy. I'll take right in case of equality.
    */
   if( _pressure_left > _pressure_right )
   {
      return vec_centralprop[ _gdbk_left ];
   }
   // take right
   return vec_centralprop[ _gdbk_right ];
}
// --- ---
ADscalar DiscreteProblem :: Flux_Term ( const ADscalar& _transmissibility,
					const ADscalar& _potential) const
{
   ADscalar ret( 0.0 );
   ret = _transmissibility * _potential;
   return ret;
}
// --- ---
void DiscreteProblem :: Evaluate_Flux ()
{
   const double DV = grid.Dx * grid.Dy * grid.Dz;
   
   std::size_t left_gdbk, right_gdbk;
   double trans_constant = interfcprop.trans_constant;
   double K_inter = interfcprop.K_interface;

   for( std::size_t clist_index = 0; clist_index < CList.Size; ++clist_index )
   {
      left_gdbk =  CList[clist_index].left;
      right_gdbk = CList[clist_index].right;
      
      // oil phase & water phase ( assume concurrent flow ONLY )
      std::size_t PoIndex_L = GetPoIndex( left_gdbk );
      std::size_t SwIndex_L = GetSwIndex( left_gdbk );
      std::size_t PoIndex_R = GetPoIndex( right_gdbk );
      std::size_t SwIndex_R = GetSwIndex( right_gdbk );
      
      const CentralProperty& prop = GetUpstream( left_gdbk,
						 mUnew[ PoIndex_L ],
						 right_gdbk,
						 mUnew[ PoIndex_R ] );

      // pp( mUnew[ PoIndex_L ] );
      // pp( mUnew[ PoIndex_R ] );

      // pp( vec_centralprop[left_gdbk].density_oil );
      // pp( vec_centralprop[right_gdbk].density_oil );
	  

      // oil
      ADscalar potential_oil = PROPERTY::Potential( mUnew[ PoIndex_R ], mUnew[ PoIndex_L ] );

      ADscalar trans_oil = PROPERTY::Transmissibility( K_inter, trans_constant,
						       prop.kro, prop.density_oil,
						       prop.viscosity_oil  );

						       
      // water
      ADscalar potential_wat = PROPERTY::Potential( mUnew[ PoIndex_R ], mUnew[ PoIndex_L ] );

      ADscalar trans_wat = PROPERTY::Transmissibility( K_inter, trans_constant,
						       prop.krw, prop.density_wat,
						       prop.viscosity_wat );

      // flux
      ADscalar flux_oil = Flux_Term ( trans_oil, potential_oil ) / DV;
      ADscalar flux_wat = Flux_Term ( trans_wat, potential_wat ) / DV;


      // p( K_inter );
      // p( trans_constant );

      // pp( prop.density_oil );
      // pp( prop.viscosity_oil );


      // printf( "Kro = %.6e\n", prop.kro.value() );
      // printf( "density_oil = %.6e\n", prop.density_oil.value() );
      // printf( "potential_oil = %.6e\n", potential_oil.value() );
      // printf( "flux_oil = %.6e\n", flux_oil.value() );
      // printf( "flux_wat = %.6e\n", flux_wat.value() );
      // n();
      // printf( "trans_oil = %.6e\n", trans_oil.value() );
      // printf( "potential_oil = %.6e\n", potential_oil.value() );
      // printf( "flux_oil = %.6e\n", flux_oil.value() );

      // n();
      // printf( "flux_oil at face# %d = %.6e\n", clist_index, flux_oil.value() );
      // printf( "flux_wat at face# %d = %.6e\n", clist_index, flux_wat.value() );

      // if( clist_index == 0 )
      // {
      // 	 printf( "%-30s%d = %20.12e\n",
      // 		 "residual before at cell# ",
      // 		 left_gdbk, residual[ PoIndex_L ].value() );

      // 	 printf( "%-30s%d = %20.12e\n",
      // 		 "flux_oil at face# ",
      // 		 clist_index, flux_oil.value() );
      // }

      
      residual[ PoIndex_L ] -= flux_oil;
      residual[ PoIndex_R ] += flux_oil;
      residual[ SwIndex_L ] -= flux_wat;
      residual[ SwIndex_R ] += flux_wat;

      
      // if( clist_index == 0 )
      // {
      // 	 printf( "%-30s%d = %14.6e\n",
      // 	      "residual after  at cell# ",
      // 		 left_gdbk, residual[ PoIndex_L ].value() );
      // }
   }
}


// ================================================================== Well
// ================================================================== Well
// ================================================================== Well
//
ADscalar DiscreteProblem :: BHP_Residual ( const Injector& _injector )
{
   const std::size_t i = _injector.I();
   const std::size_t j = _injector.J();
   const std::size_t k = _injector.K();

   std::size_t BlockIndex = grid.GetGridBlockIndex( i, j, k );
   std::size_t PoIndex    = GetPoIndex( BlockIndex );
   std::size_t SwIndex    = GetSwIndex( BlockIndex );
   std::size_t BHPIndex   = GetBHPIndex( _injector.Number() ); 
   const CentralProperty& prop = vec_centralprop[ BlockIndex ];

   // for water injector, make Krw = Krw( 1-Sor )
   ADscalar Sw( 0.0 );
   Sw = mUnew[ SwIndex ];

   ADscalar krw_inj( 0.0 );
   // krw_inj = PROPERTY :: Krw( 1 );
   krw_inj = PROPERTY :: Krw( Sw );

   ADscalar ret( 0.0 );
   ret = _injector.WI() * prop.density_wat * krw_inj / prop.viscosity_wat
     * ( mUnew[ PoIndex ] - mUnew[ BHPIndex ] ) - _injector.Rate();

     
   return ret;
}
//
ADscalar DiscreteProblem :: Sink_Term ( const Producer& _producer, int _tag )
{
   // oil rate ( _tag == 0 )
   // wat rate ( _tag == 1 )
   std::size_t i = _producer.I();
   std::size_t j = _producer.J();
   std::size_t k = _producer.K();

   std::size_t BlockIndex = grid.GetGridBlockIndex( i, j, k );
   std::size_t PoIndex    = GetPoIndex( BlockIndex );
   const CentralProperty& prop = vec_centralprop[ BlockIndex ];

   ADscalar ret( 0.0 );
   if( _tag == 0 )
   {
      ret = _producer.WI() * prop.density_oil * prop.kro / prop.viscosity_oil 
	 * ( mUnew[ PoIndex ] - _producer.BHP() );
   }
   else
   {
      ret = _producer.WI() * prop.density_wat * prop.krw / prop.viscosity_wat
	 * ( mUnew[ PoIndex ] - _producer.BHP() );
   }
   return ret;
}
// --- ---
void DiscreteProblem :: Evaluate_Source_Sink ()
{
   const double DV = grid.Dx * grid.Dy * grid.Dz;
   
   // source: water injection @ constant rate
   for( int i = 0; i < VarNumber_Added; ++i )
   {
      // 1. Add constant rate to residual equation
      const Injector& well = vec_injector[ i ];
      const std::size_t SwIndex = GetSwIndex( well.I(), well.J(), well.K() );

      ADscalar rate( 0.0 );
      rate = well.Rate() / DV;
      residual[ SwIndex ] += rate;

      // 2. Form new residual equation
      const std::size_t BHPIndex = GetBHPIndex( i );
      residual[ BHPIndex ] = BHP_Residual( vec_injector[i] ) / DV;
   }

   // sink: producing oil and water @ constant bhp
   for( std::size_t i = 0; i < vec_producer.size(); ++i )
   {
      const Producer& well = vec_producer[ i ];
      const std::size_t BlockIndex = grid.GetGridBlockIndex( well.I(), well.J(), well.K() );
      const std::size_t PoIndex = GetPoIndex( BlockIndex );
      const std::size_t SwIndex = GetSwIndex( BlockIndex );
      // oil
      ADscalar sink_oil( 0.0 );
      sink_oil = Sink_Term( well, 0 ) / DV;
      residual[ PoIndex ] += sink_oil;

      
      // water
      ADscalar sink_wat( 0.0 );
      sink_wat = Sink_Term( well, 1 ) / DV;
      residual[ SwIndex ] += sink_wat;

      // std::cout << "WI: "<< well.WI() << std::endl;
      // std::cout << "Re: "<< well.Re() << std::endl;
      // std::cout << "Producer Oil: " << Sink_Term( well, 0 ) << std::endl;
      // std::cout << "Producer2: " << Sink_Term( well, 1 ) << std::endl;
   }
}


// Auxillary Fucntions
void DiscreteProblem :: ClearCentralProperty ()
{
   for( std::size_t i = 0; i < vec_centralprop.size(); ++i )
   {
      vec_centralprop[i].Reset();
   }
}
   
const std::size_t DiscreteProblem :: GetPoIndex( std::size_t _GridBlockIndex ) const
{
   return _GridBlockIndex << 1;
}

const std::size_t DiscreteProblem :: GetSwIndex( std::size_t _GridBlockIndex ) const
{
   return GetPoIndex( _GridBlockIndex ) + 1;
}

const std::size_t DiscreteProblem :: GetPoIndex( std::size_t _i, std::size_t _j, std::size_t _k ) const 
{
   return grid.GetGridBlockIndex( _i, _j, _k ) << 1;
}

const std::size_t DiscreteProblem :: GetSwIndex( std::size_t _i, std::size_t _j, std::size_t _k ) const
{
   return GetPoIndex( _i, _j, _k ) + 1;
}

const std::size_t DiscreteProblem :: GetGBIndex( std::size_t _i, std::size_t _j, std::size_t _k ) const
{
   return grid.GetGridBlockIndex( _i, _j, _k );
}

const std::size_t DiscreteProblem :: GetBHPIndex( int _injector_number ) const
{
  return _injector_number + VarNumber_Base;
}


int DiscreteProblem :: GetKxColIndex( std::size_t _Kx_GBIndex ) const
{
   return _Kx_GBIndex;
}

int DiscreteProblem :: GetKyColIndex( std::size_t _Ky_GBIndex ) const
{
   return _Ky_GBIndex + GridNumber;
}

int DiscreteProblem :: GetKzColIndex( std::size_t _Kz_GBIndex ) const
{
   return _Kz_GBIndex + 2 * GridNumber;
}

int DiscreteProblem :: GetPSColIndex( std::size_t _ps_GBIndex ) const // porosity
{
   return _ps_GBIndex + 3 * GridNumber;
}
