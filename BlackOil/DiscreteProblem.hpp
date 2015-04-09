#pragma once

#include "FluidRockProperty.hpp"
#include "Gridding.hpp"
#include "Well.hpp"

struct CentralProperty
{
   CentralProperty () : kro(0.0), krw(0.0),
			bo(0.0), bw(0.0),
			pw(0.0),
			spec_weight_oil(0.0), spec_weight_wat(0.0),
			viscosity_oil(0.0), viscosity_wat(0.0),
			porosity(0.0)
   {
   }
   
   void Reset ()
   {
      kro = 0.0;
      krw = 0.0;
      bo = 0.0;
      bw = 0.0;
      pw = 0.0;
      spec_weight_oil = 0.0;
      spec_weight_wat = 0.0;
      viscosity_oil = 0.0;
      viscosity_wat = 0.0;
      porosity = 0.0;
   }
   
   ADscalar kro, krw;
   ADscalar bo, bw;
   ADscalar pw;
   ADscalar spec_weight_oil, spec_weight_wat;
   ADscalar viscosity_oil, viscosity_wat;
   ADscalar porosity;
};

struct InterfcProperty
{			
   double trans_constant_x;
   double trans_constant_y;
   double trans_constant_z;
   std::vector<double> K_interface;
};

class DiscreteProblem
{
public:   
   DiscreteProblem();

private: // for initialization & preparation
   void Activate_Uold ();
   void Initialize_K_interface (); 
   void Add_Producer ( double _rate, std::size_t _i, std::size_t j, std::size_t k, double _rw = 0.5 );
   void Add_Injector ( double _bhp,  std::size_t _i, std::size_t j, std::size_t k, double _rw = 0.5 );
   

public:  // for update
   void Update_Uold ();                                                     // 1
   void Update_Unew ( const std::vector<double>& _negative_newton_update ); // 2
   void Update_Property    ();                                              // 3
   void Update_Accum_Old ( double _time_step );                             // 4

public:  // for evaluation
   
   template< typename VECTOR, typename CSR >
   void Evaluate ( double _time_step,
		   CSR& _Jacobian, VECTOR& _vec_residual );
   
private: // for evaluation
   // Accumulation
   double   Accum_Constant ( double _time_step,
			     double _del_x, double _del_y, double _del_z ) const;

   ADscalar Accum_Term     ( double _accum_constant,
			     const ADscalar& _porosity,
			     const ADscalar& _saturation,
			     const ADscalar& _volume_factor ) const;

   void     Evaluate_Accumulation ( double _time_step );

   // Flux
   const CentralProperty& GetUpstream ( std::size_t _gdbk_left,
					const ADscalar& _pressure_left,
					std::size_t _gdbk_right,
					const ADscalar& _pressure_right ) const;

   ADscalar Flux_Term ( const ADscalar& _transmissibility, const ADscalar& _potential) const;
   void     Evaluate_Flux ();  
 
   // Source & Sink

   ADscalar Sink_Term ( const Producer& _producer, int _tag ); // for constant BHP wells ONLY
   void     Evaluate_Source_Sink ();

public:  // for access & information
   int GetGridNumber       () const { return GridNumber;       }
   int GetVarNumber        () const { return VarNumber; }
   int GetInjectorNumber   () const { return vec_injector.size();   }

   const ADvector& GetUold () const  { return mUold; }
   const ADvector& GetUnew () const  { return mUnew; }

private: // to help functionality
   void ClearCentralProperty ();

   std::size_t GetPoIndex( std::size_t _GridBlockIndex )                    const;
   std::size_t GetSwIndex( std::size_t _GridBlockIndex )                    const;
   std::size_t GetPoIndex( std::size_t _i, std::size_t _j, std::size_t _k ) const;
   std::size_t GetSwIndex( std::size_t _i, std::size_t _j, std::size_t _k ) const;

private:
   Cartesian3D                    grid;
   const int                      GridNumber;
   const int                      VarNumber;
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

};


// ================================================================================
// ================================================================================
// ================================================================================
// ================================================================================
// ================================================================================
DiscreteProblem :: DiscreteProblem () : grid( 200, 200, 25, 20, 30, 3 ), //20-30-3
					GridNumber      ( grid.Nx * grid.Ny * grid.Nz ),
					VarNumber_Base  ( 2 * GridNumber ),
					CList           ( grid ),
					vec_centralprop ( GridNumber, CentralProperty() ),
					interfcprop     (),
					vec_accum_old   ( VarNumber_Base, 0.0 ),
					mUnew           ( VarNumber_Base, 0.0 ),
					mUold           ( VarNumber_Base, 0.0 ),
					residual        ( VarNumber_Base, 0.0 )
{
   interfcprop.trans_constant_x = PROPERTY::Trans_Constant( grid.Dy, grid.Dz, grid.Dx );
   interfcprop.trans_constant_y = PROPERTY::Trans_Constant( grid.Nx, grid.Nz, grid.Ny );
   interfcprop.trans_constant_z = PROPERTY::Trans_Constant( grid.Nx, grid.Ny, grid.Nz );
 
   constprop.Read_From_File( "3_layer_reservoir_properties/layer_1.txt" );
   constprop.Read_From_File( "3_layer_reservoir_properties/layer_2.txt" );
   constprop.Read_From_File( "3_layer_reservoir_properties/layer_3.txt" );

   Activate_Uold();
 
   // inject water with constant rate
   Add_Injector( 100, 6,  6,  2 );
   Add_Injector( 100, 17, 26, 2 );

   // producer with constant bottom pressure
   Add_Producer( 800, 5,  4,  0 );
   Add_Producer( 800, 16, 5,  0 );
   Add_Producer( 800, 6,  23, 0 );
   Add_Producer( 800, 5,  21, 0 );

   mUnew = mUold;
   
   interfcprop.K_interface.resize( CList.Size );
   Initialize_K_interface();
   
   Update_Property();
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
      
      mUold[ SwIndex ].value() = constprop.Siw.value() + 0.25;
      mUold[ SwIndex ].make_independent( SwIndex );
   }
   // 2. Initialize BHP as unknowns
   for( int i = 0; i < GetInjectorNumber(); ++i )
   {
      int InjectorIndex = i + VarNumber;
      mUold[ InjectorIndex ].value() = constprop.Poi.value();
      mUold[ InjectorIndex ].make_independent( InjectorIndex );
   }
}
// --- ---
void DiscreteProblem :: Initialize_K_interface()
{
   std::size_t count = 0;
   for( ; count < CList.XSize; ++count )
   {
      const std::size_t left_index = CList[count].left;
      const std::size_t right_index = CList[count].right;
      interfcprop.K_interface[count] = PROPERTY::K_Interface( constprop.vec_Kx[ left_index ],
							      constprop.vec_Kx[ right_index ] );
   }
   
   for( ; count < CList.XSize + CList.YSize; ++count )
   {
      const std::size_t left_index = CList[count].left;
      const std::size_t right_index = CList[count].right;
      interfcprop.K_interface[count] = PROPERTY::K_Interface( constprop.vec_Ky[ left_index ],
							      constprop.vec_Ky[ right_index ] );
   }
   
   for( ; count < CList.Size; ++count )
   {
      const std::size_t left_index = CList[count].left;
      const std::size_t right_index = CList[count].right;
      interfcprop.K_interface[count] = PROPERTY::K_Interface( constprop.vec_Kz[ left_index ],
							      constprop.vec_Kz[ right_index ] );
   }
}
// --- ---
void DiscreteProblem :: Add_Injector ( double _rate,
				       std::size_t _i, std::size_t _j, std::size_t _k, double _rw )
{
   const std::size_t BlockIndex = grid.GetGridBlockIndex( _i, _j, _k );
   Injector well( _rate, _i, _j, _k, _rw,
		  grid.Dx, grid.Dy,
		  constprop.vec_Kx[ BlockIndex ], constprop.vec_Ky[ BlockIndex ] );

   vec_injector.push_back( well );
   // expand
   ADscalar tmp( 0.0 );
   mUold.push_back( tmp );
   residual.push_back( tmp );
}
//
void DiscreteProblem :: Add_Producer  ( double _bhp,
					std::size_t _i, std::size_t _j, std::size_t _k, double _rw )
{
   std::size_t BlockIndex = grid.GetGridBlockIndex( _i, _j, _k );
   Producer well( _bhp, _i, _j, _k, _rw,
		  grid.Dx, grid.Dy,
		  constprop.vec_Kx[ BlockIndex ], constprop.vec_Ky[ BlockIndex ] );
   
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
   for( std::size_t i = 0; i < VarNumber; )
   {
      // pressure
      mUnew[i].value() -= _negative_newton_update[i];
      ++i;
      // water saturation
      mUnew[i].value() -= _negative_newton_update[i];
      ++i;
   }
}
//3
void DiscreteProblem :: Update_Property ()
{
   for( std::size_t i = 0; i < GridNumber; ++i )
   {
      CentralProperty&  prop = vec_centralprop[ i ];
   
      const std::size_t PoIndex = GetPoIndex( i );
      const std::size_t SwIndex = GetSwIndex( i );
      const ADscalar&  _Po = mUnew[ PoIndex ];
      const ADscalar&  _Sw = mUnew[ SwIndex ];

      ADscalar S_wd( 0.0 );
      S_wd = PROPERTY::Swd( _Sw, constprop.Siw, constprop.Sor );

      prop.kro = PROPERTY::Kro( S_wd );
      prop.krw = PROPERTY::Krw( S_wd );
	    
      prop.pw = PROPERTY::Pw( _Po );

      prop.bo = PROPERTY::Bo ( _Po, constprop.Bob, constprop.Co, constprop.Pb );
      prop.bw = PROPERTY::Bw ( prop.pw, constprop.Bwb, constprop.Cw, constprop.Pb );
   
      ADscalar dens_o = PROPERTY::Density_Oil( prop.bo, constprop.DensOil_sc );
      ADscalar dens_w = PROPERTY::Density_Wat( prop.bw, constprop.DensWat_sc );
   
      prop.spec_weight_oil = PROPERTY::SpecificWeight( dens_o );
      prop.spec_weight_wat = PROPERTY::SpecificWeight( dens_w );
   
      prop.viscosity_oil = constprop.VisO;
      prop.viscosity_wat = constprop.VisW;
   
      const double porosity_b = constprop.vec_porosity_b[ i ];
      prop.porosity = PROPERTY::Porosity( _Po, porosity_b, constprop.Cr, constprop.Pb );

      // std::cout << "---------------------------------------- "<< _grid_number << std::endl;
      // std::cout << "****************"<< std::endl;
      // std::cout << "Po =\t" << _Po.value() << std::endl;
      // std::cout << "Sw =\t" << _Sw.value() << std::endl;
      // std::cout << "Pw =\t" << prop.pw.value() << std::endl;   
      // std::cout << "****************"<< std::endl;
      // std::cout << "Kro =\t" << prop.kro.value() << std::endl;
      // std::cout << "Krw =\t" << prop.krw.value() << std::endl;
      // std::cout << "Bo =\t" << prop.bo.value() << std::endl;
      // std::cout << "Bw =\t" << prop.bw.value() << std::endl;
      // std::cout << "spwto =\t" << prop.spec_weight_oil.value() << std::endl;
      // std::cout << "spwtw =\t" << prop.spec_weight_wat.value() << std::endl;
      // std::cout << "viso =\t" << prop.viscosity_oil.value() << std::endl;
      // std::cout << "visw =\t" << prop.viscosity_wat.value() << std::endl;
      // std::cout << "poro_b =\t" << porosity_b << std::endl;
      // std::cout << "poro =\t" << prop.porosity.value() << std::endl;
   }
}
// 4
void DiscreteProblem :: Update_Accum_Old ( double _time_step )
{
   std::cout << "Accum_Old Term Updated..." << std::endl;
   const double accum_const = Accum_Constant( _time_step, grid.Dx, grid.Dy, grid.Dz );
   for( int i = 0; i < GetGridNumber(); ++i )
   {      
      const CentralProperty& prop = vec_centralprop[ i ];
      const std::size_t PoIndex = GetPoIndex( i );
      const std::size_t SwIndex = GetSwIndex( i );

      // oil phase
      ADscalar So( 0.0 );
      So = 1.0 - mUnew[ SwIndex ];
      vec_accum_old[ PoIndex ] = accum_const * prop.porosity.value() * So.value() / prop.bo.value();

      // water phase
      const ADscalar& Sw = mUnew[ SwIndex ];
      vec_accum_old[ SwIndex ] = accum_const * prop.porosity.value() * Sw.value() / prop.bw.value();
   }
}




template< typename VECTOR, typename CSR >
void DiscreteProblem :: Evaluate ( double _time_step,
				   CSR& _Jacobian, VECTOR& _vec_residual )
{
   residual.resize( VarNumber, 0.0 );
   Evaluate_Accumulation( _time_step );
   Evaluate_Flux();
   Evaluate_Source_Sink();

   //std::cout << residual << std::endl;
   
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
   return _del_x * _del_y * _del_z / ( 5.615 * _time_step );
}
// --- ---
ADscalar DiscreteProblem :: Accum_Term ( double _accum_constant,
					 const ADscalar& _porosity,
					 const ADscalar& _saturation,
					 const ADscalar& _volume_factor ) const
{
   ADscalar ret( 0.0 );
   ret = _accum_constant * ( _porosity * _saturation / _volume_factor );
   return ret;
}

void DiscreteProblem :: Evaluate_Accumulation ( double _time_step )
{
   const double accum_constant = Accum_Constant( _time_step, grid.Dx, grid.Dy, grid.Dz );

   for( std::size_t i = 0; i < GridNumber; ++i )
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
      residual[ PoIndex ] -= vec_accum_old[ PoIndex ];

      // 2. update oil accumulation term, which will become the "old" one in the next iterate
      accum_oil = Accum_Term( accum_constant,
      			      prop.porosity,
      			      1.0 - mUnew[SwIndex],
      			      prop.bo );

      // 3.
      residual[ PoIndex ] += accum_oil;

      //std::cout << "accum_oil: " << accum_oil << std::endl;

      // water phase
      // 1.
      residual[ SwIndex ] -= vec_accum_old[ SwIndex ];
      // 2.
      accum_wat = Accum_Term( accum_constant,
      			      prop.porosity,
      			      mUnew[SwIndex],
      			      prop.bw );
      // 3.
      residual[ SwIndex ] += accum_wat;

      //std::cout << "accum_wat: " << accum_wat << std::endl;
      
      // std::cout << "--------------------------------------- " << i << std::endl;
      // std::cout << residual[PoIndex].value() << " --- --- " << residual[SwIndex].value() << std::endl;
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
    * I'm lazy. I'll take left in case of equality.
    */
   if( _pressure_left >= _pressure_right )
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
   std::size_t left_gdbk, right_gdbk;
   double trans_constant;
   double del_depth;
   for( std::size_t clist_index = 0; clist_index < CList.Size; ++clist_index )
   {
      // direction-dependent properties
      if( clist_index < CList.XSize )
      {
	 trans_constant = interfcprop.trans_constant_x;
	 del_depth = 0.0;
      }
      else if( clist_index >= CList.XSize && clist_index < CList.XSize + CList.YSize )
      {
	 trans_constant = interfcprop.trans_constant_y;
	 del_depth = 0.0;
      }
      else
      {
	 trans_constant = interfcprop.trans_constant_z;
	 del_depth = grid.Dz;
      }
	 
      left_gdbk =  CList[clist_index].left;
      right_gdbk = CList[clist_index].right;

      double K_inter = interfcprop.K_interface[ clist_index ];
      
      // oil phase & water phase ( assume concurrent flow ONLY )
      std::size_t PoIndex_L = GetPoIndex( left_gdbk );
      std::size_t SwIndex_L = GetSwIndex( left_gdbk );
      std::size_t PoIndex_R = GetPoIndex( right_gdbk );
      std::size_t SwIndex_R = GetSwIndex( right_gdbk );
      
      const CentralProperty& prop = GetUpstream( left_gdbk,
						 mUnew[ PoIndex_L ],
						 right_gdbk,
						 mUnew[ PoIndex_R ] );

      // oil
      ADscalar potential_oil = PROPERTY::Potential   ( mUnew[ PoIndex_R ],
						       mUnew[ PoIndex_L ],
						       prop.spec_weight_oil,
						       del_depth );
      ADscalar trans_oil = PROPERTY::Transmissibility( K_inter,
						       trans_constant,
						       prop.kro,
						       prop.viscosity_oil,
						       prop.bo );
      // water
      ADscalar potential_wat = PROPERTY::Potential   ( vec_centralprop[ right_gdbk ].pw,
						       vec_centralprop[ left_gdbk  ].pw,
						       prop.spec_weight_wat,
						       del_depth );	 
      ADscalar trans_wat = PROPERTY::Transmissibility( K_inter,
						       trans_constant,
						       prop.krw,
						       prop.viscosity_wat,
						       prop.bw );
      
      // flux
      ADscalar flux_oil = Flux_Term ( trans_oil, potential_oil );
      ADscalar flux_wat = Flux_Term ( trans_wat, potential_wat );

      // std::cout << "----------------------------------- "<< left_gdbk << " --- " << right_gdbk << std::endl;
      // std::cout << "potential_water: "<< potential_wat << std::endl;
      // std::cout << "transmiss_water: "<< trans_wat << std::endl;
      //std::cout << "flux_water: "     << flux_wat << std::endl;
      //std::cout << "potential_oil: "<< potential_oil << std::endl;
      //std::cout << "transmiss_oil: "<< trans_oil << std::endl;
      //std::cout << "flux_oil: "     << flux_oil << std::endl;
      
      //std::cout << prop.viscosity_wat << std::endl;
      //std::cout << prop.krw << std::endl;      
      //std::cout << trans_oil << std::endl;

      
      residual[ PoIndex_L ] -= flux_oil;
      residual[ PoIndex_R ] += flux_oil;
      residual[ SwIndex_L ] -= flux_wat;
      residual[ SwIndex_R ] += flux_wat;

      // ADscalar left_block_oil = -flux_oil;
      // ADscalar left_block_wat = -flux_wat;

      // std::cout << "Left Oil: "<< left_block_oil << std::endl;
      // std::cout << "Right Oil: "<< flux_oil << std::endl;

   }
}


// ================================================================== Well
// ================================================================== Well
// ================================================================== Well
//
ADscalar DiscreteProblem :: Source_Term ( const Injector& _injector )
{

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
      ret = _producer.WI() * prop.kro / ( prop.viscosity_oil * prop.bo )
	 * ( mUnew[ PoIndex ] - _producer.BHP() );
   }
   else
   {
      ret = _producer.WI() * prop.krw / ( prop.viscosity_wat * prop.bw )
	 * ( mUnew[ PoIndex ] - _producer.BHP() );
   }
   return ret;
}
// --- ---
void DiscreteProblem :: Evaluate_Source_Sink ()
{
   // source: water injection @ constant rate
   for( std::size_t i = 0; i < vec_injector.size(); ++i )
   {
      const Injector& well = vec_injector[ i ];
      const std::size_t SwIndex = GetSwIndex( well.I(), well.J(), well.K() );      
      
      residual[ SwIndex ] -= well.Rate();
   }
   // sink: producing oil and water @ constant bhp
   for( std::size_t i = 0; i < vec_producer.size(); ++i )
   {
      const Producer& well = vec_producer[ i ];
      const std::size_t BlockIndex = grid.GetGridBlockIndex( well.I(), well.J(), well.K() );
      const std::size_t PoIndex = GetPoIndex( BlockIndex );
      const std::size_t SwIndex = GetSwIndex( BlockIndex );
      // oil
      residual[ PoIndex ] += Sink_Term( well, 0 );
      // water
      residual[ SwIndex ] += Sink_Term( well, 1 );

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
   
std::size_t DiscreteProblem :: GetPoIndex( std::size_t _GridBlockIndex ) const
{
   return _GridBlockIndex << 1;
}

std::size_t DiscreteProblem :: GetSwIndex( std::size_t _GridBlockIndex ) const
{
   return GetPoIndex( _GridBlockIndex ) + 1;
}

std::size_t DiscreteProblem :: GetPoIndex( std::size_t _i, std::size_t _j, std::size_t _k ) const 
{
   return grid.GetGridBlockIndex( _i, _j, _k ) << 1;
}

std::size_t DiscreteProblem :: GetSwIndex( std::size_t _i, std::size_t _j, std::size_t _k ) const
{
   return GetPoIndex( _i, _j, _k ) + 1;
}
