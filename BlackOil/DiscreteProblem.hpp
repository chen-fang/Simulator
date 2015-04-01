#pragma once

#include "FluidRockProperty.hpp"
#include "Gridding.hpp"
#include "Well.hpp"

struct CentralProperty
{
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
   void Update_AllProperty    ();
   void Evaluate ( double _time_step );

private:
   /* ********* Initiates the discretized problem ********* */
   void     Initialize_Unknown ();
   void     Initialize_K_interface ();

   // Add a well with const BHP if tag == 0. Else, add a well with constant flow rate.
   void     Add_Well ( int _tag, double _value,
		       std::size_t _i, std::size_t _j, std::size_t _k, double _rw );
   //void     Shut_Well( std::size_t _i, std::size_t _j, std::size_t _k );
   //void     Change_Well_Mode ();


   /* ********* Calculate properties according to unknown variables ********* */
   void     Update_GridProperty   ( std::size_t _grid_number );



   /* ********* Accumulation ********* */
   //
   double   Accum_Constant ( double _time_step,
			     double _del_x, double _del_y, double _del_z ) const;
   //
   ADscalar Accum_Term     ( double _accum_constant,
			     const ADscalar& _porosity,
			     const ADscalar& _saturation,
			     const ADscalar& _volume_factor ) const;
   //
   void     Evaluate_Accumulation ( double _time_step );


   /* ********* Flux ********* */
   //
   const CentralProperty& GetUpstream ( std::size_t _gdbk_left,
					const ADscalar& _pressure_left,
					std::size_t _gdbk_right,
					const ADscalar& _pressure_right ) const;
   //
   ADscalar Flux_Term ( const ADscalar& _transmissibility, const ADscalar& _potential) const;
   //
   void     Evaluate_Flux ();  
 
   /* ********* Well ********* */
   //
   ADscalar Well_Term ( double _WI, const ADscalar& _krm,
			const ADscalar& _viscosity_m, const ADscalar& _bm,
			const ADscalar& _Po, double _Pwf );
   //
   void     Evaluate_Well ();

   
   /* ********* Auxillary Functions ********* */
   //
   void ClearCentralProperty ();
   //
   std::size_t GetPoIndex( std::size_t _GridBlockIndex )                    const;
   std::size_t GetSwIndex( std::size_t _GridBlockIndex )                    const;
   std::size_t GetPoIndex( std::size_t _i, std::size_t _j, std::size_t _k ) const;
   std::size_t GetSwIndex( std::size_t _i, std::size_t _j, std::size_t _k ) const;

private:
   Cartesian3D                    grid;
   ConnectionList                 CList;
   
   ADvector                       vec_accum_oil;
   ADvector                       vec_accum_wat;
   
   std::vector< CentralProperty > vec_centralprop;
   InterfcProperty                interfcprop;

   std::vector< Well_ConstRate >  vec_const_rate_well;
   std::vector< Well_ConstBHP >   vec_const_bhp_well;
		       
   ADvector vec_unknown; // water saturation & oil pressure
   ADvector residual;
};


// ================================================================================
// ================================================================================
// ================================================================================
// ================================================================================
// ================================================================================
DiscreteProblem :: DiscreteProblem ()
{
   grid.Set( 200, 200, 25, 20, 20, 3 );
   CList.Initialize( grid );
   interfcprop.trans_constant_x = PROPERTY::Trans_Constant( grid.Ny, grid.Nz, grid.Nx );
   interfcprop.trans_constant_y = PROPERTY::Trans_Constant( grid.Nx, grid.Nz, grid.Ny );
   interfcprop.trans_constant_z = PROPERTY::Trans_Constant( grid.Nx, grid.Ny, grid.Nz );
   
   PROPERTY :: Read_Layer_Property     ( "3_layer_reservoir_properties/layer_1.txt" );
   PROPERTY :: Read_Layer_Property     ( "3_layer_reservoir_properties/layer_2.txt" );
   PROPERTY :: Read_Layer_Property     ( "3_layer_reservoir_properties/layer_3.txt" );

   const std::size_t Grid_Number = grid.Nx * grid.Ny * grid.Nz;
   const std::size_t Var_Number = 2 * Grid_Number;

   vec_accum_oil.resize      ( Var_Number );
   vec_accum_wat.resize  ( Var_Number );
   
   vec_centralprop.resize( Grid_Number );
   vec_unknown.resize    ( Var_Number );
   residual.resize       ( Var_Number );

   interfcprop.K_interface.resize( CList.Size );
   Initialize_K_interface();
   
   Initialize_Unknown();

}


void DiscreteProblem :: Evaluate ( double _time_step )
{
   Evaluate_Accumulation( _time_step );
   Evaluate_Flux();
}


void DiscreteProblem :: Initialize_Unknown()
{
   // 1. Activation
   grid.Activate_NaturalOrdering( vec_unknown );
   // 2. Initialize Po and Sw
   std::size_t count_grid_number = 0;
   for( std::size_t k = 0; k < grid.Nz; ++k )
   {
      for( std::size_t j = 0; j < grid.Ny; ++j )
      {
	 for( std::size_t i = 0; i < grid.Nx; ++i )
	 {
	    const std::size_t PoIndex = GetPoIndex( count_grid_number );
	    const std::size_t SwIndex = GetSwIndex( count_grid_number );
	    vec_unknown[ PoIndex ] = PROPERTY :: Poi;
	    vec_unknown[ SwIndex ] = PROPERTY :: Siw;
	    ++count_grid_number;
	 }
      }
   }
}

void DiscreteProblem :: Initialize_K_interface()
{
   std::size_t count = 0;
   for( ; count < CList.XSize; ++count )
   {
      const std::size_t left_index = CList[count].left;
      const std::size_t right_index = CList[count].right;
      interfcprop.K_interface[count] = PROPERTY::K_Interface( PROPERTY::AbsK::vec_Kx[ left_index ],
							      PROPERTY::AbsK::vec_Kx[ right_index ] );
   }
   
   for( ; count < CList.XSize + CList.YSize; ++count )
   {
      const std::size_t left_index = CList[count].left;
      const std::size_t right_index = CList[count].right;
      interfcprop.K_interface[count] = PROPERTY::K_Interface( PROPERTY::AbsK::vec_Ky[ left_index ],
							      PROPERTY::AbsK::vec_Ky[ right_index ] );
   }
   
   for( ; count < CList.Size; ++count )
   {
      const std::size_t left_index = CList[count].left;
      const std::size_t right_index = CList[count].right;
      interfcprop.K_interface[count] = PROPERTY::K_Interface( PROPERTY::AbsK::vec_Kz[ left_index ],
							      PROPERTY::AbsK::vec_Kz[ right_index ] );
   }
}

void DiscreteProblem :: Add_Well( int _tag, double _value,
				  std::size_t _i, std::size_t _j, std::size_t _k, double _rw )
{
   // If tag == 0,
   // -- constant BHP;
   // Else
   // -- constant flow rate.
   std::size_t index = grid.GetGridBlockIndex( _i, _j, _k );
   if ( _tag == 0 )
   {
      Well_ConstBHP well( _value,
			  _i, _j, _k, _rw,
			  grid.Dx, grid.Dy,
			  PROPERTY::AbsK::vec_Kx[ index ],
			  PROPERTY::AbsK::vec_Ky[ index ] );
      vec_const_bhp_well.push_back( well );
   }
   else
   {
      Well_ConstRate well( _tag, _value,
			   _i, _j, _k, _rw,
			   grid.Dx, grid.Dy,
			   PROPERTY::AbsK::vec_Kx[ index ],
			   PROPERTY::AbsK::vec_Ky[ index ] );
      vec_const_rate_well.push_back( well );
   }
}

void DiscreteProblem :: Update_GridProperty( std::size_t _grid_number )
{
   vec_centralprop.clear();
   CentralProperty& prop = vec_centralprop[ _grid_number ];

   const std::size_t PoIndex = GetPoIndex( _grid_number );
   const std::size_t SwIndex = GetSwIndex( _grid_number );
   const ADscalar& _Po = vec_unknown[ PoIndex ];
   const ADscalar& _Sw = vec_unknown[ SwIndex ];

   ADscalar S_wd = PROPERTY::Swd( _Sw );
   prop.kro = PROPERTY::Kro( S_wd );
   prop.krw = PROPERTY::Krw( S_wd );
	    
   prop.pw = PROPERTY::Pw( _Po );

   prop.bo = PROPERTY::Bo ( _Po );
   prop.bw = PROPERTY::Bw ( prop.pw );
   
   ADscalar dens_o = PROPERTY::Density_Oil( prop.bo );
   ADscalar dens_w = PROPERTY::Density_Wat( prop.bw );
   
   prop.spec_weight_oil = PROPERTY::SpecificWeight( dens_o );
   prop.spec_weight_wat = PROPERTY::SpecificWeight( dens_w );
   
   prop.viscosity_oil = PROPERTY::VisO;
   prop.viscosity_wat = PROPERTY::VisW;
   
   double porosity_b = PROPERTY::vec_porosity_b[ _grid_number ];
   prop.porosity    = PROPERTY::Porosity( _Po, porosity_b );

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
   // std::cout << "poro =\t" << prop.porosity.value() << std::endl;
}

void DiscreteProblem :: Update_AllProperty ()
{
   std::size_t count_grid_number = 0;
   for( std::size_t k = 0; k < grid.Nz; ++k )
   {
      for( std::size_t j = 0; j < grid.Ny; ++j )
      {
	 for( std::size_t i = 0; i < grid.Nx; ++i )
	 {
	    Update_GridProperty( count_grid_number );
	    ++count_grid_number;
	 }
      }
   }
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
   return _accum_constant * ( _porosity * _saturation / _volume_factor );
}
// --- ---
void DiscreteProblem :: Evaluate_Accumulation ( double _time_step )
{
   const double accum_constant = Accum_Constant( _time_step, grid.Dx, grid.Dy, grid.Dz );
   std::size_t count_grid_number = 0;
   for( std::size_t k = 0; k < grid.Nz; ++k )
   {
      for( std::size_t j = 0; j < grid.Ny; ++j )
      {
	 for( std::size_t i = 0; i < grid.Nx; ++i )
	 {
	    // evaluate
	    Update_GridProperty( count_grid_number );

	    // accumulate
	    const CentralProperty& prop = vec_centralprop[ count_grid_number ];
	    const std::size_t PoIndex = GetPoIndex( count_grid_number );
	    const std::size_t SwIndex = GetSwIndex( count_grid_number );	    
	    
	    // oil phase
	    ADscalar& accum_oil = vec_accum_oil[ count_grid_number ];
	    // 1. use oil accumulation value;
	    accum_oil.make_constant();
	    residual[ PoIndex ] -= accum_oil;
	    
	    // 2. update oil accumulation term, which will become the "old" one in the next iterate
	    accum_oil = Accum_Term( accum_constant,
				    prop.porosity,
				    1.0 - vec_unknown[SwIndex],
				    prop.bo );
	    // 3.
	    residual[ PoIndex ] += accum_oil;

	    // water phase
	    ADscalar& accum_wat = vec_accum_wat[ count_grid_number ];
	    // 1.
	    accum_wat.make_constant();
	    residual[ SwIndex ] -= accum_wat;
	    // 2.
	    accum_wat = Accum_Term( accum_constant,
				    prop.porosity,
				    vec_unknown[SwIndex],
				    prop.bw );
	    // 3.
	    residual[ SwIndex ] += accum_wat;
	    
	    ++count_grid_number;
	 }
      }
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
   return _transmissibility * _potential;
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

      // std::cout << "------------------------------------- " << clist_index << std::endl;
      // std::cout << left_gdbk << "\t" << right_gdbk << std::endl;
      // std::cout << "K_inter =\t" << K_inter << std::endl;
      
      // oil phase & water phase ( assume concurrent flow ONLY )
      std::size_t PoIndex_L = GetPoIndex( left_gdbk );
      std::size_t SwIndex_L = GetSwIndex( left_gdbk );
      std::size_t PoIndex_R = GetPoIndex( right_gdbk );
      std::size_t SwIndex_R = GetSwIndex( right_gdbk );
      
      const CentralProperty& prop = GetUpstream( left_gdbk,
						 vec_unknown[ PoIndex_L ],
						 right_gdbk,
						 vec_unknown[ PoIndex_R ] );

      // oil
      ADscalar potential_oil = PROPERTY::Potential   ( vec_unknown[ PoIndex_R ],
						       vec_unknown[ PoIndex_L ],
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

      residual[ PoIndex_L ] -= flux_oil;
      residual[ PoIndex_R ] += flux_oil;
      residual[ SwIndex_L ] -= flux_wat;
      residual[ SwIndex_R ] += flux_wat;
   }
}


// ================================================================== Well
// ================================================================== Well
// ================================================================== Well
// --- ---
ADscalar DiscreteProblem :: Well_Term ( double _WI, const ADscalar& _krm,
					const ADscalar& _viscosity_m, const ADscalar& _bm,
					const ADscalar& _Po, double _Pwf )
{
   return _WI * _krm / ( _viscosity_m * _bm ) * ( _Po - _Pwf );
}
// --- ---
void DiscreteProblem :: Evaluate_Well ()
{
   // wells with constant rate
   for( std::size_t i = 0; i < vec_const_rate_well.size(); ++i )
   {
      const Well_ConstRate& well = vec_const_rate_well[ i ];
      const std::size_t gdbk_index = grid.GetGridBlockIndex( well.I(), well.J(), well.K() );
      if( well.Tag() == 0 ) // const oil rate
      {
	 const std::size_t PoIndex = GetPoIndex( gdbk_index );
	 residual[ PoIndex ] += well.Rate();
      }
      else if( well.Tag() == 1 ) // const water rate
      {
	 const std::size_t SwIndex = GetSwIndex( gdbk_index );
	 residual[ SwIndex ] += well.Rate();
      }
      else // const total rate
      {
	 // skip for now
      }
   }
   // wells with constant pressure
   for( std::size_t i = 0; i < vec_const_bhp_well.size(); ++i )
   {
      const Well_ConstBHP& well = vec_const_bhp_well[ i ];
      const std::size_t gdbk_index = grid.GetGridBlockIndex( well.I(), well.J(), well.K() );
      const std::size_t PoIndex = GetPoIndex( gdbk_index );
      const std::size_t SwIndex = GetSwIndex( gdbk_index );
      const CentralProperty& prop = vec_centralprop[ gdbk_index ];
      // oil
      residual[ PoIndex ] += Well_Term( well.WI(), prop.kro,
					prop.viscosity_oil, prop.bo,
					vec_unknown[ PoIndex ], well.BHP() );
      // water
      residual[ SwIndex ] += Well_Term( well.WI(), prop.krw,
					prop.viscosity_wat, prop.bw,
					vec_unknown[ PoIndex ], well.BHP() );
   }
}


// Auxillary Fucntions
void DiscreteProblem :: ClearCentralProperty ()
{
   for( std::size_t i = 0; i < vec_centralprop.size(); ++i )
   {
      vec_centralprop[i].kro             = 0.0;
      vec_centralprop[i].krw             = 0.0;
      vec_centralprop[i].bo              = 0.0;
      vec_centralprop[i].bw              = 0.0;
      vec_centralprop[i].spec_weight_oil = 0.0;
      vec_centralprop[i].spec_weight_wat = 0.0;
      vec_centralprop[i].viscosity_oil   = 0.0;
      vec_centralprop[i].viscosity_wat   = 0.0;
      vec_centralprop[i].porosity        = 0.0;
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
