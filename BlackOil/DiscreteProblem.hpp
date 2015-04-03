#pragma once

#include "FluidRockProperty.hpp"
#include "Gridding.hpp"
#include "Well.hpp"

struct CentralProperty
{
   void Initialize ()
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
   void Initialize ( double _first_time_step );

   void Update_Unknown ( const std::vector<double>& _negative_newton_update );
   void Update_AllProperty    ();

   //void Update_Unknown ( const std::vector<double>& _negative_newton_update );
   
   template< typename VECTOR, typename CSR >
   void Evaluate ( double _time_step,
		   CSR& _Jacobian, VECTOR& _vec_residual );

private:
   /* ********* Initiates the discretized problem ********* */
   void     Initialize_Unknown ();
   void     Initialize_K_interface ();

   // Add a well with const BHP if tag == 0. Else, add a well with constant flow rate.
   void     Add_Well ( int _tag, double _value,
		       std::size_t _i, std::size_t _j, std::size_t _k, double _rw = 0.5 );
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
   //void     Initialize_Accum      ( double _first_timp_step );
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
   std::size_t                    GridNumber;
   std::size_t                    VarNumber;
   ConnectionList                 CList;

   std::vector< double >          vec_accum_oil_old;
   std::vector< double >          vec_accum_wat_old;
   ADvector                       vec_accum_oil;
   ADvector                       vec_accum_wat;

   ConstProperty                  constprop;
   std::vector< CentralProperty > vec_centralprop;
   InterfcProperty                interfcprop;

   std::vector< Well_ConstRate >  vec_const_rate_well;
   std::vector< Well_ConstBHP >   vec_const_bhp_well;
		       
   ADvector vec_unknown; // water saturation & oil pressure
   ADvector residual;

public:

};


// ================================================================================
// ================================================================================
// ================================================================================
// ================================================================================
// ================================================================================
DiscreteProblem :: DiscreteProblem ()
{
   //grid.Set( 200, 200, 25, 20, 30, 3 );
   grid.Set( 200, 200, 25, 1, 1, 1 );
   
   GridNumber = grid.Nx * grid.Ny * grid.Nz;
   VarNumber  = 2 * GridNumber;
   
   CList.Initialize( grid );

   interfcprop.trans_constant_x = PROPERTY::Trans_Constant( grid.Ny, grid.Nz, grid.Nx );
   interfcprop.trans_constant_y = PROPERTY::Trans_Constant( grid.Nx, grid.Nz, grid.Ny );
   interfcprop.trans_constant_z = PROPERTY::Trans_Constant( grid.Nx, grid.Ny, grid.Nz );
 
   constprop.Read_From_File( "3_layer_reservoir_properties/layer_1.txt" );
   constprop.Read_From_File( "3_layer_reservoir_properties/layer_2.txt" );
   constprop.Read_From_File( "3_layer_reservoir_properties/layer_3.txt" );

   vec_accum_oil_old.resize ( GridNumber );
   vec_accum_wat_old.resize ( GridNumber );
   vec_accum_oil.resize     ( GridNumber );
   vec_accum_wat.resize     ( GridNumber );
   
   vec_centralprop.resize( GridNumber );
   for( std::size_t i = 0; i < GridNumber; ++i )
   {
      vec_centralprop[ i ].Initialize();
   }
   
   vec_unknown.resize    ( VarNumber );
   residual.resize       ( VarNumber );

   interfcprop.K_interface.resize( CList.Size );
}

void DiscreteProblem :: Initialize ( double _first_time_step )
{
   Initialize_Unknown();
   Initialize_K_interface();
   Update_AllProperty();
   //Initialize_Accum( _first_time_step );
}
   
template< typename VECTOR, typename CSR >
void DiscreteProblem :: Evaluate ( double _time_step,
				   CSR& _Jacobian, VECTOR& _vec_residual )
{       
   Update_AllProperty();
   
   Evaluate_Accumulation( _time_step );
   //Evaluate_Flux();
   // // Evaluate_Well();

   residual.extract_CSR( _vec_residual,
    			 _Jacobian.Row(), _Jacobian.Col(), _Jacobian.NZV() );
   //Jacobian.Update_Status();
}


void DiscreteProblem :: Initialize_Unknown()
{
   // 1. Activation
   grid.Activate_NaturalOrdering( vec_unknown );
   // 2. Initialize Po and Sw
   for( std::size_t i = 0; i < GridNumber; ++i )
   {
      const std::size_t PoIndex = GetPoIndex( i );
      const std::size_t SwIndex = GetSwIndex( i );
      vec_unknown[ PoIndex ].value() = constprop.Poi.value();
      vec_unknown[ SwIndex ].value() = constprop.Siw.value() + 0.1;
   }
}

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
			  constprop.vec_Kx[ index ],
			  constprop.vec_Ky[ index ] );
      vec_const_bhp_well.push_back( well );
   }
   else
   {
      Well_ConstRate well( _tag, _value,
			   _i, _j, _k, _rw,
			   grid.Dx, grid.Dy,
			   constprop.vec_Kx[ index ],
			   constprop.vec_Ky[ index ] );
      vec_const_rate_well.push_back( well );
   }
}


void DiscreteProblem :: Update_GridProperty( std::size_t _grid_number )
{
   CentralProperty& prop = vec_centralprop[ _grid_number ];
   
   const std::size_t PoIndex = GetPoIndex( _grid_number );
   const std::size_t SwIndex = GetSwIndex( _grid_number );
   const ADscalar& _Po = vec_unknown[ PoIndex ];
   const ADscalar& _Sw = vec_unknown[ SwIndex ];

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
   
   const double porosity_b = constprop.vec_porosity_b[ _grid_number ];
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

void DiscreteProblem :: Update_Unknown ( const std::vector<double>& _negative_newton_update )
{
   for( std::size_t i = 0; i < VarNumber; ++i )
   {
      //std::cout << "---------------------------------- "<< i << std::endl;
      //std::cout << vec_unknown[i].value() << " --- --- ";
      vec_unknown[i].value() -= _negative_newton_update[i];
      //std::cout << vec_unknown[i].value() << std::endl;
   }
}

void DiscreteProblem :: Update_AllProperty ()
{
   for( std::size_t i = 0; i < GridNumber; ++i )
   {
      Update_GridProperty( i );
   }
}

// void Update_Unknown ( const std::vector<double>& _negative_newton_update )
// {
//    //std::size_t GridNumber = grid.GetGridBlockNumber();
//    // for( std::size_t i = 0; i < 1; ++i )
//    // {
//    //    vec_unknown[i].value() -= _negative_newton_update[i];
//    // }
// }
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
// // --- ---
// void DiscreteProblem :: Initialize_Accum ( double _time_step )
// {
//    const double accum_constant = Accum_Constant( _time_step, grid.Dx, grid.Dy, grid.Dz );
//    for( std::size_t i = 0; i < GridNumber; ++i )
//    {
//       const CentralProperty& prop = vec_centralprop[ i ];
//       const std::size_t SwIndex = GetSwIndex( i );	    

//       vec_accum_oil[ i ] = Accum_Term( accum_constant,
// 				       prop.porosity,
// 				       1.0 - vec_unknown[ SwIndex ],
// 				       prop.bo );
//       vec_accum_wat[ i ] = Accum_Term( accum_constant,
// 				       prop.porosity,
// 				       vec_unknown[ SwIndex ],
// 				       prop.bw );
//    }
// }
// --- ---
void DiscreteProblem :: Evaluate_Accumulation ( double _time_step )
{
   const double accum_constant = Accum_Constant( _time_step, grid.Dx, grid.Dy, grid.Dz );
   for( std::size_t i = 0; i < GridNumber; ++i )
   {
      std::cout << " ------------------------------ " << i << std::endl;
      // accumulate
      const CentralProperty& prop = vec_centralprop[ i ];
      const std::size_t PoIndex = GetPoIndex( i );
      const std::size_t SwIndex = GetSwIndex( i );	    

      // oil phase
      ADscalar& accum_oil = vec_accum_oil[ i ];
      // 1. use old accumulation value;
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
      ADscalar& accum_wat = vec_accum_wat[ i ];
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

      std::cout << "----------------------------------- "<< left_gdbk << " --- " << right_gdbk << std::endl;
      // std::cout << "--- potential_water: "<< potential_wat.value() << std::endl;
      // std::cout << "--- transmiss_water: "<< trans_wat.value() << std::endl;
      // std::cout << "--- flux_water: "     << flux_wat.value() << std::endl;
      // std::cout << "--- potential_oil: "<< potential_oil.value() << std::endl;
      // std::cout << "--- transmiss_oil: "<< trans_oil.value() << std::endl;
      // std::cout << "--- flux_oil: "     << flux_oil.value() << std::endl;
      
      //std::cout << prop.viscosity_wat << std::endl;
      //std::cout << prop.krw << std::endl;      
      //std::cout << trans_oil << std::endl;

      
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
