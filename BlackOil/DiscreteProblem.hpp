#pragma once
//#include <functional>

#include "FluidRockProperty.hpp"
#include "Gridding.hpp"

struct CentralProperty
{
   ADscalar kro, krw;
   ADscalar pw, pc;
   ADscalar bo, bw;
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


   void     Update_Property( std::size_t _grid_number,
			     ADscalar (*_Pc_Function)( const ADscalar& ) );
   //std::function< ADscalar(const ADscalar&) > _Pc_Function );
      void     Update_Property( ADscalar (*_Pc_Function)( const ADscalar& ) );
   //std::function< ADscalar(const ADscalar&) > _Pc_Function );
   double   Accum_Constant ( double _time_step, double _del_x, double _del_y, double _del_z );
   ADscalar Accum_Term     ( double _accum_constant, const ADscalar& _porosity,
			     const ADscalar& _saturation, const ADscalar& _volume_factor );
   ADscalar Flux_Term      ( const ADscalar& _transmissibility, const ADscalar& _potential);

   void     Fill_HorizontalLayer ( std::size_t _Kth_layer, double _Po, double _Sw );
   void     Initialize     ();

   void     Evaluate_Accumulate ( double _time_step,
				  ADscalar (*_Pc_Function)( const ADscalar& ) );
   //std::function< ADscalar(const ADscalar&) > _Pc_Function );
   CentralProperty& GetUpstream ( const ADscalar& _pressure_left,  const ADscalar& _pressure_right,
				  std::size_t _gdbk_left, std::size_t _gdbk_right );
   void     Evaluate_Flux  ();
   void     Evaluate       ( double _time_step,
			     ADscalar (*_Pc_Function)( const ADscalar& ) );
   //std::function< ADscalar(const ADscalar&) > _Pc_Function );
   

private:
   Cartesian3D                   grid;
   ConnectionList                CList;
   
   ADvector                      vec_accum;
   std::vector<double>           vec_accum_old;
   vector< CentralProperty >     vec_centralprop;
   InterfcProperty               interfcprop;

   ADvector vec_unknown; // water saturation & oil pressure
   ADvector residual;

private:
   std::size_t GetPoIndex( std::size_t _GridBlockIndex ) { return 2 * _GridBlockIndex; }
   std::size_t GetPoIndex( std::size_t _i, std::size_t _j, std::size_t _k )
   {
      return 2 * grid.GetGridBlockIndex( _i, _j, _k );
   }
   std::size_t GetSwIndex( std::size_t _GridBlockIndex ) { return GetPoIndex( _GridBlockIndex ) + 1; }
   std::size_t GetSwIndex( std::size_t _i, std::size_t _j, std::size_t _k )
   {
      return GetPoIndex( _i, _j, _k ) + 1;
   }
   

};


// ================================================================================

DiscreteProblem :: DiscreteProblem ()
{
   grid.Set( 200, 200, 25, 20, 20, 3 );
   CList.Initialize( grid );
   interfcprop.trans_constant_x = PROPERTY::Trans_Constant( grid.Ny, grid.Nz, grid.Nx );
   interfcprop.trans_constant_y = PROPERTY::Trans_Constant( grid.Nx, grid.Nz, grid.Ny );
   interfcprop.trans_constant_z = PROPERTY::Trans_Constant( grid.Nx, grid.Ny, grid.Nz );
   
   PROPERTY :: Generate_Capillary_Table( "3_layer_reservoir_properties/capillary.txt" );
   PROPERTY :: Read_Layer_Property     ( "3_layer_reservoir_properties/layer_1.txt" );
   PROPERTY :: Read_Layer_Property     ( "3_layer_reservoir_properties/layer_2.txt" );
   PROPERTY :: Read_Layer_Property     ( "3_layer_reservoir_properties/layer_3.txt" );

   const std::size_t Grid_Number = grid.Nx * grid.Ny * grid.Nz;
   const std::size_t Var_Number = 2 * Grid_Number;
   
   vec_centralprop.resize( Grid_Number );
   vec_unknown.resize    ( Var_Number );
   residual.resize       ( Var_Number );
}

void DiscreteProblem :: Update_Property( std::size_t _grid_number,
					 ADscalar (*_Pc_Function)( const ADscalar& ) )
//std::function< ADscalar(const ADscalar&) > _Pc_Function )
{
   std::cout << "---------------------------------------- "<< _grid_number << std::endl;
   // evaluate properties at the center of grid blocks
   CentralProperty& prop = vec_centralprop[ _grid_number ];
   const std::size_t PoIndex = GetPoIndex( _grid_number );
   const std::size_t SwIndex = GetSwIndex( _grid_number );
   const ADscalar& _Po = vec_unknown[ PoIndex ];
   const ADscalar& _Sw = vec_unknown[ SwIndex ];

   ADscalar S_wd = PROPERTY::Swd( _Sw );
   prop.kro = PROPERTY::Kro( S_wd );
   prop.krw = PROPERTY::Krw( S_wd );
	    
   prop.pc  = _Pc_Function( _Sw );
   // std::cout << "wtf_Sw =\t" << _Sw.value() << std::endl;
   // std::cout << "wtf_Po =\t" << _Po.value() << std::endl;
   // std::cout << "wtf_Pc =\t" << prop.pc.value() << std::endl;
   ADscalar p_w = PROPERTY::Pw( _Po, prop.pc );
   //std::cout << "wft_pw =\t" << p_w.value() << std::endl;
   prop.pw = p_w;
   ADscalar b_o = PROPERTY::Bo ( _Po );
   ADscalar b_w = PROPERTY::Bw ( prop.pw );
   prop.bo = b_o;
   prop.bw = b_w;
   ADscalar dens_o = PROPERTY::Density_Oil( b_o );
   ADscalar dens_w = PROPERTY::Density_Wat( b_w );
   ADscalar spwt_o = PROPERTY::SpecificWeight( dens_o );
   ADscalar spwt_w = PROPERTY::SpecificWeight( dens_w );	    
   prop.spec_weight_oil = spwt_o;
   prop.spec_weight_wat = spwt_w;
   prop.viscosity_oil = PROPERTY::VisO;
   prop.viscosity_wat = PROPERTY::VisW;
   double porosity_b = PROPERTY::vec_porosity_b[ _grid_number ];
   prop.porosity    = PROPERTY::Porosity( _Po, porosity_b );


   std::cout << "****************"<< std::endl;
   std::cout << "Po =\t" << _Po.value() << std::endl;
   std::cout << "Sw =\t" << _Sw.value() << std::endl;
   std::cout << "Pw =\t" << prop.pw.value() << std::endl;   
   std::cout << "Pc =\t" << prop.pc.value() << std::endl;
   std::cout << "****************"<< std::endl;
   std::cout << "Kro =\t" << prop.kro.value() << std::endl;
   std::cout << "Krw =\t" << prop.krw.value() << std::endl;
   std::cout << "Bo =\t" << prop.bo.value() << std::endl;
   std::cout << "Bw =\t" << prop.bw.value() << std::endl;
   std::cout << "spwto =\t" << prop.spec_weight_oil.value() << std::endl;
   std::cout << "spwtw =\t" << prop.spec_weight_wat.value() << std::endl;
   std::cout << "viso =\t" << prop.viscosity_oil.value() << std::endl;
   std::cout << "visw =\t" << prop.viscosity_wat.value() << std::endl;
   std::cout << "poro =\t" << prop.porosity.value() << std::endl;
}

void DiscreteProblem :: Update_Property ( ADscalar (*_Pc_Function)( const ADscalar& ) )
//std::function< ADscalar(const ADscalar&) > _Pc_Function )
{
   std::size_t count_grid_number = 0;
   for( std::size_t k = 0; k < grid.Nz; ++k )
   {
      for( std::size_t j = 0; j < grid.Ny; ++j )
      {
	 for( std::size_t i = 0; i < grid.Nx; ++i )
	 {
	    Update_Property( count_grid_number, _Pc_Function );
	    ++count_grid_number;
	 }
      }
   }
}
	 
double DiscreteProblem :: Accum_Constant ( double _time_step,
					   double _del_x, double _del_y, double _del_z )
{
   return _del_x * _del_y * _del_z / 5.615 / _time_step;
}

ADscalar DiscreteProblem :: Accum_Term ( double _accum_constant, const ADscalar& _porosity,
					 const ADscalar& _saturation, const ADscalar& _volume_factor )
{
   return _accum_constant * ( _porosity * _saturation / _volume_factor );
}

ADscalar DiscreteProblem :: Flux_Term ( const ADscalar& _transmissibility, const ADscalar& _potential)
{
   return _transmissibility * _potential;
}


void DiscreteProblem :: Fill_HorizontalLayer ( std::size_t _Kth_layer, double _Po, double _Sw )
{
   for( std::size_t j = 0; j < grid.Ny; ++j )
   {
      for( std::size_t i = 0; i < grid.Nx; ++i )
      {
	 double var_index = GetPoIndex( i, j, _Kth_layer );
	 vec_unknown[var_index].value() = _Po;
	 vec_unknown[var_index+1].value() = _Sw;
      }
   }
}

void DiscreteProblem :: Initialize ()
{
   // 1. Activation
   grid.Activate_NaturalOrdering( vec_unknown );
   // 2. Find Po & Sw from Pc (drainage) table
   /* Note
    * This part is hard-coded:
    * WOC is located at the bottom boundary, where
    * oil pressure is given as 3500 [psi]
    */
   double Po_woc = PROPERTY::Poi.value();
   double Pc_woc = PROPERTY::Pc_Drainage( 1.0 - PROPERTY::Sor ).value();
   double Pw_woc = Po_woc - Pc_woc;
   const double Pc_Max = PROPERTY::Pc_Drainage( PROPERTY::Siw ).value();
   int track_Kth_layer;
   // 2.1 Bottom Layer
   track_Kth_layer = grid.Nz - 1;
   double Po_bottom = PROPERTY :: Find_Po( Po_woc,
   					   PROPERTY::SpecificWeight( Po_woc ).value(),
   					   grid.Dz/2.0 );
   double Pw_bottom = PROPERTY :: Find_Pw( Pw_woc,
   					   PROPERTY::SpecificWeight( Pw_woc ).value(),
   					   grid.Dz/2.0 );
   double Pc_bottom = Po_bottom - Pw_bottom;

   if( Pc_bottom > Pc_Max )
   {
      Fill_HorizontalLayer( track_Kth_layer, Po_bottom, PROPERTY::Siw.value() );
   }
   else
   {
      double Sw_bottom = PROPERTY :: Sw_FromPc_Drainage( Po_bottom - Pw_bottom );
      Fill_HorizontalLayer( track_Kth_layer, Po_bottom, Sw_bottom );
   }
   
   // 2.2 (Bottom Layer, Top of Transition Zone)

   double Pc_transition;
   double Po_transition = Po_woc;
   double Pw_transition = Pw_woc;
   double Sw_transition;
   --track_Kth_layer;

   while( track_Kth_layer >= 0 )
   {
      Po_transition = PROPERTY :: Find_Po( Po_transition,
   					   PROPERTY::SpecificWeight( Po_transition ).value(),
   					   grid.Dz );
      Pw_transition = PROPERTY :: Find_Pw( Pw_transition,
   					   PROPERTY::SpecificWeight( Pw_transition ).value(),
   					   grid.Dz );
      Pc_transition = Po_transition - Pw_transition;

      if( Pc_transition >= Pc_Max )
      {
   	 // reach the top of transition zone
   	 Fill_HorizontalLayer( track_Kth_layer, Po_transition, PROPERTY::Siw.value() );
      }
      else
      {
   	 Sw_transition = PROPERTY :: Sw_FromPc_Drainage( Po_transition - Pw_transition );
   	 Fill_HorizontalLayer( track_Kth_layer, Po_transition, Sw_transition );
      }
      --track_Kth_layer;
   }

   /*
   for( std::size_t i = 0; i < grid.Nx*grid.Ny*grid.Nz; ++i )
   {
      const std::size_t PoIndex = GetPoIndex( i );
      const std::size_t SwIndex = GetSwIndex( i );
      const ADscalar& Po = vec_unknown[ PoIndex ];
      const ADscalar& Sw = vec_unknown[ SwIndex ];

      
      std::cout << "---------------------------------------- "<< i << std::endl;
      std::cout << "Po =\t" << Po.value() << std::endl;
      std::cout << "Sw =\t" << Sw.value() << std::endl;
      ADscalar pc;
      pc = PROPERTY::Pc_Drainage( Sw );
      std::cout << "Pc =\t" << pc.value() << std::endl;   
      std::cout << "Pw =\t" << PROPERTY::Pw( Po, pc ).value() << std::endl;
      
   }
   */
}

void DiscreteProblem :: Evaluate_Accumulate ( double _time_step,
					      ADscalar (*_Pc_Function)( const ADscalar& ) )
//std::function< ADscalar(const ADscalar&) > _Pc_Function )
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
	    Update_Property( count_grid_number, _Pc_Function );

	    // accumulate
	    const CentralProperty& prop = vec_centralprop[ count_grid_number ];
	    const std::size_t PoIndex = GetPoIndex( count_grid_number );
	    const std::size_t SwIndex = GetSwIndex( count_grid_number );	    
	    
	    // oil phase
	    residual[PoIndex] += Accum_Term( accum_constant,
					     prop.porosity,
					     vec_unknown[SwIndex],
					     prop.bo );
	    // water phase
	    residual[SwIndex] += Accum_Term( accum_constant,
					     prop.porosity,
					     vec_unknown[SwIndex],
					     prop.bw );
	    ++count_grid_number;
	 }
      }
   }
}

CentralProperty& DiscreteProblem :: GetUpstream ( const ADscalar& _pressure_left,
						  const ADscalar& _pressure_right,
						  std::size_t _gdbk_left,
						  std::size_t _gdbk_right )
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

      CentralProperty& prop = GetUpstream( vec_unknown[ PoIndex_L ],
					   vec_unknown[ PoIndex_R ],
					   left_gdbk,
					   right_gdbk );
      

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
      ADscalar flux_oil = Flux_Term                   ( trans_oil, potential_oil );
      ADscalar flux_wat = Flux_Term                   ( trans_wat, potential_wat );

      residual[ PoIndex_L ] -= flux_oil;
      residual[ PoIndex_R ] += flux_oil;
      residual[ SwIndex_L ] -= flux_wat;
      residual[ SwIndex_R ] += flux_wat;
   }
}


void DiscreteProblem :: Evaluate ( double _time_step,
				   ADscalar (*_Pc_Function)( const ADscalar& ) )
//std::function< ADscalar(const ADscalar&) > _Pc_Function )
{
   Evaluate_Accumulate( _time_step, _Pc_Function );
   Evaluate_Flux();
}
