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
   DiscreteProblem ();
   DiscreteProblem ( const std::vector<double>& _data );

private: // for initialization & preparation
   void Add_Wells ();
   void Activate_Uold ();
   void ActivateProblem ();
   void Initialize_K_interface (); 
   void Add_Producer ( double _rate, std::size_t _i, std::size_t j, std::size_t k, double _rw = 0.5 );
   void Add_Injector ( double _bhp,  std::size_t _i, std::size_t j, std::size_t k, double _rw = 0.5 );
   

public:  // for update
   void Activate_U ( const std::vector<double>& _data );
   void Update_Uold ();                                                     // 1
   void Update_Unew ( const std::vector<double>& _negative_newton_update ); // 2
   void Update_Property    ();                                              // 3
   void Update_Accum_Old ( double _time_step );                             // 4

   
public:
   // for history matching
   template< typename CSR >
   void Get_Jacobian_Dm ( double _time_step, CSR& _Jacobian_Dm );

   
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
   ADscalar GetQ ( const Producer& _producer, int _tag )
   {
      ADscalar tmp( 0.0 );
      tmp = Sink_Term( _producer, _tag );
      return tmp;
   }

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
   // to assist Get_Jacobian_Dm
   void Get_Dm ( std::size_t _GBIndex_Central, std::size_t _GBIndex_Other,
		 double _trans_constant, int _is_Zdirection, const std::vector<double>& _vec_K,
		 double& _DO, double & _DW );

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

};


// ================================================================================
// ================================================================================
// ================================================================================
// ================================================================================
// ================================================================================
DiscreteProblem :: DiscreteProblem () : grid( 200, 200, 25, 10, 1, 1 ), //20-30-3
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
   interfcprop.trans_constant_x = PROPERTY::Trans_Constant( grid.Dy, grid.Dz, grid.Dx );
   interfcprop.trans_constant_y = PROPERTY::Trans_Constant( grid.Nx, grid.Nz, grid.Ny );
   interfcprop.trans_constant_z = PROPERTY::Trans_Constant( grid.Nx, grid.Ny, grid.Nz );
 
   constprop.Read_From_File( "3_layer_reservoir_properties/layer_1.txt" );
   constprop.Read_From_File( "3_layer_reservoir_properties/layer_2.txt" );
   constprop.Read_From_File( "3_layer_reservoir_properties/layer_3.txt" );

   interfcprop.K_interface.resize( CList.Size );
   Initialize_K_interface();

   Add_Wells();

   Activate_Uold();
   mUnew = mUold;
   
   Update_Property();
}

DiscreteProblem :: DiscreteProblem ( const std::vector<double>& _data )
   : grid( 200, 200, 25, 10, 1, 1 ), //20-30-3
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
   interfcprop.trans_constant_x = PROPERTY::Trans_Constant( grid.Dy, grid.Dz, grid.Dx );
   interfcprop.trans_constant_y = PROPERTY::Trans_Constant( grid.Nx, grid.Nz, grid.Ny );
   interfcprop.trans_constant_z = PROPERTY::Trans_Constant( grid.Nx, grid.Ny, grid.Nz );
 
   constprop.Read_From_File( "3_layer_reservoir_properties/layer_1.txt" );
   constprop.Read_From_File( "3_layer_reservoir_properties/layer_2.txt" );
   constprop.Read_From_File( "3_layer_reservoir_properties/layer_3.txt" );

   interfcprop.K_interface.resize( CList.Size );
   Initialize_K_interface();

   Add_Wells();
      
   Activate_U( _data );
   mUnew = mUold;
   Update_Property();
}     

// --- ---
void DiscreteProblem :: Add_Wells()
{
   // inject water with constant rate
   // Add_Injector( -100, 6,  6,  2 );
   // Add_Injector( -100, 17, 26, 2 );
   


   // producer with constant bottom pressure
   // Add_Producer( 800, 5,  4,  0 );
   // Add_Producer( 800, 16, 5,  0 );
   // Add_Producer( 800, 6,  23, 0 );
   // Add_Producer( 800, 5,  21, 0 );

   Add_Injector( -200, 0, 0, 0 );
   Add_Producer( 800,  9, 0, 0 );

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
void DiscreteProblem :: Add_Producer  ( double _bhp,
					std::size_t _i, std::size_t _j, std::size_t _k, double _rw )
{
   std::size_t BlockIndex = grid.GetGridBlockIndex( _i, _j, _k );
   Producer well( _bhp, _i, _j, _k, _rw,
		  grid.Dx, grid.Dy,
		  constprop.vec_Kx[ BlockIndex ], constprop.vec_Ky[ BlockIndex ] );

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
   Sw.value() = ( 1.0 - constprop.Sor ).value();
   ADscalar Swd( 0.0 );
   Swd = PROPERTY :: Swd( Sw, constprop.Siw, constprop.Sor );
   ADscalar krw_inj( 0.0 );
   krw_inj = PROPERTY :: Krw( Swd );

   ADscalar ret( 0.0 );
   ret = _injector.WI() * krw_inj / ( prop.viscosity_wat * prop.bw )
     * ( mUnew[ PoIndex ] - mUnew[ BHPIndex ] ) - _injector.Rate();

   // std::cout << "Evaluating BHP Residual -------------------" << std::endl;
   // std::cout << "Swd: " << Swd << std::endl;
   // std::cout << "krw: " << krw_inj << std::endl;
   // std::cout << "vis: " << prop.viscosity_wat << std::endl;
   // std::cout << "bw:  " << prop.viscosity_wat << std::endl;
   // std::cout << "Po:  " << mUnew[ PoIndex ] << std::endl;
   // std::cout << "BHP: " << mUnew[ BHPIndex ] << std::endl;
   // std::cout << "Rate: " << _injector.Rate() << std::endl;
   // std::cout << "R:   " << ret << std::endl;
     
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
   for( int i = 0; i < VarNumber_Added; ++i )
   {
      // 1. Add constant rate to residual equation
      const Injector& well = vec_injector[ i ];
      const std::size_t SwIndex = GetSwIndex( well.I(), well.J(), well.K() );      
      residual[ SwIndex ] += well.Rate();

      // 2. Form new residual equation
      const std::size_t BHPIndex = GetBHPIndex( i );
      residual[ BHPIndex ] = BHP_Residual( vec_injector[i] );
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



void DiscreteProblem :: Get_Dm ( std::size_t _GBIndex_Central, std::size_t _GBIndex_Other, // location
				 double _trans_constant, int _is_Zdirection, // direction related
				 const std::vector<double>& _vec_K,
				 double& _DO, double & _DW ) // result
{
   std::size_t PoIndex_C = GetPoIndex( _GBIndex_Central );
   std::size_t PoIndex_O = GetPoIndex( _GBIndex_Other );

   const CentralProperty& prop = GetUpstream( _GBIndex_Central,
					      mUnew[ PoIndex_C ],
					      _GBIndex_Other,
					      mUnew[ PoIndex_O ] );
	       
   double PO = PROPERTY ::Potential( mUnew[ PoIndex_C ].value(),
				     mUnew[ PoIndex_O ].value(),
				     prop.spec_weight_oil.value(),
				     grid.Dz * _is_Zdirection );

   double PW = PROPERTY :: Potential( mUnew[ PoIndex_C ].value(),
				      mUnew[ PoIndex_O ].value(),
				      prop.spec_weight_wat.value(),
				      grid.Dz * _is_Zdirection );

   double MO = prop.kro.value() / ( prop.viscosity_oil.value() * prop.bo.value() );
   double MW = prop.krw.value() / ( prop.viscosity_wat.value() * prop.bw.value() );

   double K_C = _vec_K[ _GBIndex_Central ];
   double K_O = _vec_K[ _GBIndex_Other   ];

   double K_Sum = K_C + K_O;
   double K_Term = 2 * K_C * K_C * K_O / ( K_Sum * K_Sum );
	       
   _DO = PO * _trans_constant * MO * K_Term;
   _DW = PW * _trans_constant * MW * K_Term;
}

template< typename CSR >
void DiscreteProblem ::Get_Jacobian_Dm ( double _time_step, CSR& _Jacobian_Dm )
{
   int track_nzv(0);
   _Jacobian_Dm.Row().clear();
   _Jacobian_Dm.Col().clear();
   _Jacobian_Dm.NZV().clear();
   
   std::size_t GBIndex_C(0);
   std::size_t GBIndex_IL(0);
   std::size_t GBIndex_IR(0); 
   std::size_t GBIndex_JL(0); 
   std::size_t GBIndex_JR(0); 
   std::size_t GBIndex_KL(0); 
   std::size_t GBIndex_KR(0); 
   
   for( std::size_t k = 0; k < grid.Nz; ++k )
   {
      for( std::size_t j = 0; j < grid.Ny; ++j )
      {
	 for( std::size_t i = 0; i < grid.Nx; ++i )
	 {
	    // 1. Ln( Kx ) ++++++++++++++++++++++++++++++++++++++++++++
	    
	    // Indicator
	    bool IL = false, IR = false;
	    bool JL = false, JR = false;
	    bool KL = false, KR = false;

	    // Absoluate Permeability
	    double K_Cx, K_Cy, K_Cz;
	    double K_IL(0.0), K_IR(0.0);
	    double K_JL(0.0), K_JR(0.0);
	    double K_KL(0.0), K_KR(0.0);

	    // Derivatives   
	    double DO_Cx, DO_Cy, DO_Cz;
	    double DO_IL(0.0), DO_IR(0.0); // partial derivative ( oil equation )
	    double DO_JL(0.0), DO_JR(0.0);
	    double DO_KL(0.0), DO_KR(0.0);
	    double DW_Cx, DW_Cy, DW_Cz;
	    double DW_IL(0.0), DW_IR(0.0); // partial derivative ( wat equation )
	    double DW_JL(0.0), DW_JR(0.0);
	    double DW_KL(0.0), DW_KR(0.0);
	    
	    GBIndex_C = GetGBIndex( i, j, k );
	    K_Cx = constprop.vec_Kx[ GBIndex_C ];
	    K_Cy = constprop.vec_Ky[ GBIndex_C ];
	    K_Cz = constprop.vec_Kz[ GBIndex_C ];

	    // i-1
	    if( i != 0 )
	    {
	       IL = true;
	       GBIndex_IL = GetGBIndex( i-1, j, k );
	       K_IL = constprop.vec_Kx[ GBIndex_IL ];
	       Get_Dm( GBIndex_C, GBIndex_IL,
		       interfcprop.trans_constant_x, 0, constprop.vec_Kx,
		       DO_IL, DW_IL );
	    }
	    
	    // i+1
	    if( i != grid.Nx-1 )
	    {
	       IR = true;
	       GBIndex_IR = GetGBIndex( i+1, j, k );
	       K_IR = constprop.vec_Kx[ GBIndex_IR ];
	       Get_Dm( GBIndex_C, GBIndex_IR,
		       interfcprop.trans_constant_x, 0, constprop.vec_Kx,
		       DO_IR, DW_IR );
	    }
	    // j-1
	    if( j != 0 )
	    {
	       JL = true;
	       GBIndex_JL = GetGBIndex( i, j-1, k );
	       K_JL = constprop.vec_Ky[ GBIndex_JL ];
	       Get_Dm( GBIndex_C, GBIndex_JL,
		       interfcprop.trans_constant_y, 0, constprop.vec_Ky,
		       DO_JL, DW_JL );
	    }

	    // j+1
	    if( j != grid.Ny-1 )
	    {
	       JR = true;
	       GBIndex_JR = GetGBIndex( i, j+1, k );
	       K_JR = constprop.vec_Ky[ GBIndex_JR ];
	       Get_Dm( GBIndex_C, GBIndex_JR,
		       interfcprop.trans_constant_y, 0, constprop.vec_Ky,
		       DO_JR, DW_JR );
	    }
	    
	    // k-1
	    if( k != 0 )
	    {
	       KL = true;
	       GBIndex_KL = GetGBIndex( i, j, k-1 );
	       K_KL = constprop.vec_Kz[ GBIndex_KL ];
	       Get_Dm( GBIndex_C, GBIndex_KL,
		       interfcprop.trans_constant_z, 1, constprop.vec_Kz,
		       DO_KL, DW_KL );
	    }

	    // k+1
	    if( k != grid.Nz-1 )
	    {
	       KR = true;
	       GBIndex_KR = GetGBIndex( i, j, k+1 );
	       K_KR = constprop.vec_Kz[ GBIndex_KR ];
	       Get_Dm( GBIndex_C, GBIndex_KR,
		       interfcprop.trans_constant_z, 1, constprop.vec_Kz,
		       DO_KR, DW_KR );
	    }

	    // i,j,k --> Central
	    DO_Cx = ( DO_IL * K_IL + DO_IR * K_IR ) / K_Cx;
	    DO_Cy = ( DO_JL * K_JL + DO_JR * K_JR ) / K_Cy;
	    DO_Cz = ( DO_KL * K_KL + DO_KR * K_KR ) / K_Cz;

	    DW_Cx = ( DW_IL * K_IL + DW_IR * K_IR ) / K_Cx;
   	    DW_Cy = ( DW_JL * K_JL + DW_JR * K_JR ) / K_Cy;
	    DW_Cz = ( DW_KL * K_KL + DW_KR * K_KR ) / K_Cz;

	    // 2. porosity ++++++++++++++++++++++++++++++++++++++++++++
	    double accum_const = Accum_Constant( _time_step, grid.Dx, grid.Dy, grid.Dz );
	    const CentralProperty& prop = vec_centralprop[ GBIndex_C ];
	    std::size_t SwIndex_C = GetSwIndex( GBIndex_C );
	    double Sw_C = mUnew[ SwIndex_C ].value();
	    double So_C = 1.0 - Sw_C;
	    // derivatives w.r.t porosity
	    double DOp = accum_const * prop.porosity.value() * So_C / prop.bo.value();
	    double DWp = accum_const * prop.porosity.value() * Sw_C / prop.bw.value();

	    
	    // 3. place derivatives +++++++++++++++++++++++++++++++++++++++++++++++++++
	    std::vector<int> tmp_wat_placement_col;
	    std::vector<double> tmp_wat_placement_nzv;
	    // 3.1 oil equation w.r.t Ln(k)
	    _Jacobian_Dm.Row().push_back( track_nzv );
	    if( true == KL ) // k-1, o
	    {
	       int KL_ColIndex = GetKzColIndex( GBIndex_KL );
	       _Jacobian_Dm.Col().push_back( KL_ColIndex );
	       _Jacobian_Dm.NZV().push_back( DO_KL );
	       track_nzv++;
	       
	       tmp_wat_placement_col.push_back( KL_ColIndex );
	       tmp_wat_placement_nzv.push_back( DW_KL );
	    }
	    if( true == JL ) // j-1, o
	    {
	       int JL_ColIndex = GetKyColIndex( GBIndex_JL );
	       _Jacobian_Dm.Col().push_back( JL_ColIndex );
	       _Jacobian_Dm.NZV().push_back( DO_JL );
	       track_nzv++;

	       tmp_wat_placement_col.push_back( JL_ColIndex );
	       tmp_wat_placement_nzv.push_back( DW_JL );	       
	    }
	    if( true == IL ) // i-1, o
	    {
	       int IL_ColIndex = GetKxColIndex( GBIndex_IL );
	       _Jacobian_Dm.Col().push_back( IL_ColIndex );
	       _Jacobian_Dm.NZV().push_back( DO_IL );
	       track_nzv++;

	       tmp_wat_placement_col.push_back( IL_ColIndex );
	       tmp_wat_placement_nzv.push_back( DW_IL );	       	       
	    }
	    if( true ) // central, o
	    {
	       // central-kx
	       int Cx_ColIndex = GetKxColIndex( GBIndex_C );
	       _Jacobian_Dm.Col().push_back( Cx_ColIndex );
	       _Jacobian_Dm.NZV().push_back( DO_Cx );
	       track_nzv++;

	       tmp_wat_placement_col.push_back( Cx_ColIndex );
	       tmp_wat_placement_nzv.push_back( DW_Cx );
	       
	       // central-ky
	       int Cy_ColIndex = GetKyColIndex( GBIndex_C );
	       _Jacobian_Dm.Col().push_back( Cy_ColIndex );
	       _Jacobian_Dm.NZV().push_back( DO_Cy );
	       track_nzv++;

	       tmp_wat_placement_col.push_back( Cy_ColIndex );
	       tmp_wat_placement_nzv.push_back( DW_Cy );
	       
	       // central-kz
	       int Cz_ColIndex = GetKzColIndex( GBIndex_C );
	       _Jacobian_Dm.Col().push_back( Cz_ColIndex );
	       _Jacobian_Dm.NZV().push_back( DO_Cz );
	       track_nzv++;

	       tmp_wat_placement_col.push_back( Cz_ColIndex );
	       tmp_wat_placement_nzv.push_back( DW_Cz );
	    }
	    
	    if( true == IR ) // i+1, o
	    {
	       int IR_ColIndex = GetKxColIndex( GBIndex_IR );
	       _Jacobian_Dm.Col().push_back( IR_ColIndex );
	       _Jacobian_Dm.NZV().push_back( DO_IR );
	       track_nzv++;

	       tmp_wat_placement_col.push_back( IR_ColIndex );
	       tmp_wat_placement_nzv.push_back( DW_IR );
	    }
	    if( true == JR ) // j+1, o
	    {
	       int JR_ColIndex = GetKyColIndex( GBIndex_JR );
	       _Jacobian_Dm.Col().push_back( JR_ColIndex );
	       _Jacobian_Dm.NZV().push_back( DO_JR );
	       track_nzv++;

	       tmp_wat_placement_col.push_back( JR_ColIndex );
	       tmp_wat_placement_nzv.push_back( DW_JR );
	    }

	    if( true == KR ) // k+1, o
	    {
	       int KR_ColIndex = GetKzColIndex( GBIndex_KR );
	      
	       _Jacobian_Dm.Col().push_back( KR_ColIndex );
	       _Jacobian_Dm.NZV().push_back( DO_KR );
	       track_nzv++;

	       tmp_wat_placement_col.push_back( KR_ColIndex );
	       tmp_wat_placement_nzv.push_back( DW_KR );
	    }


	    // 3.2 oil equation w.r.t porosity
	    int PSO_ColIndex = GetPSColIndex( GBIndex_C );
	    _Jacobian_Dm.Col().push_back( PSO_ColIndex );
	    _Jacobian_Dm.NZV().push_back( DOp );
	    track_nzv++;

	    
	    
	    // 3.3 water equation w.r.t Ln(k)
	    _Jacobian_Dm.Row().push_back( track_nzv );

	    for( std::size_t i = 0; i < tmp_wat_placement_col.size(); ++i )
	    {
	       _Jacobian_Dm.Col().push_back( tmp_wat_placement_col[ i ] );
	       _Jacobian_Dm.NZV().push_back( tmp_wat_placement_nzv[ i ] );
	       track_nzv++;
	    }

	    // 3.4 water equation w.r.t porosity
	    int PSW_ColIndex = GetPSColIndex( GBIndex_C );
	    _Jacobian_Dm.Col().push_back( PSW_ColIndex );
	    _Jacobian_Dm.NZV().push_back( DWp );
	    track_nzv++;
	 }
      }
   }
   _Jacobian_Dm.Row().push_back( track_nzv );
   _Jacobian_Dm.nRow() = VarNumber_Total;
   _Jacobian_Dm.nCol() = 4 * GridNumber;
   _Jacobian_Dm.nNZV() = track_nzv;

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
