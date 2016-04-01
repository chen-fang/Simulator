#pragma once
#include <vector>

#include "PropertyCalculator.hpp"
#include "Drift-Flux/EPRI.hpp"

typedef std::vector<double> Vec;

class DiscreteProblem
{
//private:
public:
   enum PhaseID { L=0, G=1 };

   typedef struct
   {
      ADs HG;
      ADs P;
   }                                         CellStateElement;
   typedef std::vector< CellStateElement >   CellStateVector;

   typedef struct
   {
      ADs VL;
      ADs VG;
   }                                         FaceStateElement;
   typedef std::vector< FaceStateElement >   FaceStateVector;

   typedef struct
   {
      ADs HG;
      ADs P;
      ADs VL;
      ADs VG;
   }                                         StateElement;
   typedef std::vector< StateElement >       StateVector;
   
   // < staggerd grid > variables defined in CV center
   typedef struct
   {
      ADs Den[2]; // CellProp_PropCal
      ADs DenM;   // CellProp_PropCal
      ADs VelM;   // CellProp_Upwind
   }                                         CellProps;
   typedef std::vector< CellProps >          CellVector;

   // < staggerd grid > variables defined in CV face
   typedef struct
   {
      ADs H_[2];
      ADs Den_[2];
      ADs Vis_[2];
      ADs H_Den_[2]; // FaceProp_Upwind
      ADs H_Vis_[2]; // FaceProp_Upwind
      ADs DenM_;     // FaceProp_Upwind
      ADs VisM_;
      ADs VelM_;     // FaceProp_Upwind
   }                                        FaceProps;
   typedef std::vector< FaceProps >         FaceVector;

   typedef struct
   {
      double MatBal[2];
      double NormSat[2];
   }                                        ConvergenceInfo;

   enum SSMode{ FixQ=0, FixP=1 };
   typedef struct
   {
      std::size_t loc;
      bool        is_source;
      int         mode;       // [0]-const rate / [1]-pressure
      double      Q[2];       // for: const rate
      double      P;          // for: const pressure
   }                                        SS;
      
//private:
public:
   // Data from previous time step
   //
   // MeshType mesh;
   // PropertyCalculator PropCalc;

   const int PHASE_NUM;
   const int TOTAL_CELL_NUM;
   const int TOTAL_FACE_NUM;
   const int dX;
   const int D;
   const double g;
   const double surface_tension;
   const double theta;
   
   CellVector Cprop;
   FaceVector Fprop;

   CellStateVector Cstate;
   FaceStateVector Fstate;

   // for mass
   ADv mass_accumL;
   ADv mass_accumG;
   ADv mass_transL;
   ADv mass_transG;
   Vec mass_accumL_old; 
   Vec mass_accumG_old;
   // for momentum
   ADv momt_accum;
   Vec momt_accum_old;
   ADv momt_trans;
   ADv momt_press;
   ADv momt_frict;
   ADv momt_gravt;
   // for source_sink
   std::vector< SS > ss;
   
   ADv residual;

   FluidProperty PropCalc;
   EPRI epri;

public:
   DiscreteProblem( int _PHASE_NUM, int _TOTAL_CELL_NUM, int _dX, int _D, double _theta );

   // Node and face indices
   // '|' indicates CV face
   // '.' indicates CV nodal points
   // 'X' indicates CV faces that are NOT indexed
   //
   // Node Index: | 0 | 1 | 2 | ... | N-1 |   --- --- Total: N
   // Face Index: X . 0 . 1 . 2 ... N-2 . X   --- --- Total: N-1
   //
   // Left and Right Cell Index
   int Get_LC_Index( int face_index )   { return face_index; }
   int Get_RC_Index( int face_index )   { return face_index + 1; }
   // Left and Right Face Index
   int Get_LF_Index( int cell_index )   { return cell_index - 1; }
   int Get_RF_Index( int cell_index )   { return cell_index;   }

   std::size_t max_num_eqns() const     { return PHASE_NUM * ( TOTAL_CELL_NUM + TOTAL_CELL_NUM ); }
   std::size_t max_num_nnz()  const     { return 44 * max_num_eqns() - 28; }

   int MassEqnID( int c, int phaseID )  { return 4 * c + phaseID; }
   int MomtEqnID( int f )               { return 4 * f + 2; }
   int DFlxEqnID( int f )               { return 4 * f + 3; }

   void Transfer_StateVector( const StateVector& state );
   void initiate_state( StateVector& state );
   void bind_to_old_state( const StateVector& old_state );
   bool discretize( const StateVector& state, double dT );

   template< typename R >
   void update_state( StateVector& state, const R& update, bool do_safeguard );
   
   template< typename V, typename M >
   void extract_R_J( V&r, M& m, std::size_t offset );


   // Section: Property
   void Compute_CellProp_PropCacl();
   void Compute_FaceProp_Upwind();
   void Compute_CellProp_Upwind();
   void Compute_Properties(); // summary

   // Section: Mass
   void Compute_Mass_Accum();
   void Compute_Mass_Trans();
   void Compute_Mass_Resid( double dT );

   // Section: Momentum
   void Compute_Momt_Accum();
   void Compute_Momt_Trans();
   void Compute_Momt_Frict();
   void Compute_Momt_Press();
   void Compute_Momt_Gravt();
   void Compute_Momt_Resid( double dT );

   // Section: Drift-Flux
   void Compute_DFlx_Resid();

   // Source & Sink
   void Setup_SS();
   void Compute_SS();
};

DiscreteProblem::DiscreteProblem( int _PHASE_NUM, int _TOTAL_CELL_NUM, int _DX, int _D, double _theta ) :
   PHASE_NUM( _PHASE_NUM ),
   TOTAL_CELL_NUM( _TOTAL_CELL_NUM ),
   TOTAL_FACE_NUM( _TOTAL_CELL_NUM - 1 ),
   dX( _DX ),
   D( _D ),
   g( 9.80665 ), // [ m/s2 ]
   surface_tension( 72.75E-03 ), // [ N/m ]
   theta( _theta ),
   Cprop(  TOTAL_CELL_NUM ),
   Fprop(  TOTAL_FACE_NUM ),
   Cstate( TOTAL_CELL_NUM ),
   Fstate( TOTAL_FACE_NUM )
{
   Setup_SS();
}

void DiscreteProblem::
Transfer_StateVector( const StateVector& state )
{
   for( int c = 0; c < TOTAL_CELL_NUM; ++c )
   {
      Cstate[c].HG = state[c].HG;
      Cstate[c].P  = state[c].P;
   }
   for( int f = 0; f < TOTAL_FACE_NUM; ++f )
   {
      Fstate[f].VL = state[f].VL;
      Fstate[f].VG = state[f].VG;
   }
}

void DiscreteProblem::
initiate_state( StateVector& state )
{
   for( int c = 0; c < TOTAL_CELL_NUM; ++c )
   {
      state[c].HG = 0.4;
      state[c].P  = 101.325E+03; // Pa
      //
      state[c].HG.make_independent( MassEqnID(c,PhaseID::L) );
      //state[c].P.make_independent(  MassEqnID(c,PhaseID::G) );
   }
   for( int f = 0; f < TOTAL_FACE_NUM; ++f )
   {
      state[f].VL = 0.0;
      state[f].VG = 0.0;
      //
      state[f].VL.make_independent( MomtEqnID(f) );
      state[f].VG.make_independent( DFlxEqnID(f) );
   }
}

void DiscreteProblem::
bind_to_old_state( const StateVector& old_state )
{
   Transfer_StateVector( old_state );
   Compute_Properties();
   Compute_Mass_Accum();
   Compute_Momt_Accum();
   for( int c = 0; c < TOTAL_CELL_NUM; ++c )
   {
      mass_accumL_old[c] = mass_accumL[c].value();
      mass_accumG_old[c] = mass_accumG[c].value();
   }
   for( int f = 0; f < TOTAL_FACE_NUM; ++f )
   {
      momt_accum_old[f] = momt_accum[f].value();
   }
}

bool DiscreteProblem::
discretize( const StateVector& state, double dT )
{
   bool is_badvalue = false;
   Transfer_StateVector( state );
   Compute_Properties();
   Compute_Mass_Resid( dT );
   Compute_SS();
   Compute_Momt_Resid( dT );
   Compute_DFlx_Resid();

   std::size_t eqn_massL, eqn_massG;
   for( int c = 0; c < TOTAL_CELL_NUM; ++c )
   {
      eqn_massL = MassEqnID( c, PhaseID::L );
      eqn_massG = MassEqnID( c, PhaseID::G );
      if( !std::isfinite( residual[ eqn_massL ].value() ) ) is_badvalue = true;
      if( !std::isfinite( residual[ eqn_massG ].value() ) ) is_badvalue = true;
   }
   std::size_t eqn_momt, eqn_dflx;
   for( int f = 0; f < TOTAL_FACE_NUM; ++f )
   {
      eqn_momt = MomtEqnID( f );
      eqn_dflx = DFlxEqnID( f );
      if( !std::isfinite( residual[ eqn_momt ].value() ) ) is_badvalue = true;
      if( !std::isfinite( residual[ eqn_dflx ].value() ) ) is_badvalue = true;
   }
   return is_badvalue;
}

template< typename R >
void DiscreteProblem::
update_state( StateVector& state, const R& update, bool do_safeguard )
{
   std::size_t eqn_massL, eqn_massG;
   for( int c = 0; c < TOTAL_CELL_NUM; ++c )
   {
      eqn_massL = MassEqnID( c, PhaseID::L );
      eqn_massG = MassEqnID( c, PhaseID::G );

      state[c].HG += update[ eqn_massL ];
      state[c].P  += update[ eqn_massG ];
   }
   std::size_t eqn_momt, eqn_dflx;
   for( int f = 0; f < TOTAL_FACE_NUM; ++f )
   {
      eqn_momt = MomtEqnID( f );
      eqn_dflx = DFlxEqnID( f );

      state[f].VL += update[ eqn_momt ];
      state[f].VG += update[ eqn_dflx ];
   }   
}

template< typename V, typename M >
void DiscreteProblem::
extract_R_J( V&r, M& m, std::size_t offset )
{
   residual.extract_CSR( r, m.rowptr(), m.colind(), m.value() );
   for( std::size_t i = 0; i < m.rowptr().size(); ++i ) m.rowptr()[i] += offset;
   for( std::size_t i = 0; i < m.colind().size(); ++i ) m.colind()[i] += offset;
   m.check_size();
   if( r.size() != residual.size() )
      std::cout << "BUG IN JACOBIAN\t ZERO ROW FOUND" << r.size() <<"!="<< residual.size() << std::endl;
   m.ainfo.elliptic_varr_id = 0;
   m.ainfo.n_vars_per_block = 4;
   m.ainfo.tile_offset      = 0;
   m.ainfo.n_blocks         = TOTAL_CELL_NUM; // I doubt it
   m.ainfo.n_unblocked_vars = 0;
}


// -----------------------------------------------------
// ------------- Staggered grid scheme -----------------
// -----------------------------------------------------
// Section: Property
void DiscreteProblem::
Compute_CellProp_PropCacl()
{
   // DenL, DenG, DenM
   for( int c = 0; c < TOTAL_CELL_NUM; ++c )
   {
      Cprop[c].Den[PhaseID::L] = PropCalc.Density_Wat( Cstate[c].P );
      Cprop[c].Den[PhaseID::G] = PropCalc.Density_Air( Cstate[c].P, 50.0 );
      Cprop[c].DenM = PropCalc.Density_Mix( Cstate[c].HG * Cprop[c].Den[PhaseID::L],
					    Cstate[c].HG * Cprop[c].Den[PhaseID::G] );
   }
}

// Section: Property
void DiscreteProblem::
Compute_FaceProp_Upwind()
{
   // H_Den_, DenM_, VelM_
   ADs HL( 0.0 ), HG( 0.0 );

   for( int f = 0; f < TOTAL_FACE_NUM; ++f )
   {
      int LC = Get_LC_Index( f );
      int RC = Get_RC_Index( f );
      int C;
      // liquid
      if( Fstate[f].VL.value() >= 0.0 )   C = LC;
      else                                C = RC;

      HL = 1.0 - Cstate[C].HG;
      Fprop[f].H_[PhaseID::L]     = HL;
      Fprop[f].Den_[PhaseID::L]   = Cprop[C].Den[PhaseID::L];
      Fprop[f].Vis_[PhaseID::L]   = PropCalc.Viscosity_Wat();
      Fprop[f].H_Den_[PhaseID::L] = HL * Cprop[C].Den[PhaseID::L];

      // gas
      if( Fstate[f].VG.value() >= 0.0 )   C = LC;
      else                                C = RC;
      HG = Cstate[C].HG;
      Fprop[f].H_[PhaseID::G]     = HG;
      Fprop[f].Den_[PhaseID::G]   = Cprop[C].Den[PhaseID::G];
      Fprop[f].Vis_[PhaseID::L]   = PropCalc.Viscosity_Wat();
      Fprop[f].H_Den_[PhaseID::G] = HG * Cprop[C].Den[PhaseID::G];

      // mixture: DenM_, VelM_
      Fprop[f].DenM_ = PropCalc.Density_Mix(   Fprop[f].H_Den_[PhaseID::L],
					       Fprop[f].H_Den_[PhaseID::G] );
      Fprop[f].VisM_ = PropCalc.Viscosity_Mix( Fprop[f].H_[PhaseID::L], Fprop[f].Vis_[PhaseID::L],
					       Fprop[f].H_[PhaseID::G], Fprop[f].Vis_[PhaseID::G] );

      Fprop[f].VelM_ = ( Fprop[f].H_Den_[PhaseID::L] * Fstate[f].VL +
			 Fprop[f].H_Den_[PhaseID::G] * Fstate[f].VG ) / Fprop[f].DenM_;
   }
}

// Section: Property
void DiscreteProblem::
Compute_CellProp_Upwind()
{
   // VelM
   for( int c = 0; c < TOTAL_CELL_NUM; ++c )
   {
      // 1st cell
      if( c == 0 )
      {
	 int RF = Get_RF_Index(c);
	 Cprop[c].VelM = Fprop[RF].VelM_; // / 2.0;
      }
      // last cell
      else if( c == TOTAL_CELL_NUM -1 )
      {
	 int LF = Get_LF_Index(c);
	 Cprop[c].VelM = Fprop[LF].VelM_; // / 2.0;
      }
      // other cells
      else
      {
	 int LF = Get_LF_Index(c);
	 int RF = Get_RF_Index(c);

	 // liquid
	 // < cocurrent right >
	 if( Fprop[LF].VelM_.value() > 0.0 && Fprop[RF].VelM_.value() > 0.0 )
	 {
	    Cprop[c].VelM = Fprop[LF].VelM_;
	 }
	 // < cocurrent left >
	 else if( Fprop[LF].VelM_.value() > 0.0 && Fprop[RF].VelM_.value() > 0.0 )
	 {
	    Cprop[c].VelM = Fprop[RF].VelM_;
	 }
	 // < counter-current >
	 else
	 {
	    Cprop[c].VelM = (Fprop[LF].VelM_ + Fprop[RF].VelM_) / 2.0;
	 }
      }
   }
}

// Section: Property
void DiscreteProblem::
Compute_Properties()
{
   Compute_CellProp_PropCacl();
   Compute_FaceProp_Upwind();
   Compute_CellProp_Upwind();
}

// Section: Mass
void DiscreteProblem::
Compute_Mass_Accum()
{
   for( int c = 0; c < TOTAL_CELL_NUM; ++c )
   {
      mass_accumL[c] = (1.0 - Cstate[c].HG) * Cprop[c].Den[PhaseID::L];
      mass_accumG[c] =        Cstate[c].HG  * Cprop[c].Den[PhaseID::G];
   }
}

// Section: Mass
void DiscreteProblem::
Compute_Mass_Trans()
{
   for( int c = 0; c < TOTAL_CELL_NUM; ++c )
   {
      mass_transL[c] = 0.0;
      mass_transG[c] = 0.0;
   }

   ADs tmpL( 0.0 ), tmpG( 0.0 );
   for( int f = 0; f < TOTAL_FACE_NUM; ++f )
   {
      int LC = Get_LC_Index(f);
      int RC = Get_RC_Index(f);

      tmpL = Fprop[f].H_Den_[PhaseID::L] * Fstate[f].VL / dX;
      tmpG = Fprop[f].H_Den_[PhaseID::G] * Fstate[f].VG / dX;

      mass_transL[LC] += tmpL;
      mass_transG[LC] += tmpG;

      mass_transL[RC] -= tmpL;
      mass_transG[RC] -= tmpG;
   }
}

// Section: Mass
void DiscreteProblem::
Compute_Mass_Resid( double dT )
{
   Compute_Mass_Accum();
   Compute_Mass_Trans();

   // !!!!!!!!!!!!!!!!!!!!!!!
   // add source term here...

   int MassL_ID, MassG_ID;
   for( int c = 0; c < TOTAL_CELL_NUM; ++c )
   {
      MassL_ID = MassEqnID( c, PhaseID::L );
      MassG_ID = MassEqnID( c, PhaseID::G );
      
      residual[ MassL_ID ] = ( mass_accumL[c] - mass_accumL_old[c] ) / dT + mass_transL[c];
      residual[ MassG_ID ] = ( mass_accumL[c] - mass_accumG_old[c] ) / dT + mass_transG[c];
   }


}


// Section: Momentum
void DiscreteProblem::
Compute_Momt_Accum()
{
   for( int f = 0; f < TOTAL_FACE_NUM; ++f )
   {
      momt_accum[f] = Fprop[f].DenM_ * Fprop[f].VelM_;
   }
}

// Section: Momentum
void DiscreteProblem::
Compute_Momt_Trans()
{
   for( int f = 0; f < TOTAL_FACE_NUM; ++f )
   {
      int LC = Get_LC_Index( f );
      int RC = Get_RC_Index( f );
      momt_trans[f] = (Cprop[RC].DenM * Cprop[RC].VelM * Cprop[RC].VelM - 
		       Cprop[LC].DenM * Cprop[LC].VelM * Cprop[LC].VelM ) / dX;
   }
}

// Section: Momentum
void DiscreteProblem::
Compute_Momt_Press()
{
   for( int f = 0; f < TOTAL_FACE_NUM; ++f )
   {
      int LC = Get_LC_Index( f );
      int RC = Get_RC_Index( f );
      momt_press[f] = (Cstate[RC].P - Cstate[LC].P) / dX;
   }
}

// Section: Momentum
void DiscreteProblem::
Compute_Momt_Frict()
{
   // Units
   // D        [ m ]
   // density  [ Kg/m3 ]
   // velocity [ m/s ]
   // fM - Moody
   ADs re_( 0.0 );
   ADs fM_( 0.0 );
   ADs VisM_( 0.0 );
   for( int f = 0; f < TOTAL_FACE_NUM; ++f )
   {
      re_ = PropCalc.ReynoldsNumber( Fprop[f].DenM_, Fprop[f].VelM_, Fprop[f].VisM_, D );
      fM_ = PropCalc.Moody_Friction_Factor( re_ );
      momt_frict[f] = fM_ * Fprop[f].DenM_ * Fprop[f].VelM_ * fabs(Fprop[f].VelM_) / (2.0*D);
   }
}

// Section: Momentum
void DiscreteProblem::
Compute_Momt_Gravt()
{
   for( int f = 0; f < TOTAL_FACE_NUM; ++f )
   {
      momt_gravt[f] = -Fprop[f].DenM_ * g * std::cos(theta);
   }
}

// Section: Momentum
void DiscreteProblem::
Compute_Momt_Resid( double dT )
{
   Compute_Momt_Accum();
   Compute_Momt_Trans();
   Compute_Momt_Frict();
   Compute_Momt_Press();
   Compute_Momt_Gravt();

   for( int f = 0; f < TOTAL_FACE_NUM; ++f )
   {
      int Momt_ID = MomtEqnID( f );
      residual[ Momt_ID ] = (momt_accum[f] - momt_accum[f])/dT + momt_trans[f] + momt_press[f]
	 + momt_frict[f] - momt_gravt[f];
   }
}

// Section: Drift-Flux
void DiscreteProblem::
Compute_DFlx_Resid()
{

   ADs c0( 0.0 ), vgj( 0.0 );
   ADs reG( 0.0 ), reL( 0.0 ), re( 0.0 );
   ADs L( 0.0 ), a1( 0.0 ), b1( 0.0 ), k0( 0.0 ), r( 0.0 );
   ADs c1( 0.0 ), c2( 0.0 ), c3( 0.0 );
   ADs c1x( 0.0 ), c5( 0.0 ), c10( 0.0 );
   ADs vsL( 0.0 ), vsG( 0.0 ), vs( 0.0 );
   ADs gama( 0.0 ), jfrx( 0.0 );
   double c4;

   int status;
   double cos_theta = std::cos( theta );
   for( int f = 0; f < TOTAL_FACE_NUM; ++f )
   {
      const ADs& VL = Fstate[f].VL;
      const ADs& VG = Fstate[f].VG;

      const ADs& HL = Fprop[f].H_[ PhaseID::L ];
      const ADs& HG = Fprop[f].H_[ PhaseID::G ];

      const ADs& denL = Fprop[f].Den_[ PhaseID::L ];
      const ADs& denG = Fprop[f].Den_[ PhaseID::G ];

      const ADs& visL = Fprop[f].Vis_[ PhaseID::L ];
      const ADs& visG = Fprop[f].Vis_[ PhaseID::G ];

      if( cos_theta * VL.value() >= 0.0 && cos_theta * VG.value() >= 0.0 )
	 status = 0; // co-down
      else if( cos_theta * VL.value() < 0.0 && cos_theta * VG.value() < 0.0 )
	 status = 1; // co-up
      else
	 status = 2; // counter-current

      reL = epri.ReyNum( denL, VL, visL, D );
      reG = epri.ReyNum( denG, VG, visG, D );

      // < b1 > in all cases for convenience (avoid one more switch-statement)
      switch( status )
      {
	 case 0:    // co-down
	    L   = epri.L_CD( HG );
	    re  = epri.Re_CD( reG );
	    b1  = epri.B1( re );
	    c1x = epri.C1x_CD();
	    break;

	 case 1:    // co-up
	    L   = epri.L_CU( HG );
	    re  = epri.Re_CU( reG, reL );
	    b1  = epri.B1( re );
	    c1x = epri.C1x_CU( b1 );
	    break;

	 case 2:    // counter-current
	    b1 = epri.B1( re );
	    break;
      }

      k0 = epri.K0( denG, denL, b1 );
      r  = epri.R(  denG, denL, b1 );

      c1 = epri.C1( HG,   c1x );
      c2 = epri.C2( denG, denL );
      c4 = epri.C4( D );

      vsL   = HL * VL;
      vsG   = HG * VG;
      vs    = vsL + vsG;

      // Note
      // compute < vgj > first because it's needed in < C0_CD >
      switch( status )
      {
	 case 0:    // co-down
	    jfrx  = epri.Jfrx_CD();
	    gama = epri.Gamma_CD( reG, reL, jfrx, D );
	    c10   = epri.C10( reL, gama, jfrx,  D );
	    c3    = epri.C3_CD( reL, c10 );
	    //
	    vgj = epri.Vgj( denG, denL, c1, c2, c3, c4, D, g, surface_tension );
	    c0  = epri.C0_CD( L, k0, r, HG, vgj, c1, vsG, vsL );
	    break;

	 case 1:    // co-up
	    c3  = epri.C3_CU( reL );
	    //
	    vgj = epri.Vgj( denG, denL, c1, c2, c3, c4, D, g, surface_tension );
	    c0  = epri.C0_CU( L, k0, r, HG );
	    break;

	 case 2:    // counter-current
	    break;
      }

      int DFlx_ID = DFlxEqnID( f );
      residual[ DFlx_ID ] = Fstate[f].VG - c0*vs - vgj;
   }
}

void DiscreteProblem::
Setup_SS()
{
   ss.resize( 2 );
   
   // Top: source
   ss[0].loc             = 0;
   ss[0].is_source       = true;
   ss[0].mode            = SSMode::FixQ;
   ss[0].Q[ PhaseID::G ] = -1.0;
   ss[0].Q[ PhaseID::L ] = -1.0;
   

   // Bottom: sink
   ss[1].loc             = TOTAL_CELL_NUM - 1;
   ss[1].is_source       = false;
   ss[1].mode            = SSMode::FixQ;
   ss[1].Q[ PhaseID::G ] = 1.0;
   ss[1].Q[ PhaseID::L ] = 1.0;
}

void DiscreteProblem::
Compute_SS()
{
   const double dV = 3.1415926 / 4.0 * D * D * dX;
   for( std::size_t i = 0; i < ss.size(); ++i )
   {
      std::size_t loc = ss[i].loc;
      
      if( ss[i].mode == SSMode::FixQ )
      {
	 ADs QL( 0.0 ), QG( 0.0 );
	 QL = Cprop[loc].Den[PhaseID::L] * ss[i].Q[PhaseID::L] / dV;
	 residual[ MassEqnID( loc, PhaseID::L ) ] -= QL;

	 QG = Cprop[loc].Den[PhaseID::G] * ss[i].Q[PhaseID::G] / dV;
	 residual[ MassEqnID( loc, PhaseID::G ) ] -= QG;
      }
      else // FixP
      {
	 // ... ...
      }
   }
}
