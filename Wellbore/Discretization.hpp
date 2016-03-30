#pragma once
#include <vector>
#include "Drift-Flux/EPRI.hpp"

typedef std::vector<double> Vec;

class Discretization
{
private:
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
   
   ADv residual;

   FluidProperty PropCalc;

public:
   Discretization( int _PHASE_NUM, int _TOTAL_CELL_NUM, int _dX, int _D );

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
   int DflxEqnID( int f )               { return 4 * f + 3; }

   void Transfer_StateVector( const StateVector& state );
   void Initiate_State( StateVector& state );
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
   void Compute_Dflx_Resid();
};

Discretization::Discretization( int _PHASE_NUM, int _TOTAL_CELL_NUM, int _DX, int _D ) :
   PHASE_NUM( _PHASE_NUM ),
   TOTAL_CELL_NUM( _TOTAL_CELL_NUM ),
   TOTAL_FACE_NUM( _TOTAL_CELL_NUM - 1 ),
   dX( _DX ),
   D( _D ),
   Cprop(  TOTAL_CELL_NUM ),
   Fprop(  TOTAL_FACE_NUM ),
   Cstate( TOTAL_CELL_NUM ),
   Fstate( TOTAL_FACE_NUM )
{

}

void Discretization::
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

void Discretization::
Initiate_State( StateVector& state )
{
   for( int c = 0; c < TOTAL_CELL_NUM; ++c )
   {
      state[c].HG = 0.6;
      state[c].P  = 2.0E+07; // Pa
      //
      state[c].HG.make_independent( MassEqnID(c,PhaseID::L) );
      state[c].P.make_independent(  MassEqnID(c,PhaseID::G) );
   }
   for( int f = 0; f < TOTAL_FACE_NUM; ++f )
   {
      state[f].VL = 0.0;
      state[f].VG = 0.0;
      //
      state[f].VL.make_independent( MomtEqnID(f) );
      state[f].VG.make_independent( DflxEqnID(f) );
   }
}

void Discretization::
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

bool Discretization::
discretize( const StateVector& state, double dT )
{
   bool is_badvalue = false;
   Transfer_StateVector( state );
   Compute_Properties();
   Compute_Mass_Resid( dT );
   Compute_Momt_Resid( dT );
   Compute_Dflx_Resid();

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
      eqn_dflx = DflxEqnID( f );
      if( !std::isfinite( residual[ eqn_momt ].value() ) ) is_badvalue = true;
      if( !std::isfinite( residual[ eqn_dflx ].value() ) ) is_badvalue = true;
   }
   return is_badvalue;
}

template< typename R >
void Discretization::
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
      eqn_dflx = DflxEqnID( f );

      state[f].VL += update[ eqn_momt ];
      state[f].VG += update[ eqn_dflx ];
   }   
}

template< typename V, typename M >
void Discretization::
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
void Discretization::
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
void Discretization::
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
      Fprop[f].H_Den_[PhaseID::L] = HL * Cprop[C].Den[PhaseID::L];
      Fprop[f].H_Vis_[PhaseID::L] = HL * PropCalc.Viscosity_Wat();

      // gas
      if( Fstate[f].VG.value() >= 0.0 )   C = LC;
      else                                C = RC;
      HG = Cstate[C].HG;
      Fprop[f].H_[PhaseID::G]     = HG;
      Fprop[f].Den_[PhaseID::G]   = Cprop[C].Den[PhaseID::G];
      Fprop[f].H_Den_[PhaseID::G] = HG * Cprop[C].Den[PhaseID::G];
      Fprop[f].H_Vis_[PhaseID::G] = HG * PropCalc.Viscosity_Air();

      // mixture: DenM_, VelM_
      Fprop[f].DenM_ = PropCalc.Density_Mix(   Fprop[f].H_Den_[PhaseID::L],
					       Fprop[f].H_Den_[PhaseID::G] );
      Fprop[f].VisM_ = PropCalc.Viscosity_Mix( Fprop[f].H_Vis_[PhaseID::L],
					       Fprop[f].H_Vis_[PhaseID::G] );

      Fprop[f].VelM_ = ( Fprop[f].H_Den_[PhaseID::L] * Fstate[f].VL +
			 Fprop[f].H_Den_[PhaseID::G] * Fstate[f].VG ) / Fprop[f].DenM_;
   }
}

// Section: Property
void Discretization::
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
void Discretization::
Compute_Properties()
{
   Compute_CellProp_PropCacl();
   Compute_FaceProp_Upwind();
   Compute_CellProp_Upwind();
}

// Section: Mass
void Discretization::
Compute_Mass_Accum()
{
   for( int c = 0; c < TOTAL_CELL_NUM; ++c )
   {
      mass_accumL[c] = (1.0 - Cstate[c].HG) * Cprop[c].Den[PhaseID::L];
      mass_accumG[c] =        Cstate[c].HG  * Cprop[c].Den[PhaseID::G];
   }
}

// Section: Mass
void Discretization::
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
void Discretization::
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

   int Cell_ID;
   const double dV = 3.1415926 / 4.0 * D * D * dX;

   // well top: air source at const Q = 5 [m3/s]
   Cell_ID = 0;
   const double Qg_inj = -5.0;
   const double Qw_inj = -2.0;
   residual[ MassEqnID( Cell_ID, PhaseID::L ) ] -= Cprop[0].Den[PhaseID::L] * Qw_inj / dV;
   residual[ MassEqnID( Cell_ID, PhaseID::G ) ] -= Cprop[0].Den[PhaseID::G] * Qg_inj / dV;

   // well bottom: sink
   Cell_ID = TOTAL_CELL_NUM - 1;
   const double Qg_prd = 5.0;
   const double Qw_prd = 2.0;
   residual[ MassEqnID( Cell_ID, PhaseID::L ) ] -= Cprop[0].Den[PhaseID::L] * Qw_inj / dV;
   residual[ MassEqnID( Cell_ID, PhaseID::G ) ] -= Cprop[0].Den[PhaseID::G] * Qg_inj / dV; 
}


// Section: Momentum
void Discretization::
Compute_Momt_Accum()
{
   for( int f = 0; f < TOTAL_FACE_NUM; ++f )
   {
      momt_accum[f] = Fprop[f].DenM_ * Fprop[f].VelM_;
   }
}

// Section: Momentum
void Discretization::
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
void Discretization::
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
void Discretization::
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
      momt_frict[f] = fM * Fprop[f].DenM_ * Fprop[f].VelM_ * abs(Fprop[f].VelM_) / (2.0*D);
   }
}

// Section: Momentum
void Discretization::
Compute_Momt_Gravt()
{
   for( int f = 0; f < TOTAL_FACE_NUM; ++f )
   {
      momt_gravt[f] = -Fprop[f].DenM_ * g * std::cos(theta);
   }
}

// Section: Momentum
void Discretization::
Compute_Momt_Resid( double dT )
{
   Compute_Momt_Accum();
   Compute_Momt_Trans();
   Compute_Momt_Frict();
   Compute_Momt_Press();
   Compute_Momt_Gravt();

   int Momt_ID, DF_ID;
   for( int f = 0; f < TOTAL_FACE_NUM; ++f )
   {
      Momt_ID = MomtEqnID( f );
      residual[ Momt_ID ] = (momt_accum[f] - momt_accum[f])/dT + momt_trans[f] + momt_press[f]
	 + momt_frict[f] - momt_gravt[f];
   }
}

// Section: Drift-Flux
void Discretization::
Compute_Dflx_Resid()
{

   ADs c0( 0.0 ), vgj( 0.0 );
   ADs reG( 0.0 ), reL( 0.0 ), re( 0.0 );
   ADs L( 0.0 ), a1( 0.0 ), b1( 0.0 ), k0( 0.0 ), r( 0.0 );

   // L = L_CD( HG );

   // reG = ReyNum( denG, VG, visG, diameter );
   // reL = ReyNum( denL, VL, visL, diameter );

   // re = Re_CU( reG, reL );
   // a1 = A1( re );
   // b1 = B1( a1 );
   // k0 = K0( denG, denL, b1 );
   // r = R( denG, denL, b1 );

   for( int f = 0; f < TOTAL_FACE_NUM; ++f )
   {
   //    c0 = C0_CD( Fprop[f].H_[PhaseID::G], Fprop[f].Den_[PhaseID::L], Fprop[f].Den_[PhaseID::G],
   // 		  Fstate[f].VL, Fstate[f].VG );
      residual[ DF_ID ] = 1.0;
   }
}
