#pragma once
#include <vector>

typedef std::vector<double> VEC;

class Discretization
{
private:
   int PHASE_NUM;
   int TOTAL_CELL_NUM;
   int TOTAL_FACE_NUM;

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
      CellStateVector                        C;
      FaceStateVector                        F;
   }                                         State;
   
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
      ADs H_Den_[2]; // FaceProp_Upwind
      ADs H_Vis_[2]; // FaceProp_Upwind
      ADs DenM_;     // FaceProp_Upwind
      ADs VelM_;     // FaceProp_Upwind
   }                                        FaceProps;
   typedef std::vector< FaceProps >         FaceVector;


//private:
public:
   // Data from previous time step
   //
   // MeshType mesh;
   // PropertyCalculator PropCalc;
   
   CellVector Cprop;
   FaceVector Fprop;

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



public:
   Discretization( int _PHASE_NUM, int _TOTAL_CELL_NUM );

   // Node and face indices
   // '|' indicates CV face
   // '.' indicates CV nodal points
   // 'X' indicates CV faces that are NOT indexed
   //
   // Node Index: | 0 | 1 | 2 | ... | N-1 |   --- --- Total: N
   // Face Index: X . 0 . 1 . 2 ... N-2 . X   --- --- Total: N-1
   //
   // Left and Right Cell Index
   int Get_LC_Index( int face_index )    { return face_index; }
   int Get_RC_Index( int face_index )    { return face_index + 1; }
   // Left and Right Face Index
   int Get_LF_Index( int cell_index )    { return cell_index - 1; }
   int Get_RF_Index( int cell_index )    { return cell_index;   }

   int MassEqnID( int c, int phaseID )   { return 4 * c + phaseID; }
   int MomtEqnID( int f )                { return 4 * f + 2; }
   int DFlxEqnID( int f )                { return 4 * f + 3; }

   // Section: Property
   void Compute_CellProp_PropCacl( const StateVector& state );
   void Compute_FaceProp_Upwind( const StateVector& state );
   void Compute_CellProp_UpWind()
   void Compute_Properties( const StateVector& state ); // summary

   // Section: Mass
   void Compute_Mass_Accum( const State& state );
   void Compute_Mass_Trans( const State& state );
   void Compute_Mass_Resid( const State& state );

   // Section: Momentum
   void Compute_Momt_Accum();
   void Compute_Momt_Trans();
   void Compute_Momt_Frict();
   void Compute_Momt_Press( const State& state );
   void Compute_Momt_Gravt();
   void Compute_Momt_Resid( const State& state );


};

Discretization::Discretization( int _PHASE_NUM, int _TOTAL_CELL_NUM )
{
   PHASE_NUM = _PHASE_NUM;
   TOTAL_CELL_NUM = _TOTAL_CELL_NUM;
   TOTAL_FACE_NUM = _TOTAL_CELL_NUM - 1;
   //
   DenL.resize( TOTAL_CELL_NUM, 0.0 );
   DenG.resize( TOTAL_CELL_NUM, 0.0 );
   DenM.resize( TOTAL_CELL_NUM, 0.0 );
   //
   DenL_.resize( TOTAL_FACE_NUM, 0.0 );
   DenG_.resize( TOTAL_FACE_NUM, 0.0 );
   DenM_.resize( TOTAL_FACE_NUM, 0.0 );
   HL_.resize(   TOTAL_FACE_NUM, 0.0 );
   HG_.resize(   TOTAL_FACE_NUM, 0.0 );
   VelM_.resize( TOTAL_FACE_NUM, 0.0 );
}


// -----------------------------------------------------
// ------------- Staggered grid scheme -----------------
// -----------------------------------------------------
// Section: Property
void Discretization::
Compute_CellProp_PropCacl( const StateVector& state )
{
   // DenL, DenG, DenM
   for( int c = 0; c < TOTAL_CELL_NUM; ++c )
   {
      Cprop[c].Den[PhaseID::L] = FluidProperty::Density_Wat( state.C[c].P );
      Cprop[c].Den[PhaseID::G] = FluidProperty::Density_Air( state.C[c].P );
      Cprop[c].DenM = FluidProperty::Density_Mix( Cprop[c].Den[PhaseID::L],
						  Cprop[c].Den[PhaseID::G],
						  state.C[c].Hg );
   }
}

// Section: Property
void Discretization::
Compute_FaceProp_Upwind( const StateVector& state )
{
   // H_Den_, DenM_, VelM_
   ADscalar HL( 0.0 );

   for( int f = 0; f < TOTAL_FACE_NUM; ++f )
   {
      int LC = Get_LC_Index( f );
      int RC = Get_RC_Index( f );
      // liquid
      if( state.F[f].VL.value() >= 0.0 )
      {
	 HL = 1.0 - state.C[LC].HG;
	 Fprop[f].H_Den_[PhaseID::L] = HL * Cprop[LC].Den[PhaseID::L];
	 Fprop[f].H_Vis_[PhaseID::L] = HL * FluidProperty::Viscosity_Wat();
      }
      else
      {
	 HL = 1.0 - state.C[RC].HG;
	 Fprop[f].H_Den_[PhaseID::L] = HL * Cprop[RC].Den[PhaseID::L];
	 Fprop[f].H_Vis_[PhaseID::L] = HL * FluidProperty::Viscosity_Wat();
      }
      // gas
      if( state.F[f].VG.value() >= 0.0 )
      {
	 Fprop[f].H_Den_[PhaseID::G] = state.C[LC].HG * Cprop[LC].Den[PhaseID::G];
	 Fprop[f].H_Vis_[PhaseID::G] = state.C[LC].HG * FluidProperty::Viscosity_Air();
      }
      else
      {
	 Fprop[f].H_Den_[PhaseID::G] = state.C[RC].HG * Cprop[RC].Den[PhaseID::G];
	 Fprop[f].H_Vis_[PhaseID::G] = state.C[RC].HG * FluidProperty::Viscosity_Air();
      }

      // mixture: DenM_, VelM_
      Fprop[f].DenM_ = FluidProperty::Density_Mix(   Fprop[f].H_Den_[PhaseID::L],
						     Fprop[f].H_Den_[PhaseID::G] );
      Fprop[f].VisM_ = FluidProperty::Viscosity_Mix( Fprop[f].H_Vis_[PhaseID::L],
						     Fprop[f].H_Vis_[PhaseID::G] );

      Fprop[f].VelM_[i] = ( Fprop[f].H_Den_[PhaseID::L] * state.F[f].VL +
			    Fprop[f].H_Den_[PhaseID::G] * state.F[f].VG ) / Fprop[f].DenM_;
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
	 Cprop[c].VelM = Cprop[RF].VelM_; // / 2.0;
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
Compute_Properties( const StateVector& state )
{
   Compute_CellProp_PropCacl( state );
   Compute_FaceProp_Upwind( state );
   Compute_CellProp_Upwind();
}

// Section: Mass
void Discretization::
void Compute_Mass_Accum( const State& state, double dT )
{
   for( int c = 0; c < TOTAL_CELL_NUM; ++c )
   {
      mass_accumL[c] = (1.0 - state.C[c].HG) * Cprop[c].Den[PhaseID::L] / dT;
      mass_accumG[c] =        state.C[c].HG  * Cprop[c].Den[PhaseID::G] / dT;
   }
}

// Section: Mass
void Discretization::
void Compute_Mass_Trans( const State& state, double dX )
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

      tmpL = Fprop[f].H_Den_[PhaseID::L] * state.F[f].VL / dX;
      tmpG = Fprop[f].H_Den_[PhaseID::G] * state.F[f].VG / dX;

      mass_transL[LC] += tmpL;
      mass_transG[LC] += tmpG;

      mass_transL[RC] -= tmpL;
      mass_transG[RC] -= tmpG;
   }
}

// Section: Mass
void Discretization::
void Compute_Mass_Resid( const State& state )
{
   Compute_Mass_Accum( state );
   Compute_Mass_Trans( state );

   // !!!!!!!!!!!!!!!!!!!!!!!
   // add source term here...

   int MassL_ID, MassG_ID;
   for( int c = 0; c < TOTAL_CELL_NUM; ++c )
   {
      MassL_ID = MassEqnID( c, PhaseID::L );
      MassG_ID = MassEqnID( c, PhaseID::G );
      
      residual[ MassL_ID ] = mass_accumL[c] + mass_transL[c];
      residual[ MassG_ID ] = mass_accumG[c] + mass_transG[c];
   }

   int Cell_ID;
   const double dV = 3.1415926 / 4.0 * D * D / dX;

   // well top: air source at const Q = 5 [m3/s]
   Cell_ID = 0;
   const double Qg_inj = -5.0;
   residual[ MassEqnID( Cell_ID, PhaseID::G ) ] -= Cprop[0].Den[PhaseID::G] * Qg_inj / dV;

   // well bottom: sink
   c = TOTAL_CELL_NUM - 1;
   double Qg = 
   
}


// Section: Momentum
void Discretization::
Compute_Momt_Accum( double dT )
{
   for( int f = 0; f < TOTAL_FACE_NUM; ++f )
   {
      momt_accum[f] = (Fprop[f].DenM_ * Fprop[f].VelM_) / dT;
   }
}

// Section: Momentum
void Discretization::
Compute_Momt_Trans( double dX )
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
Compute_Momt_Press( const State& state, double dX )
{
   for( int f = 0; f < TOTAL_FACE_NUM; ++f )
   {
      int LC = Get_LC_Index( f );
      int RC = Get_RC_Index( f );
      momt_press[f] = (state.C[RC].P - state.C[LC].P) / dX;
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
      re_ = FluidProperty::ReynoldsNumber( Fprop[f].DenM_, Fprop[f].VelM_, Fprop[f].VisM_, D );
      fM_ = FluidProperty::Moody_Friction_Factor( re_ );
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
Compute_Momt_Resid( const State& state )
{
   Compute_Momt_Accum();
   Compute_Momt_Trans();
   Compute_Momt_Frict();
   Compute_Momt_Press( state );
   Compute_Momt_Gravt();

   int Momt_ID, DF_ID;
   for( int f = 0; f < TOTAL_FACE_NUM; ++f )
   {
      Momt_ID = MomtEqnID( f );
      DF_ID   = DflxEqnID( f );

      residual[ Momt_ID ] = 
   }
}
