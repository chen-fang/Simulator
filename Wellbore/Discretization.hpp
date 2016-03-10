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
   }                                         StateVector
   
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
   
   ADv mass_accum;
   Vec mass_accum_old; 
   ADv momt_accum;
   Vec momt_accum_old;
   ADV residual;


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
   int Get_LC_Index( int face_index ) { return face_index; }
   int Get_RC_Index( int face_index ) { return face_index + 1; }
   // Left and Right Face Index
   int Get_LF_Index( int cell_index ) { return cell_index - 1; }
   int Get_RF_Index( int cell_index ) { return cell_index;   }
   // Face Index
   int Get_F_Index( int left_cell_index, int right_cell_index ) { return left_cell_index; }


   // Retrieve information from [ VAR ]
   //
   int Get_P_Index_in_VAR(  int cell_index );
   int Get_HG_Index_in_VAR( int cell_index );

   int Get_VL_Index_in_VAR( int left_cell_index,  int right_cell_index );
   int Get_VG_Index_in_VAR( int left_cell_index,  int right_cell_index );

   void Compute_CellProp_PropCacl( const StateVector& VAR );
   void Compute_FaceProp_Upwind(   const StateVector& VAR );
   void Compute_CellProp_UpWind(   const StateVector& VAR );
   void Compute_Prop( const StateVector& VAR ); // summary

   void Compute_Mass_Accum();
   void Compute_Mass_Trans();
   ADs Mass_GL( double dX, double dT, const ADv& VAR, ADv& OUTPUT_Residual );

   void Compute_Momt_Accum();
   void Compute_Momt_Trans();
   void Compute_Momt_Frict();
   void Compute_Momt_Press();
   void Compute_Momt_Gravt();
   ADs Momentum_GL( double dX, double dT, const ADv& VAR, ADv& OUTPUT_Residual );

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


int Discretization::
Get_P_Index( int cell_index, int CELL_VAR_SIZE );
{
   return ( CELL_VAR_SIZE * cell_index + 0 );
}

int Discretization::
Get_Hg_Index( int cell_index, int CELL_VAR_SIZE )
{
   return ( CELL_VAR_SIZE * cell_index + 1 );
}

int Discretization::
Get_VL_Index( int left_cell_index, int right_cell_index, int CELL_VAR_SIZE )
{
   // Assumption: left & right cell indices are valid; no boundary check here.
   //
   return ( CELL_VAR_SIZE * left_cell_index + 2 );
}

int Discretization::
Get_VG_Index( int left_cell_index, int right_cell_index, int CELL_VAR_SIZE )
{
   // Assumption: left & right cell indices are valid; no boundary check here.
   //
   return ( CELL_VAR_SIZE * left_cell_index + 3 );
}


// -----------------------------------------------------
// ------------- Staggered grid scheme -----------------
// -----------------------------------------------------
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

void Discretization::
Compute_FaceProp_Upwind( const StateVector& state )
{
   // H_Den_, DenM_, VelM_
   ADscalar HL( 0.0 );

   for( int f = 0; f < TOTAL_FACE_NUM; ++f )
   {
      int LC = Get_LC_Index( f );
      int RC = Get_RC_Index( f );

      if( state.F[f].VL.value() >= 0.0 )
      {
	 HL = 1.0 - state.C[LC].HG;
	 Fprop[f].H_Den_[PhaseID::L] = HL * Cprop[LC].Den[PhaseID::L];
      }
      else
      {
	 HL = 1.0 - state.C[RC].HG;
	 Fprop[f].H_Den_[PhaseID::L] = HL * Cprop[RC].Den[PhaseID::L];
      }
      // gas
      if( state.F[f].VG.value() >= 0.0 )
      {
	 Fprop[f].H_Den_[PhaseID::G] = state.C[LC].HG * Cprop[LC].Den[PhaseID::G];
      }
      else
      {
	 Fprop[f].H_Den_[PhaseID::G] = state.C[RC].HG * Cporp[RC].Den[PhaseID::G];
      }

      // mixture: DenM_, VelM_
      Fprop[f].DenM_ = Fprop[f].H_Den_[PhaseID::L] + Fprop[f].H_Den_[PhaseID::G];

      Fprop[f].VelM_[i] = ( Fprop[f].H_Den_[PhaseID::L] * state.C[f].VL +
			    Fprop[f].H_Den_[PhaseID::G] * state.C[f].VG ) / Fprop[f].DenM_;
   }
}

void Discretization::
Compute_CellProp_Upwind( const StateVector& state )
{
   // VelM
   for( int c = 0; c < TOTAL_CELL_NUM; ++c )
   {
      // 1st cell
      if( c == 0 )
      {
	 int RF = Get_RF_Index(c);
	 state.C[c].VelM = state.F[RF].VelM_ / 2.0;
      }
      // last cell
      else if( c == TOTAL_CELL_NUM -1 )
      {
	 int LF = Get_LF_Index(c);
	 state.C[c].VelM = state.F[LF].VelM_ / 2.0;
      }
      // other cells
      else
      {
	 int LF = Get_LF_Index(c);
	 int RF = Get_RF_Index(c);

	 // liquid
	 // < cocurrent right >
	 if( state.F[LF].VelM_.value() > 0.0 && state.F[RF].VelM_.value() > 0.0 )
	 {
	    state.C[c].VelM = state.F[LF].VelM_;
	 }
	 // < cocurrent left >
	 else if( state.F[LF].VelM_.value() > 0.0 && state.F[RF].VelM_.value() > 0.0 )
	 {
	    state.C[c].VelM = state.F[RF].VelM_;
	 }
	 // < counter-current >
	 else
	 {
	    state.C[c].VelM = (state.F[LF].VelM_ + state.F[RF].VelM_) / 2.0;
	 }
      }
   }

   // VelM2
   CpV[c].VelM = ( CpV[c].DenL[c] * CpV[c].HL[c] * VL[c]*VL[c] + DenG[c] * HG[c] * VG[c]*VG[c]) / DenM[c];
}


// Description:
// Discretization of mass equation for two-phase (gas & liquid) system.
// Staggered grid applied.
// Mass equation is discretized at nodal points (center of CV).
ADs Discretization::
Mass_GL( double dX, double dT, const ADv& VAR, ADv& OUTPUT_Residual )
{
   for( int i = 0; i < TOTAL_CELL_NUM; ++i )
   {
      ADs R1_L( 0.0 ), R1_G( 0.0 );

      int HG_idx = Get_HG_Index_in_VAR( i );

      R1_L = ( (1.0 - VAR[HG_idx]) * DenL[i] - (1.0 - HG_Old[i]) * DenL_Old[i] ) / dT;
      R1_G = ( (      VAR[HG_idx]) * DenG[i] - (      HG_Old[i]) * DenG_Old[i] ) / dT;

      //
      ADs R2_L( 0.0 ), R2_G( 0.0 );
      if( i == 0 ) // left boundary
      {
	 int VL_idx_RF = Get_VL_Index_in_VAR( i, i+1 );
	 int VG_idx_RF = Get_VG_Index_in_VAR( i, i+1 );

	 int RF = Get_RF_Index( i );

	 R2_L = ( DenL_[RF] * HL_[RF] * VAR[VL_idx_RF] - 0.0 ) / dX;
	 R2_G = ( DenG_[RF] * HG_[RF] * VAR[VG_idx_RF] - 0.0 ) / dX;
      }
      if( i == TOTAL_CELL_NUM - 1 ) // right boundary
      {
	 int VL_idx_LF = Get_VL_Index_in_VAR( i-1, i );
	 int VG_idx_LF = Get_VG_Index_in_VAR( i-1, i );

	 int LF = Get_LF_Index( i );

	 R2_L = ( 0.0 - DenL_[LF] * HL_[LF] * VAR[VL_idx_LF] ) / dX;
	 R2_R = ( 0.0 - DenG_[LF] * HG_[LF] * VAR[VG_idx_LF] ) / dX;
      }
      if( i != 0 && i != TOTAL_CELL_NUM - 1 )
      {
	 int VL_idx_LF = Get_VL_Index( i-1, i );
	 int VG_idx_LF = Get_VG_Index( i-1, i );
	 int VL_idx_RF = Get_VL_Index( i, i+1 );
	 int VG_idx_RF = Get_VG_Index( i, i+1 );

	 int LF = Get_LF_Index( i );
	 int RF = Get_RF_Index( i );
	 
	 R2_L = (DenL_[RF] * HL_[RF] * VAR[VL_idx_RF] - DenL_[LF] * HL_[LF] * VAR[VL_idx_LF]) / dX;
	 R2_G = (DenG_[RF] * HG_[RF] * VAR[VG_idx_RF] - DenG_[LF] * HG_[LF] * VAR[VG_idx_LF]) / dX;
      }

      // Source

      // Assemble
      int index_for_MassL = Get_P_Index_in_VAR( i );
      int index_for_MassG = Get_Hg_Index_in_VAR( i );
      OUTPUT_Residual[index_for_MassL] = R1_L + R2_L; // + Source Term
      OUTPUT_Residual[index_for_MassG] = R1_G + R2_G; // + Source Term
   }
}


// Description:
// Discretization of momentum equation for two-phase (gas & liquid) system.
// Staggered grid applied.
// Momentum equation is discretized at boundary surface of CV.
ADs Discretization::
Momentum_GL( double dX, double dT, const ADv& VAR, ADv& OUTPUT_Residual )
{
   for( int i = 0; i < TOTAL_FACE_NUM; ++i )
   {
      ADs R1_A( 0.0 );
      R1 = (DenM_[i] * VelM_[i] - DenM_Old_[i] * VelM_Old_[i]) / dT;

      int LC = Get_LC_Index( i );
      int RC = Get_RC_Index( i );

      ADs R2_T( 0.0 );
      R2 = (DenM[RC] * VelM2[RC] - DenM[LC] * VelM2[LC]) / dX;

      ADs R3_P( 0.0 );
      
      int P_idx_LC = Get_P_Index_in_VAR( LC );
      int P_idx_RC = Get_P_Index_in_VAR( RC );

      R3_P = (VAR[P_idx_RC] - VAR[P_idx_LC]) / dX;

      ADs R4_F( 0.0 );
      R4_F = f_M_ * DenM_[i] * VelM_[i] * abs(VelM_[i]) / (2.0 * d);

      ADs R5_G( 0.0 );
      R5_G[i] = DenM_[i] * g * Cos(theta);

      int index_for_Momentum_Mix = Get_VL_Index_in_VAR( LC, RC );
      OUTPUT_Residual[index_for_Momentum_Mix] = R1_A + R2_T + R3_P + R4_F - R5_G;
   }
}
