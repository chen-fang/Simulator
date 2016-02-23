#pragma once
#include <vector>

typedef std::vector<double> VEC;

class Discretization
{
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

   void Evaluate_Variables( const ADv& VAR );

   ADs Mass_GL( double dX, double dT, const ADv& VAR, ADv& OUTPUT_Residual );

   ADs Momentum_GL( double dX, double dT, const ADv& VAR, ADv& OUTPUT_Residual );

private:
   int PHASE_NUM;
   int TOTAL_CELL_NUM;
   int TOTAL_FACE_NUM;

   // Central variables in staggerd grid scheme
   ADv DenL;
   ADv DenG;
   ADv DenM;
   ADv VL;
   ADv VG;
   ADv VelM2; // VelM^2

   // Facial variables in staggerd grid scheme

   // Central variables needed at boundary surface
   ADv DenL_;
   ADv DenG_;
   ADv DenM_;
   ADv HL_;
   ADv HG_;
   ADv VelM_;

   // Data from previous time step
   // mass
   Vec Hg_Old;
   Vec DenL_Old;
   Vec DenG_Old;
   // momentum
   Vec DenM_Old_;
   Vec VelM_Old_;
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
Evaluate_Variables( const ADv& VAR )
{
   //
   // Central variables: DenL, DenG, DenM
   //
   for( int i = 0; i < TOTAL_CELL_NUM; ++i )
   {
      // naturally (staggered grid)
      int P_idx = Get_P_Index_in_VAR( i );
      int Hg_idx = Get_Hg_Index_in_VAR( i );

      DenL[i] = FluidProperty::Density_Wat( VAR[P_idx] );
      DenG[i] = FluidProperty::Density_Air( VAR[P_idx] );
      DenM[i] = FluidProperty::Density_Mix( DenL[i], DenG[i], VAR[Hg_idx] );
   }

   //
   // Face variables: DenL_, DenG_, DenM_, HL_, HG_, VelM_
   //
   for( int i = 0; i < TOTAL_FACE_NUM; ++i )
   {
      int LC = Get_LC_Index( i );
      int RC = Get_RC_Index( i );

      int HG_idx_LC = Get_HG_Index_in_VAR( LC );
      int HG_idx_RC = Get_HG_Index_in_VAR( RC );

      int VL_idx = Get_VL_Index_in_VAR( LC, RC );
      int VG_idx = Get_VG_Index_in_VAR( LC, RC );

      // upstream
      // liquid: DenL_, HL_
      if( VAR[ VL_idx ].value() >= 0.0 )
      {
	 DenL_[i] = DenL[ LC ];
	 HL_[i] = 1.0 - VAR[ HG_idx_LC ];
      }
      else
      {
	 DenL_[i] = DenL[ RC ];
	 HL_[i] = 1.0 - VAR[ HG_idx_RC ];
      }
      // gas: DenG_, Hg_
      if( VAR[ VG_idx ].value() >= 0.0 )
      {
	 DenG_[i] = DenG[ LC ];
	 HG_[i] = VAR[ HG_idx_LC ];
      }
      else
      {
	 DenG_[i] = DenG[ RC ];
	 HG_[i] = VAR[ HG_idx_RC ];
      }
      // mixture: DenM_, VelM_
      DenM_[i] = DenL_[i] * HL_[i] + DenG_[i] * HG_[i];
      VelM_[i] = ( DenL_[i] * HL_[i] * VAR[VL_idx] + DenG_[i] * HG_[i] * VAR[VG_idx] ) / DenM_[i];
   }

   // Central (dependent on above variables): VL, VG, VelM2
   for( int i = 0; i < TOTAL_CELL_NUM; ++i )
   {
      // VL & VG: leftmost
      if( i == 0 )
      {
	 int VL_idx_RF = Get_VL_Index_in_VAR( i, i+1 );
	 int VG_idx_RF = Get_VG_Index_in_VAR( i, i+1 );

	 VL[i] = VAR[VL_idx_RF] / 2.0;
	 VG[i] = VAR[VG_idx_RF] / 2.0;
      }
      // VL & VG: rightmost
      if( i == TOTAL_CELL_NUM -1 )
      {
	 int VL_idx_LF = Get_VL_Index_in_VAR( i-1, i );
	 int VG_idx_LF = Get_VG_Index_in_VAR( i-1, i );

	 VL[i] = VAR[VL_idx_LF] / 2.0;
	 VG[i] = VAR[VG_idx_LF] / 2.0;
      }
      // VL & VG: other
      if( i != 0 && i != TOTAL_CELL_NUM -1 )
      {
	 int VL_idx_LF = Get_VL_Index_in_VAR( i-1, i );
	 int VG_idx_LF = Get_VG_Index_in_VAR( i-1, i );

	 int VL_idx_RF = Get_VL_Index_in_VAR( i, i+1 );
	 int VG_idx_RF = Get_VG_Index_in_VAR( i, i+1 );

	 // liquid
	 if( VAR[VL_idx_LF].value() > 0.0 && VAR[VL_idx_RF].value() > 0.0 )
	 {
	    VL[i] = VAR[VL_idx_LF];
	 }
	 else if( VAR[VL_idx_LF].value() < 0.0 && VAR[VL_idx_RF].value() < 0.0 )
	 {
	    VL[i] = VAR[VL_idx_RF];
	 }
	 else
	 {
	    VL[i] = (VAR[VL_idx_LF] + VAR[VL_idx_RF]) / 2.0;
	 }
	 // gas
	 if( VAR[VG_idx_LF].value() > 0.0 && VAR[VG_idx_RF].value() > 0.0 )
	 {
	    VG[i] = VAR[VG_idx_LF];
	 }
	 else if( VAR[VG_idx_LF].value() < 0.0 && VAR[VG_idx_RF].value() < 0.0 )
	 {
	    VG[i] = VAR[VG_idx_RF];
	 }
	 else
	 {
	    VG[i] = (VAR[VG_idx_LF] + VAR[VG_idx_RF]) / 2.0;
	 }
      }
   }

   // VelM2
   VelM2[i] = (DenL[i] * HL[i] * VL[i]*VL[i] + DenG[i] * HG[i] * VG[i]*VG[i]) / DenM[i];
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
      ADs R1( 0.0 );
      R1 = (DenM_[i] * VelM_[i] - DenM_Old_[i] * VelM_[i]) / dT;

      int LC = Get_LC_Index( i );
      int RC = Get_RC_Index( i );

      ADs R2( 0.0 );
      R2 = (DenM[RC] * VelM2[RC] - DenM[LC] * VelM2[LC]) / dX;

      ADs R3_P( 0.0 );
      
      int P_idx_LC = Get_P_Index_in_VAR( LC );
      int P_idx_RC = Get_P_Index_in_VAR( RC );

      R3_P = (VAR[P_idx_RC] - VAR[P_idx_LC]) / dX;



      ADs R4_friction( 0.0 );
      //
      //
      //

      ADs R5_G( 0.0 );
      R5_G[i] = DenM_[i] * g * Cos(theta);

      int index_for_Momentum_Mix = Get_VL_Index_in_VAR( LC, RC );
      OUTPUT_Residual[index_for_Momentum_Mix] = R1 + R2 + R3_P + R4_friction - R5_G;
   }
}
