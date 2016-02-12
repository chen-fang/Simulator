#pragma once
#include <vector>

typedef std::vector<double> VEC;

class Discretization
{
public:
   Discretization( int _TOTAL_CELL_NUM );

   // Retrive information from [ Den ]
   int Get_DenL_Index( int cell_index );
   int Get_DenG_Index( int cell_index );

   // Retrieve information from [ VAR ]
   //
   int Get_P_Index(  int cell_index, int CELL_VAR_SIZE );
   int Get_Hg_Index( int cell_index, int CELL_VAR_SIZE );

   int Get_VsL_Index( int left_cell_index,  int right_cell_index, int CELL_VAR_SIZE );
   int Get_VsG_Index( int left_cell_index,  int right_cell_index, int CELL_VAR_SIZE );

   bool Is_Flow_Positive( int vel_index_in_VAR, const ADv& VAR );

   int Get_DenIndex_at_Face ( int left_cell_index,  int right_cell_index, bool is_flow_positive );

   void Evaluate_Variables( const ADv& VAR );

   ADs Mass( double timestep, const ADv& Density_Field_old,
	     ADv    Density_Field_new );

private:
   int TOTAL_CELL_NUM;
   int TOTAL_FACE_NUM;

   // Variables evaluated at the surface boundary of CV at the CURRENT time -----------------
   ADv DenL_;
   ADv DenG_;
   ADv DenM_;
   ADv Hg_;
   ADv VelM_;
};

Discretization::Discretization( int _TOTAL_CELL_NUM )
{
   TOTAL_CELL_NUM = _TOTAL_CELL_NUM;
   TOTAL_FACE_NUM = _TOTAL_CELL_NUM - 1;
   
   DenL_.resize( total_interface_num, 0.0 );
   DenG_.resize( total_interface_num, 0.0 );
   DenM_.resize( total_interface_num, 0.0 );
   Hg_.resize( total_interface_num, 0.0 );
   VelM_.resize( total_interface_num, 0.0 );
}

// Hard-coded for 2-phase
int Discretization::
Get_DenL_Index( int cell_index )
{
   return ( cell_index * 2 + 0 );
}

// Hard-coded for 2-phase
int Discretization::
Get_DenG_Index( int cell_index )
{
   return ( cell_index * 2 + 1 );
}

int Discretization::
Get_P_Index( int cell_index, int CELL_VAR_SIZE );
{
   return ( CELL_VAR_SIZE * cell_index + 0 );
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
Get_VsL_Index( int left_cell_index, int right_cell_index, int CELL_VAR_SIZE )
{
   // Assumption: left & right cell indices are valid; no boundary check here.
   //
   return ( CELL_VAR_SIZE * left_cell_index + 2 );
}

int Discretization::
Get_VsG_Index( int left_cell_index, int right_cell_index, int CELL_VAR_SIZE )
{
   // Assumption: left & right cell indices are valid; no boundary check here.
   //
   return ( CELL_VAR_SIZE * left_cell_index + 3 );
}

bool Discretization::
Is_Flow_Positive( int vel_index_in_VAR, const ADv& VAR )
{
   bool is_positive = true;
   if( VAR[ vel_index_in_VAR ].value() < 0.0 )
   {
      is_positive = false;
   }
   return is_positive;
}

int Discretization::
Get_DenIndex_at_Face( int left_cell_index, int right_cell_index, bool is_flow_positive )
{
   // Assumption: left & right cell indices are valid; no boundary check here.
   // Scheme: Upstream
   //
   int den_idx;
   if( is_positive_direction == true ) // flowing from left to right, then take left
   {
      den_idx = left_cell_index;
   }
   else
   {
      den_idx = right_cell_index;
   }
   return den_idx;
}


void Discretization::
Evaluate_Variables( const ADv& VAR )
{
   for( int i = 0; i < TOTAL_FACE_NUM; ++i )
   {
      int cell_idx_left = i;
      int cell_idx_right = i+1;

      // upstream
   }
}


// For two-phase only: gas & liquid
ADs Discretization::
Mass_GL( int i, double dX, double dT, 
	 int CELL_VAR_SIZE, int TOTAL_CELL_NUM,
	 const Vec& Den_Old, const Vec& Hg_Old,
	 const ADv& Den_New, const ADv& VAR, 
	 ADv OUTPUT_Residual_Mass )
{
   ADs R1_L( 0.0 ), R1_G( 0.0 );

   int DenL_idx = Get_DenL_Index( i );
   int DenG_idx = Get_DenG_Index( i );

   int Hg_old_idx = i; // Checkout momentum and decide for old data then come back  here!!!

   int Hg_idx = Get_Hg_Index( i, CELL_VAR_SIZE );
   ADs Hl( 0.0 ); // liquid holdup
   Hl = 1.0 - VAR[ Hg_idx ];

   R1_L = ( Hl * Den_New[ DenL_idx ] - Den_Old[i] ) / dT;

      ADs R2( 0.0 );
      if( i == 0 ) // left boundary
      {
	 int VsL_idx = Get_VsL_Index( i, i+1, CELL_VAR_SIZE );
	 int VsG_idx = Get_VsG_Index( i, i+1, CELL_VAR_SIZE );

	 bool is_VsL_pos = Is_Flow_Positive( VsL_idx, VAR );
	 bool is_VsG_pos = Is_Flow_Positive( VsG_idx, VAR );

	 int denL_idx = Get_DenIndex_at_Face( i, i+1, is_VsL_pos );
	 int denG_idx = Get_DenIndex_at_Face( i, i+1, is_VsG_pos );
	 //
	 R2 = ( Den_New[den_idx] * VAR[vel_idx] - 0.0 ) / dX;
      }
      else if( i == TOTAL_CELL_NUM - 1 ) // right boundary
      {
	 int vel_idx = Get_VelIndex( i-1, i, CELL_VAR_SIZE );
	 bool is_pos = Is_Flow_Positive( vel_idx, VAR );
	 int den_idx = Get_DenIndex_Face( i-1, i, is_pos );
	 //
	 R2 = ( 0.0 - Den_New[dens_idx] * VAR[vel_idx] ) / dX;
      }
      else
      {
	 // left face
	 int left_vel_idx = Get_VelIndex( i-1, i, CELL_VAR_SIZE );
	 bool left_is_pos = Is_Flow_Positive( left_vel_idx, VAR );
	 int left_den_idx = Get_DenIndex( i-1, i, left_is_pos );
	 // right face
	 int right_vel_idx = Get_VelIndex( i, i+1, CELL_VAR_SIZE );
	 bool right_is_pos = Is_Flow_Positive( right_vel_idx, VAR );
	 int right_den_idx = Get_DenIndex( i, i+1, right_is_pos );
	 //
	 R2 = ( Den_New[right_den_idx] * VAR[right_vel_idx] - Den_New[left_den_idx] * VAR[left_vel_idx] ) / dX;
      }

      // Source

      // Assenble
      OUTPUT_Residual_Mass[i] = R1 + R2; // + Source Term
}
