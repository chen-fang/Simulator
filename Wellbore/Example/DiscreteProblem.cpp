#include "DiscreteProblem.hpp"
#include "MeshBuilder.hpp"
#include "Point.hpp"
#include <fstream>

DiscreteProblem::DiscreteProblem( std::size_t NX, std::size_t NY, std::size_t NZ,
				  double LX, double LY, double LZ,
				  const std::vector<double> &vKX, const std::vector<double> &vKY, const std::vector<double> &vKZ,
				  const std::vector<double> &Phi_ref ) :
   mMesh( MeshBuilder<UniformCartesian> ( NX,NY,NZ,LX,LY,LZ ) ),
   mPropCalc(  ),
   mCells( mMesh.size_cells() ),
   mFaces( mMesh.size_faces() ),
   mWells( ),
   mPhi_ref( Phi_ref ),
   mAccum_old( 2 * mMesh.size_cells() ),
   mAccum( 2 * mMesh.size_cells(), 0.0 ),
   mFlow( 2 * mMesh.size_cells(), 0.0 ),
   mResidual( 2 * mMesh.size_cells(), 0.0 )
{
   initialize_transmissibility( vKX, vKY, vKZ );
F   setup_wells( NX,NY,NZ, LX,LY,LZ, vKX, vKY, vKZ );
};

void
DiscreteProblem::setup_wells(std::size_t NX, std::size_t NY, std::size_t NZ, 
			     double LX, double LY, double LZ,
			     const std::vector<double> &vKX, 
			     const std::vector<double> &vKY, 
			     const std::vector<double> &vKZ  )
{
   mWells.resize( 2 );

   mWells[0].is_producer = false;
   mWells[0].loc  = 0;
   mWells[0].QINJ[PhaseID::W] = 0.75;
   mWells[0].QINJ[PhaseID::O] = 0.0;


   mWells[1].is_producer = true;
   mWells[1].loc = mMesh.size_cells() - 1;
   const double DX = LX/static_cast<double>(NX);
   const double DY = LY/static_cast<double>(NY);
   const double DZ = LZ/static_cast<double>(NZ);
   const double Kh = DZ * sqrt( vKY[mWells[1].loc] * vKX[mWells[1].loc] );
   const double r1 = vKY[mWells[1].loc]/vKX[mWells[1].loc];
   const double r2 = vKX[mWells[1].loc]/vKY[mWells[1].loc];
   const double ro = 0.28 * std::sqrt( std::sqrt(r1)*DX*DX + std::sqrt(r2)*DY*DY )/ ( std::pow(r1,0.25) + std::pow(r2,0.25) );
   mWells[1].WI  = 0.001127 * 6.2832 * Kh / log( ro / 0.5 );
   mWells[1].Pbh = 2450.0;

   std::cout << "KH = " << Kh << std::endl;
   std::cout << "WI = " << mWells[1].WI << std::endl;   
}

void 
DiscreteProblem::initialize_state( DiscreteProblem::StateVector &state )
{
   state.resize( mMesh.size_cells( ) );
   for (std::size_t c = 0; c < mMesh.size_cells(); ++c ) 
   {
      state[c].status = StatusID::OW;
      state[c].Po     = 2500.0;
      state[c].Po.make_independent( eqnID(c,PhaseID::O) ) ;
      state[c].Sw     = 0.25;
      state[c].Sw.make_independent( eqnID(c,PhaseID::W) ) ;
   }
}

void 
DiscreteProblem::bind_to_old_state( const DiscreteProblem::StateVector &old_state )
{
   compute_cell_properties( old_state );
   compute_accumulation( );
   for ( std::size_t c=0; c<mMesh.size_cells(); ++c )
   {
      std::size_t eqn1 = eqnID(c,PhaseID::O);
      std::size_t eqn2 = eqnID(c,PhaseID::W);
      mAccum_old[ eqn1 ] = mAccum[ eqn1 ].value();
      mAccum_old[ eqn2 ] = mAccum[ eqn2 ].value();
   }
}

bool
DiscreteProblem::discretize( const DiscreteProblem::StateVector &state, double DT )
{
   bool is_badvalue = false;
   mDT = DT;
   compute_cell_properties( state );
   compute_accumulation( );
   compute_flow( );
   for ( std::size_t c=0; c<mMesh.size_cells(); ++c )
   {
      std::size_t eqn1 = eqnID(c,PhaseID::O);
      std::size_t eqn2 = eqnID(c,PhaseID::W);
      mResidual[eqn1] = (mAccum[eqn1] - mAccum_old[eqn1]) + mDT * mFlow[eqn1];
      mResidual[eqn2] = (mAccum[eqn2] - mAccum_old[eqn2]) + mDT * mFlow[eqn2];
      if ( !std::isfinite( mResidual[eqn1].value() ) ) is_badvalue = true;
      if ( !std::isfinite( mResidual[eqn2].value() ) ) is_badvalue = true;
   }   
   return is_badvalue;
}

bool
DiscreteProblem::is_converged ( ConvergenceInfo & nrm )
{
   bool is_MATBAL_converged = true;
   bool is_NRMSAT_converged = true;

   double tot_PV = 0.0;
   for ( std::size_t c=0; c<mMesh.size_cells( ); ++c )
      tot_PV += mCells[c].phi.value() * mMesh.cell_measure( c )*0.1781076;
   
   for ( std::size_t phs = 0; phs < 2; ++phs )
   {
      double sum_R = 0.0;
      double avg_B = 0.0;
      double max_R_PV = 0.0;
      for ( std::size_t c=0; c<mMesh.size_cells( ); ++c )
      {
	 sum_R    += mResidual[ eqnID( c, phs ) ].value();
	 avg_B    += 1.0 / mCells[c].b[ phs ].value();
	 double R_PV = std::abs(mResidual[ eqnID( c, phs ) ].value() / (mCells[c].phi.value() * mMesh.cell_measure( c )*0.1781076 ));
	 if (R_PV > max_R_PV) max_R_PV = R_PV;
      }
      avg_B /= mMesh.size_cells( );
      nrm.MatBal[ phs ]  = std::abs( avg_B * sum_R / tot_PV );
      nrm.NormSat[ phs ] = avg_B * max_R_PV;
      if ( nrm.MatBal[ phs ] > 1.0e-7 )  is_MATBAL_converged = false;
      if ( nrm.NormSat[ phs ] > 0.001 )  is_NRMSAT_converged = false;
   }
   return ( is_MATBAL_converged && is_NRMSAT_converged );
}

void 
DiscreteProblem::initialize_transmissibility( const std::vector<double> & KX,
					     const std::vector<double> & KY,
					     const std::vector<double> & KZ )
{
   Point KL, KR;
   for ( std::size_t f=0; f < mMesh.size_faces( ) ; ++f )
   {
      const std::size_t c1 = mMesh.face_adj_cell_1( f );
      const std::size_t c2 = mMesh.face_adj_cell_2( f );
      KL.p[0] = KX[c1];
      KL.p[1] = KY[c1];
      KL.p[2] = KZ[c1];
      KR.p[0] = KX[c2];
      KR.p[1] = KY[c2];
      KR.p[2] = KZ[c2];
      const Point nrml = mMesh.unit_normal( f );
      const double kl = dot_product( nrml, KL );
      const double kr = dot_product( nrml, KR );
      const double DX = norm( mMesh.cell_coord(c1) - mMesh.cell_coord( c2 ) );
      const double A  = mMesh.face_measure( f );
      const double beta = DX*0.5/A*( 1.0/kl + 1.0/kr );
      mFaces[f].T = 0.00112712 / beta;
   }
};

void 
DiscreteProblem::compute_cell_properties( const DiscreteProblem::StateVector &state )
{
   for (std::size_t c = 0; c < mMesh.size_cells(); ++c ) 
   {
      mPropCalc.calculate( mPhi_ref[c], state[c], mCells[c] );
   }
}

void 
DiscreteProblem::compute_accumulation( )
{
   for (std::size_t c = 0; c < mMesh.size_cells(); ++c ) 
   {
      const double VOLUME = mMesh.cell_measure( c )*0.1781076;
      mAccum[ eqnID(c,PhaseID::O) ] = VOLUME * mCells[c].phi * mCells[c].S[PhaseID::O] * mCells[c].b[PhaseID::O];
      mAccum[ eqnID(c,PhaseID::W) ] = VOLUME * mCells[c].phi * mCells[c].S[PhaseID::W] * mCells[c].b[PhaseID::W];
   }
}

void 
DiscreteProblem::compute_flow( )
{
   compute_face_properties( );

   for (std::size_t c = 0; c < (2*mMesh.size_cells()); ++c ) 
      mFlow[c] = 0.0;

   for (std::size_t f = 0; f < mMesh.size_faces(); ++f ) 
   {
      std::size_t c1 = mMesh.face_adj_cell_1( f );
      std::size_t c2 = mMesh.face_adj_cell_2( f );

      mTmpVars[PhaseID::W] = mFaces[f].T * mFaces[f].L[PhaseID::W] * mFaces[f].Pot[PhaseID::W];
      mTmpVars[PhaseID::O] = mFaces[f].T * mFaces[f].L[PhaseID::O] * mFaces[f].Pot[PhaseID::O];
      
      mFlow[ eqnID(c1,PhaseID::O) ] -= mTmpVars[ PhaseID::O ];
      mFlow[ eqnID(c1,PhaseID::W) ] -= mTmpVars[ PhaseID::W ];

      mFlow[ eqnID(c2,PhaseID::O) ] += mTmpVars[ PhaseID::O ];
      mFlow[ eqnID(c2,PhaseID::W) ] += mTmpVars[ PhaseID::W ];
   }
   compute_wells( );
}

void 
DiscreteProblem::compute_face_properties( )
{
   for (std::size_t f = 0; f < mMesh.size_faces(); ++f ) 
   {
      std::size_t c1 = mMesh.face_adj_cell_1( f );
      std::size_t c2 = mMesh.face_adj_cell_2( f );
      double      dz = ( mMesh.cell_coord(c2) - mMesh.cell_coord(c1) ).p[2];

      for ( std::size_t ph=0; ph < 2; ++ph )
      {
	 mFaces[f].Pot[ph] = (mCells[c2].P[ph] - mCells[c1].P[ph]) + 
	    0.00694 * 0.5 * (mCells[c2].Rho[ph] + mCells[c1].Rho[ph]) * dz;
	    
	 mFaces[f].L[ph] = 0.5 * (mCells[c2].bmu[ph] + mCells[c1].bmu[ph]);
	    
	 if ( mFaces[f].Pot[ph].value() > 0.0 )
	    mFaces[f].L[ph] *= mCells[c2].Kr[ph];
	 else
	    mFaces[f].L[ph] *= mCells[c1].Kr[ph];
      }
   }
}

void 
DiscreteProblem::compute_wells( )
{
   for ( std::size_t w = 0; w<mWells.size(); ++w )
   {
      std::size_t c = mWells[w].loc;
      if ( mWells[w].is_producer )
      {
	 mTmpVars[ 0 ] = mCells[c].P[PhaseID::O] - mWells[w].Pbh;
//	 if ( mTmpVars[ 0 ].value() >= 0.0 )
	 {
	    mFlow[eqnID(c,PhaseID::O)] += mWells[w].WI * mTmpVars[0] * mCells[c].Kr[PhaseID::O] * mCells[c].bmu[PhaseID::O];
	    mFlow[eqnID(c,PhaseID::W)] += mWells[w].WI * mTmpVars[0] * mCells[c].Kr[PhaseID::W] * mCells[c].bmu[PhaseID::W];
         }
      }
      else
      {
	 if ( mWells[w].QINJ[PhaseID::W] > 0.0 ) 
	    mFlow[eqnID(c,PhaseID::W)] -= mWells[w].QINJ[PhaseID::W];

	 if ( mWells[w].QINJ[PhaseID::O] > 0.0 ) 
	 {
	    mFlow[eqnID(c,PhaseID::O)] -= mWells[w].QINJ[PhaseID::O];
	 }
      }
   }
}

double
DiscreteProblem::safeguard_MAC( double upd )
{
   if ( std::abs( upd ) > 0.2 )
      return upd/std::abs(upd) * 0.2;
   else
      return upd;

}

std::ostream & operator << ( std::ostream & ostr, 
			     const DiscreteProblem::ConvergenceInfo & _out )
{
   ostr << _out.NormSat[0]  << "\t" 
	<< _out.NormSat[1]  << "\t" 
	<< _out.MatBal[0]  << "\t" 
	<< _out.MatBal[1]  << std::endl; 

   return ostr;
}

