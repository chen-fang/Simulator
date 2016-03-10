#pragma once

template< typename V, typename M >
void 
DiscreteProblem::extract_R_J( V &r, M &m, std::size_t offset )
{
   mResidual.extract_CSR( r, m.rowptr(), m.colind(), m.value() );
   for ( std::size_t i=0 ; i< m.rowptr().size() ; ++ i ) m.rowptr()[i] += offset;
   for ( std::size_t i=0 ; i< m.colind().size() ; ++ i ) m.colind()[i] += offset;
   m.check_size();
   if ( r.size() != mResidual.size() ) std::cout << "BUG IN JACOBIAN\t zero row found" << r.size() << "!=" << mResidual.size() << std::endl;
   m.ainfo.elliptic_var_id  = 0;
   m.ainfo.n_vars_per_block = 2;
   m.ainfo.n_blocks         = mMesh.size_cells();
   m.ainfo.n_unblocked_vars = 0;
}

template< typename V >
void 
DiscreteProblem::extract_dR_dDT( V &r )
{
   r.resize( mFlow.size() );
   for ( std::size_t i=0 ; i< mFlow.size() ; ++ i ) r[i] = mFlow[i].value();
}

template< typename R >
void 
DiscreteProblem::update_state( StateVector &state, const R & update, bool do_safeguard )
{
   for ( std::size_t c = 0; c < mMesh.size_cells(); ++c )
   {
      state[c].Po  += update[ eqnID(c,PhaseID::O) ];
      if (do_safeguard)
         state[c].Sw  += safeguard_MAC( update[ eqnID(c,PhaseID::W) ] );
      else
         state[c].Sw  += update[ eqnID(c,PhaseID::W) ];
   }
}
   
