/*
 * Copyright (c) 2012, Robert Rueger <rueger@itp.uni-frankfurt.de>
 *
 * This file is part of MFHUB.
 *
 * MFHUB is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * MFHUB is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with MFHUB.  If not, see <http://www.gnu.org/licenses/>.
*/


#include "scc_calc.hpp"

SCCResults run_scc( const GlobalSettings& settings, const int& id )
{

  // ----- INITIALIZATION -----

  SCCResults results;

  // define short names for the most used settings:
  int const& s = settings.s;
  fptype const& t = settings.t;
  fptype const& t_prime = settings.t_prime;
  fptype const& U = settings.U;
  fptype const& m_prec = settings.m_prec;

  // create a new random number generator
  gsl_rng* rng;
  rng = gsl_rng_alloc( gsl_rng_mt19937 );
  gsl_rng_set( rng, rand() );

  // initialize mean field parameter <n_i,sigma>
  Array<fptype, Dynamic, 1> n_up( s * s, 1 );
  Array<fptype, Dynamic, 1> n_down( s * s, 1 );
  if ( settings.init == 0 ) {
    for ( int i = 0; i < s * s; ++i ) {
      n_up( i ) = gsl_rng_uniform_pos( rng );
      n_down( i ) = gsl_rng_uniform_pos( rng );
    }
  } else if ( settings.init == 1 ) {
    for ( int i = 0; i < s * s; ++i ) {
      n_up( i ) =   ( ( i + i / s ) % 2 == 0 ? 1.0 : 0.0 );
      n_down( i ) = ( ( i + i / s ) % 2 == 1 ? 1.0 : 0.0 );
    }
  } else if ( settings.init == 2 ) {
    n_up   = Array<fptype, Dynamic, 1>::Constant( s * s, 1, 0.5 );
    n_down = Array<fptype, Dynamic, 1>::Constant( s * s, 1, 0.5 );
  } else {
    #pragma omp critical (output)
    { cerr << id << ": ERROR -> unknown initialization!" << endl; }
    gsl_rng_free( rng );
    return results;
  }

  // construct the tight-binding part of H_sigma
  // (it doesn't change with the iterations, so
  //  we only need to calculate its matrix once)
  Matrix<fptype, Dynamic, Dynamic> H_tb
                       = Matrix<fptype, Dynamic, Dynamic>::Zero( s * s, s * s );
  for ( int i = 0; i < s * s; ++i ) {
    // calculate the position of atom i in the lattice
    const int x = idx2x( i, s );
    const int y = idx2y( i, s );

    // nearest neighbour hopping
    H_tb( i, xy2idx( x - 1, y, s ) ) -= t;
    H_tb( i, xy2idx( x + 1, y, s ) ) -= t;
    H_tb( i, xy2idx( x, y - 1, s ) ) -= t;
    H_tb( i, xy2idx( x, y + 1, s ) ) -= t;

    // diagonal hopping
    H_tb( i, xy2idx( x - 1, y + 1, s ) ) -= t_prime;
    H_tb( i, xy2idx( x + 1, y - 1, s ) ) -= t_prime;
  }

  // save the old mean field parameters
  Array<fptype, Dynamic, 1> n_up_old = n_up;
  Array<fptype, Dynamic, 1> n_down_old = n_down;

#ifdef _VERBOSE
  cout << endl << "Starting self consistency cycle ..." << endl;

  cout << "Iteration 0: " << endl;
  cout << n_up.transpose().head( 5 ) << endl;
  cout << n_down.transpose().head( 5 ) << endl;
#endif


  // ----- SELF CONSISTENCY CYCLE -----

  // forward declare variables needed in the SCC

  Matrix<fptype, Dynamic, Dynamic> H_up;
  Matrix<fptype, Dynamic, Dynamic> H_down;

  SelfAdjointEigenSolver< Matrix<fptype, Dynamic, Dynamic> > solver_H_up;
  SelfAdjointEigenSolver< Matrix<fptype, Dynamic, Dynamic> > solver_H_down;

  // iteration counter
  int iter = 0;

  do {
    ++iter;

    // construct H_up and H_down from the mean field parameters <n_i,sigma>
    H_up = H_tb;
    H_up += ( U * n_down ).matrix().asDiagonal();
    H_down = H_tb;
    H_down += ( U * n_up ).matrix().asDiagonal();

    // diagonalize H_up and H_down
    solver_H_up.compute( H_up );
    solver_H_down.compute( H_down );
    if ( solver_H_up.info() == NoConvergence ||
         solver_H_down.info() == NoConvergence ) {
      #pragma omp critical (output)
      { cerr << id << ": ERROR -> diagonalization did not converge!" << endl; }
      gsl_rng_free( rng );
      return results;
    }

    // save old mean field parameters
    n_up_old = n_up;
    n_down_old = n_down;



    if ( iter == 1 && settings.init == 2 ) {

      // calculate the fermi energy
      fptype E_fermi = 0.5 * ( solver_H_up.eigenvalues()( ( s * s / 2 ) - 1 ) +
                               solver_H_down.eigenvalues()( ( s * s / 2 ) - 1 ) );

      // create arrays to store which states are occupied
      vector<bool> occupied_up( s * s, false );
      vector<bool> occupied_down( s * s, false );

      // find occupied states according to the fermi distribution
      while ( ( int ) count( occupied_up.begin(), occupied_up.end(), true )
                                                                != s * s / 2 ) {
        for ( int i = 0; i < s * s; ++i ) {
          fptype fdi_up  = fermifunc( solver_H_up.eigenvalues()( i ),
                                                         E_fermi, settings.kT );
          occupied_up[i] = ( fdi_up == 1.0 || gsl_rng_uniform( rng ) < fdi_up );
        }
      }
      while ( ( int ) count( occupied_down.begin(), occupied_down.end(), true )
                                                                != s * s / 2 ) {
        for ( int i = 0; i < s * s; ++i ) {
          fptype fdi_down  = fermifunc( solver_H_up.eigenvalues()( i ),
                                        E_fermi, settings.kT );
          occupied_down[i] = ( fdi_down == 1.0 ||
                               gsl_rng_uniform( rng ) < fdi_down );
        }
      }

#ifdef _VERBOSE
      cout << "Initial occupied states according to FD-statistics:" << endl;
      for ( int i = 0; i < s * s; ++i ) {
        occupied_up[i] ? cout << '1' : cout << '0';
      }
      cout << endl;
      for ( int i = 0; i < s * s; ++i ) {
        occupied_down[i] ? cout << '1' : cout << '0';
      }
      cout << endl << endl;
#endif

      // reset mean field parameters to zero
      n_up   = Array<fptype, Dynamic, 1>::Constant( s * s, 1, 0.0 );
      n_down = Array<fptype, Dynamic, 1>::Constant( s * s, 1, 0.0 );

      // add the contributions of the individual eigenstates
      for ( int alpha = 0; alpha < s * s; ++alpha ) {
        if ( occupied_up[alpha] ) {
          n_up += solver_H_up.eigenvectors().col( alpha ).array().square();
        }
        if ( occupied_down[alpha] ) {
          n_down += solver_H_down.eigenvectors().col( alpha ).array().square();
        }
      }
    } else {
      // update mean field parameters with mixing
      fptype mix = 0.5 * gsl_rng_uniform_pos( rng );

      n_up   = ( 0.25 + mix ) * solver_H_up.eigenvectors()
               .array().block( 0, 0, s * s, s * s / 2 ).square().rowwise().sum()
               + ( 0.75 - mix ) * n_up;
      n_down = ( 0.25 + mix ) * solver_H_down.eigenvectors()
               .array().block( 0, 0, s * s, s * s / 2 ).square().rowwise().sum()
               + ( 0.75 - mix ) * n_down;
    }

#ifdef _VERBOSE
    cout << "Iteration " << iter << ": "
         << ( n_up - n_up_old ).square().sum() << ' '
         << ( n_down - n_down_old ).square().sum() << ' '
         << ( solver_H_up.eigenvalues() + solver_H_down.eigenvalues() )
            .head( s * s / 2 ).sum() << endl;
    cout << n_up.transpose().head( 5 ) << endl;
    cout << n_down.transpose().head( 5 ) << endl;
    cout << endl;
    cout.flush();
#endif

  } while ( ( ( n_up - n_up_old ).array().abs().maxCoeff() > m_prec
              || ( n_down - n_down_old ).array().abs().maxCoeff() > m_prec )
            && iter < settings.max_iterations );

  // delete random number generator
  gsl_rng_free( rng );

#ifdef _VERBOSE
  cout << "Converged after " << iter << " iterations!" << endl << endl;
#endif


  // ----- RESULT OUTPUT -----

  results.converged = ( n_up - n_up_old ).array().abs().maxCoeff() < m_prec
                      && ( n_down - n_down_old ).array().abs().maxCoeff() < m_prec;
  results.iterations_to_convergence = iter;
  results.Delta_n_up = ( n_up - n_up_old ).array().abs().maxCoeff();
  results.Delta_n_down = ( n_down - n_down_old ).array().abs().maxCoeff();

  results.energy = ( solver_H_up.eigenvalues() + solver_H_down.eigenvalues() )
                   .head( s * s / 2 ).sum();
  results.gap = min( solver_H_up.eigenvalues()( ( s * s / 2 ) + 1 )
                                       - solver_H_up.eigenvalues()( s * s / 2 ),
                     solver_H_down.eigenvalues()( ( s * s / 2 ) + 1 )
                                   - solver_H_down.eigenvalues()( s * s / 2 ) );
  results.m_z = n_up.sum() - n_down.sum();
  results.filling =   ( n_up.sum() + n_down.sum() )
                    / static_cast<fptype>( s * s * 2 );

  results.n_up = n_up;
  results.n_down = n_down;
  results.epsilon_up = solver_H_up.eigenvalues();
  results.epsilon_down = solver_H_down.eigenvalues();
  results.Q_up = solver_H_up.eigenvectors();
  results.Q_down = solver_H_down.eigenvectors();

  results.exit_code = 0;
  return results;
}

fptype fermifunc( fptype const& E, fptype const& E_fermi, fptype const& kT )
{
  // the Fermi-Dirac distribution

  if ( kT == 0.0 ) {
    return E <= E_fermi ? 1.0 : 0.0;
  } else {
    return 1.0 / ( exp( ( E - E_fermi ) / kT ) + 1.0 );
  }
}
