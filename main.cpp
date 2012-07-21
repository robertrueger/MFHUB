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


#include <iostream>
#include <string>
#include <fstream>
using namespace std;

#include "typedefs.hpp"
#include "settings.hpp"
#include "lattice.hpp"
#include "scc_inout.hpp"
#include "scc_calc.hpp"
#include "plot.hpp"


int main( int argc, char* argv[] )
{

  cout << "HUBBARD MODEL in MEAN FIELD APPROXIMATION" << endl;
  cout << "-----------------------------------------" << endl << endl;

  // properly initialize Eigen for use with OMP
  Eigen::initParallel();

  // initialize c's random number gen (used to generate seeds to the real rngs)
  srand( time( NULL ) );

  // load settings for the simulations ...
  GlobalSettings settings;
  if ( argc != 11 ) {
    cout << "Using precompiled simulation settings ..." << endl;
    settings = get_precompiled_settings();
  } else {
    cout << "Reading the settings from the command line ..." << endl;

    settings.s = atoi( argv[1] );
    settings.t = atof( argv[2] );
    settings.t_prime = atof( argv[3] ) * settings.t;
    settings.U = atof( argv[4] ) * settings.t;

    settings.N_SCC = atoi( argv[5] );

    settings.m_prec = atof( argv[6] ) * settings.s * settings.s;
    settings.max_iterations = atoi( argv[7] );
    
    settings.init = atoi( argv[8] );
    settings.kT = atof( argv[9] );

    settings.plotmode = atoi( argv[10] );
  }

  // prepare the output folder
  string dir;
  {
    stringstream tmp;
    tmp << setfill( '0' );
    tmp << "output_s"         << settings.s
        << "_t"  << setw( 5 ) << int( settings.t * 1000 )
        << "_tp" << setw( 5 ) << int( settings.t_prime * 1000 )
        << "_U"  << setw( 5 ) << int( settings.U * 1000 );
    dir = tmp.str();
    tmp.str() = "";
  }
  if ( system( ( "test -e " + dir + " && rm -r ./" + dir + "/*" ).c_str() ) != 0
       && system( ( "test -e " + dir + " || mkdir " + dir ).c_str() ) != 0 ) {
    cerr << "ERROR: unable to clean/create the output directory!" << endl;
    return 1;
  }

  // "best" simulation results = ground state
  bool some_gsc_found = false; // will be set to true as soon as
                               // one calculation finishes converged ...
  SCCResults gs_candidate;

  // launch N_SCC independent calculations
  #pragma omp parallel for shared(some_gsc_found, gs_candidate) \
                           firstprivate(settings, dir) schedule(dynamic)
  for ( int id = 0; id < settings.N_SCC; ++id ) {

    #pragma omp critical (output)
    { cout << id << ": Calculation started!" << endl; }

    SCCResults results = run_scc( settings, id );

    if ( results.exit_code != 0 ) {
      #pragma omp critical (output)
      { cout << id << ": Calculation failed!" << endl; }
      exit( 1 );
    } else {
      #pragma omp critical (output)
      { cout << id << ": Calculation finished!" << endl; }
      if ( !results.converged ) {
        #pragma omp critical (output)
        { cout << id << ": Calculation did not converge!" << endl; }
      } else {
        #pragma omp critical (output)
        {
          cout << id << ": Calculation converged!" << endl;

          // output simulation results
          cout << id << ": iterations_to_convergence = "
                     << results.iterations_to_convergence << endl;
          cout << id << ": Delta_n_up = " << results.Delta_n_up << endl;
          cout << id << ": Delta_n_down = " << results.Delta_n_down << endl;
          cout << id << ": energy = " << results.energy << endl;
          cout << id << ": gap = " << results.gap << endl;
          cout << id << ": m_z = " << results.m_z << endl;
          cout << id << ": filling = " << results.filling << endl;
        }

        #pragma omp critical (gsupdate)
        {
          // check if this is an improvement over our best estimate of the gs
          if ( !some_gsc_found ||
               ( some_gsc_found && results.energy < gs_candidate.energy ) ) {
            #pragma omp critical (output)
            { cout << id << ": Best estimate of the ground state!" << endl; }
            some_gsc_found = true;
            gs_candidate = results;
          }
        }

        if ( settings.plotmode == 2 ) {
          #pragma omp critical (output)
          { cout << id << ": Plotting started!" << endl; }

          if ( plot( settings, results, dir, id ) != 0 ) {
            #pragma omp critical (output)
            { cerr << id << ": ERROR while plotting the results!" << endl; }
            exit( 1 );
          }
          #pragma omp critical (output)
          {
            cout << id << ": Plotting finished!" << endl;
          }
        }
      }
    }
  }

  // show results on stdout

  cout << endl;
  cout << "All calculations finished!" << endl;
  cout << endl;
  cout << "Best ground state estimate:" << endl;
  cout << "iterations_to_convergence = "
       << gs_candidate.iterations_to_convergence << endl;
  cout << "Delta_n_up = " << gs_candidate.Delta_n_up << endl;
  cout << "Delta_n_down = " << gs_candidate.Delta_n_down << endl;
  cout << "energy = " << gs_candidate.energy << endl;
  cout << "gap = " << gs_candidate.gap << endl;
  cout << "m_z = " << gs_candidate.m_z << endl;
  cout << "filling = " << gs_candidate.filling << endl;

  // output results to file

  cout << "Outputting results in machine readable form to results.log ..." << endl;
  ofstream results_log( ( "./" + dir + "/results.log" ).c_str() );
  if ( !results_log.is_open() ) {
    cerr << "ERROR: unable to open results output file?" << endl;
    return 1;
  }
  results_log << setiosflags( ios::scientific );
  results_log.setf( ios::showpos );
  results_log.precision( numeric_limits<fptype>::digits10 + 1 );

  results_log        << settings.s
              << ' ' << settings.t
              << ' ' << settings.t_prime
              << ' ' << settings.U
              << ' ' << gs_candidate.energy
              << ' ' << gs_candidate.gap
              << ' ' << gs_candidate.m_z
              << ' ' << gs_candidate.filling << endl;

  results_log.close();

  // plot the results

  if ( settings.plotmode >= 1 ) {
    cout << "Plotting ..." << endl;
    if ( plot( settings, gs_candidate, dir ) != 0 ) {
      cerr << "ERROR while plotting the results!" << endl;
    }
    cout << "Plotting finished!" << endl;
  }

  return 0;
}
