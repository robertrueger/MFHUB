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


#include "plot.hpp"

int plot( const GlobalSettings& settings, const SCCResults& results,
          const string& root_dir, int id )
{

  // ----- RESULT OUTPUT TO FILE AND PLOTTING -----

  int const& s = settings.s;

  // make an output directory

  string dir;
  if ( id != -1 ) {
    stringstream tmp;
    tmp << setfill( '0' );
    tmp << "./" << root_dir << "/" << id << '/';
    dir = tmp.str();
    if ( system( ( "test -e " + dir + " || mkdir " + dir ).c_str() ) != 0 ) {
      #pragma omp critical (output)
      { cerr << id << ": ERROR -> unable to create the output directory!"; }
      return 1;
    }
  } else {
    dir = "./" + root_dir + "/";
  }

  ofstream n_log( ( dir + "n.log" ).c_str() );
  if ( !n_log.is_open() ) {
    #pragma omp critical (output)
    { cerr << "ERROR: unable to open occupation output file?" << endl; }
    return 1;
  }
  n_log << setiosflags( ios::scientific );
  n_log.setf( ios::showpos );
  n_log.precision( numeric_limits<fptype>::digits10 + 1 );

#ifdef _VERBOSE
  cout << "Final mean field parameters:" << endl;
#endif
  for ( int i = 0; i < s * s; ++i ) {
#ifdef _VERBOSE
    cout << i << ' ' << idx2x( i, s ) << ' ' << idx2y( i, s ) << ' ' << results.n_up( i ) << ' ' << results.n_down( i ) << endl;
#endif
    n_log << i << ' ' << idx2x( i, s ) << ' ' << idx2y( i, s ) << ' ' << results.n_up( i ) << ' ' << results.n_down( i ) << endl;
  }

  n_log.close();

#ifdef _VERBOSE
  cout << endl << "Plotting ...";
  cout.flush();
#endif

  ofstream gnuplot( ( dir + "plot.gnu" ).c_str() );
  if ( !gnuplot.is_open() ) {
    #pragma omp critical (output)
    { cerr << id << ": ERROR -> unable to create gnuplot script?" << endl; }
    return 1;
  }
  gnuplot << "\
	set terminal pngcairo size 1000,600 \n\
	set size ratio 2/3 \n\
	set xrange [0:" << 1.5 * ( s - 1 ) << "] \n\
	set yrange [0:" << s - 1 << "] \n\
	set tics out \n\
	set cbtics in \n\
	set cbtics \n\
	set dgrid3d " << s * 10 << "," << s * 10 << ",3 \n\
	set pm3d map \n\
	set arrow from 0,0 to " << 0.5 * ( s - 1 ) << "," << s - 1 << " nohead front \n\
	set arrow from " << s - 1 << ",0 to " << 1.5 * ( s - 1 ) << "," << s - 1 << " nohead front \n\
	set output 'm_plot.png' \n\
	set cblabel \"m_z\" \n\
	splot 'n.log' using ($2+0.5*$3):3:($4-$5) notitle \n\
	set output 'n_up_plot.png' \n\
	set cblabel \"n_up\" \n\
	splot 'n.log' using ($2+0.5*$3):3:4 notitle \n\
	set output 'n_down_plot.png' \n\
	set cblabel \"n_down\" \n\
	splot 'n.log' using ($2+0.5*$3):3:5 notitle";
  gnuplot.close();

  if ( system( ( "cd " + dir + " ; gnuplot plot.gnu" ).c_str() ) != 0 ) {
    cerr << id << ": WARNING -> gnuplot call returned exit code != 0" << endl;
  }

#ifdef _VERBOSE
  cout << " done!" << endl;
#endif

  return 0;
}
