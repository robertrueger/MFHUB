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


#ifndef __SCC_INOUT_H_INCLUDED__
#define __SCC_INOUT_H_INCLUDED__

#include <eigen3/Eigen/Core>
using namespace Eigen;

#include "typedefs.hpp"


struct SCCResults {

  int exit_code;

  // convergence information
  bool converged;
  int iterations_to_convergence;
  fptype Delta_n_up, Delta_n_down;

  // physical observables
  fptype energy;
  fptype gap;
  fptype m_z;
  fptype filling;

  // final mean field parameters
  Array<fptype, Dynamic, 1> n_up;
  Array<fptype, Dynamic, 1> n_down;

  // final eigenvalues
  Array<fptype, Dynamic, 1> epsilon_up;
  Array<fptype, Dynamic, 1> epsilon_down;

  // final eigenvectors
  Matrix<fptype, Dynamic, Dynamic> Q_up;
  Matrix<fptype, Dynamic, Dynamic> Q_down;

  SCCResults() : exit_code( 1 ) { }
};

#endif //__SCC_INOUT_H_INCLUDED__
