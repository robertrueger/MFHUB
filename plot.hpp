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


#ifndef __PLOT_H_INCLUDED__
#define __PLOT_H_INCLUDED__

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
using namespace std;

#include "typedefs.hpp"
#include "settings.hpp"
#include "lattice.hpp"
#include "scc_inout.hpp"


int plot( const GlobalSettings& settings, const SCCResults& results,
          const string& root_dir, int id = -1 );

#endif //__PLOT_H_INCLUDED__
