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


#ifndef __SETTINGS_H_INCLUDED__
#define __SETTINGS_H_INCLUDED__

#include "typedefs.hpp"


struct GlobalSettings {
  int s;
  fptype t;
  fptype t_prime;
  fptype U;

  int N_SCC;
  fptype m_prec;
  int max_iterations;
  int init;
  fptype kT;

  int plotmode;
};

GlobalSettings get_precompiled_settings();

#endif //__SETTINGS_H_INCLUDED__
