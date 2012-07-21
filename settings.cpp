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


#include "settings.hpp"

GlobalSettings get_precompiled_settings()
{
  GlobalSettings settings;


  // ----------- PHYSICAL SETTINGS -----------

  // size s of the system
  settings.s = 32;

  // hopping parameters and onsite energies
  settings.t = 1.0;
  settings.t_prime = 0.5 * settings.t;
  settings.U = 1.0 * settings.t;

  // ----------- NUMERICAL SETTINGS -----------

  // number of independent self consistency cycles to run
  settings.N_SCC = 8;

  // desired precision for the mean field parameters
  settings.m_prec = 1e-5;

  // maximum number of iterations
  settings.max_iterations = 1000;

  // initialization:
  // 0: random (0,1)
  // 1: checkerboard
  // 2: paramagnetic + initial FD dist
  settings.init = 2;
  settings.kT = 0.25;

  // ----------- OTHER SETTINGS -----------

  // plotting
  // 0: don't plot anything
  // 1: plot final ground state candidate
  // 2: plot everything
  settings.plotmode = 2;


  return settings;
}
