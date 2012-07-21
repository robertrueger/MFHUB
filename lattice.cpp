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


#include "lattice.hpp"

int xy2idx( int x, int y, const int& s )
{
  // calculate the index from the xy position in the lattice
  x = ( x + s ) % s;
  y = ( y + s ) % s;
  return y * s + x;
}

int idx2x( const int& i, const int& s )
{
  // calculate the x position from the index
  return i % s;
}

int idx2y( const int& i, const int& s )
{
  // calculate the y position from the index
  return i / s;
}
