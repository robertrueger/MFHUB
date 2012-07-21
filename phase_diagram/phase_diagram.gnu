# Copyright (c) 2012, Robert Rueger <rueger@itp.uni-frankfurt.de>
#
# This file is part of MFHUB.
#
# MFHUB is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# MFHUB is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with MFHUB.  If not, see <http://www.gnu.org/licenses/>.


set xlabel "U/t" 
set ylabel "t'/t"

print "Showing the ground state energy ..."
print "(Press ENTER to continue)"
set zlabel 'ground state energy'
splot 'results.dat' using 4:3:5 notitle
pause -1

print "Showing the gap ..."
print "(Press ENTER to continue)"
set zlabel 'gap'
splot 'results.dat' using 4:3:6 notitle
pause -1
