#!/bin/bash

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


# abort on any errors
set -e

# delete leftover files from former runs
rm -rf $(ls . | grep -v phase_diagram)

# link the mfhub executable
ln -s ../mfhub mfhub

# constant parameters
s=10;
t=1.0;

N_SCC=2000;
m_prec=0.0000001;
max_iterations=100;
init=2;
kT=0.25;

plotmode=1;

# iterate over parameter space
for t_prime in 0.00 0.05 0.10 0.15 0.20 0.25 0.30 0.35 0.40 0.45 0.50 \
                    0.55 0.60 0.65 0.70 0.75 0.80 0.85 0.90 0.95 1.00
do
    for U in 0 0.1 0.2 0.4 0.6 0.8 1.0 1.25 1.5 1.75 2.0 2.5 3 4 5 6 8 10 12 14 16
    do
        ./mfhub $s $t $t_prime $U $N_SCC $m_prec $max_iterations $init $kT $plotmode
    done
done

# collect results
for dir in `ls -l | grep ^d | awk '{print $9}'`; do
	cat ./$dir/results.log >> ./results.raw
    echo >> ./results.raw
done
sed -e 's/inf/0.0000000e+00/g' results.raw | sed -e 's/nan/0.0000000e+00/g' > results.dat
