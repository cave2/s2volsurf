#!/bin/csh
#
# build.csh
#
# Copyright 2012 David G. Barnes
#
# This file is part of S2VOLSURF.
#
# S2VOLSURF is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published
# by the Free Software Foundation, either version 3 of the License,
# or (at your option) any later version.
#
# S2VOLSURF is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with S2VOLSURF.  If not, see <http://www.gnu.org/licenses/>. 
#
# We would appreciate it if research outcomes using S2VOLSURF would
# provide the following acknowledgement:
#
# "Three-dimensional visualisation was conducted with the S2PLOT
# progamming library"
#
# and a reference to
#
# D.G.Barnes, C.J.Fluke, P.D.Bourke & O.T.Parry, 2006, Publications
# of the Astronomical Society of Australia, 23(2), 82-93.
#

if (!(${?S2EXTRAINC})) then
  setenv S2EXTRAINC ""
endif
if (!(${?S2EXTRALIB})) then
  setenv S2EXTRALIB ""
endif

setenv S2EXTRAINC "${S2EXTRAINC} -fopenmp -I/opt/local/include"

setenv S2EXTRALIB "${S2EXTRALIB} -fopenmp -lz -lpng"

clbuild.csh xrw libxrw.c
clbuild.csh obj libobj.c

cbuild.csh makexrw
cbuild.csh objrange
cbuild.csh ushortraw2xrw
cbuild.csh 3dcheckerboard2xrw
cbuild.csh tgastack2xrw
cbuild.csh plutodblstack2xrw
cbuild.csh xrw2tgastack
cbuild.csh xrw2xrw
cbuild.csh xrw2pdf
cbuild.csh xrw2points
cbuild.csh xrwinfo
cbuild.csh s2stl
cbuild.csh xrw2stl


if (${?PGPLOT_DIR}) then
    setenv saveinc "$S2EXTRAINC"
    setenv savelib "$S2EXTRALIB"
    setenv S2EXTRAINC "${S2EXTRAINC} -I$PGPLOT_DIR"
    setenv S2EXTRALIB "${S2EXTRALIB} -L$PGPLOT_DIR -lcpgplot -lpgplot -L/opt/local/lib -lX11"
    cbuild.csh xrwhisto
    setenv S2EXTRALIB "$savelib"
    setenv S2EXTRAINC "$saveinc"
else
  echo "PGPLOT_DIR not set, cannot build xrwhist"
endif

setenv S2EXTRAINC "-I/usr/X11/include ${S2EXTRAINC} -I/usr/X11/include"
setenv S2EXTRALIB "${S2EXTRALIB} -L/opt/local/lib"

cbuild.csh xrw2pngmos

setenv S2EXTRAINC "${S2EXTRAINC} -I/usr/local/niftilib/2.0.0/include/nifti -I/opt/local/include/nifti"
setenv S2EXTRALIB "${S2EXTRALIB} -L/usr/local/niftilib/2.0.0/lib -lniftiio -lznz"
cbuild.csh nifti2xrw

setenv S2EXTRALIB "${S2EXTRALIB} -lcfitsio"
cbuild.csh fits2xrw

