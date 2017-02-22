/* xrw2tgastack.c
 *
 * Copyright 2017 David G. Barnes
 *
 * This file is part of S2VOLSURF.
 *
 * S2VOLSURF is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published
 * by the Free Software Foundation, either version 3 of the License,
 * or (at your option) any later version.
 *
 * S2VOLSURF is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with S2VOLSURF.  If not, see <http://www.gnu.org/licenses/>. 
 *
 * We would appreciate it if research outcomes using S2VOLSURF would
 * provide the following acknowledgement:
 *
 * "Three-dimensional visualisation was conducted with the S2PLOT
 * progamming library"
 *
 * and a reference to
 *
 * D.G.Barnes, C.J.Fluke, P.D.Bourke & O.T.Parry, 2006, Publications
 * of the Astronomical Society of Australia, 23(2), 82-93.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <time.h>
#include <ctype.h>
#include "libxrw.c"

void usage(char *exename) {
  fprintf(stderr, "usage: %s [options] -f xrwfilename -o outfileformat.tga\n", exename);
  fprintf(stderr, " * * * (outfileformat should include integer specifier e.g. %%04d)\n");
  fprintf(stderr, "-s s1 s2 s3\t\t set data averaging along each axis (1 1 1)\n");
}

int main(int argc, char *argv[]) {

  #define FNAMELEN 400
  char ifname[FNAMELEN+1];
  int haveifname = 0;
  char ofname[FNAMELEN+1];
  int haveofname = 0;

  int stride[3] = {1,1,1}; //{3,3,1};

  if (argc < 5) {
    usage(argv[0]);
    return -1;
  }

  int ic = 1;
  while (ic < argc) {
    if (!strcmp(argv[ic], "-f")) {
      strncpy(ifname, argv[++ic], FNAMELEN);
      haveifname = 1;
    } else if (!strcmp(argv[ic], "-o")) {
      strncpy(ofname, argv[++ic], FNAMELEN);
      haveofname = 1;
    } else if (!strcmp(argv[ic], "-s")) {
      stride[0] = atoi(argv[++ic]);
      stride[1] = atoi(argv[++ic]);
      stride[2] = atoi(argv[++ic]);
      if (stride[0] < 1 || stride[1] < 1 || stride[2] < 1) {
	fprintf(stderr, "stride must be at least 1 in each direction!\n");
	exit(-1);
      }
    }
    ic++;
  }

  if (!haveifname || !haveofname) {
    usage(argv[0]);
    return -1;
  }
    
  XRAW_STRUCT *xr = loadXraw(ifname);

  if (!xr) {
    fprintf(stderr, "Failed to open or read '%s'.\n", argv[1]);
    return -1;
  }

  showXraw(xr);
  fprintf(stdout, "- - - - - - - - - - - - - - - - - - - - - - - - -\n");
  
  VOL_STRUCT *vol = Xraw2Xvol(xr, stride);
  if (!vol) {
    fprintf(stderr, "Failed to parse data volume.\n");
    return -1;
  }

  showXvol(vol);
  
  //writeTGAstack(ofname, vol);
  int nx, ny, nz;
  nx = vol->nx;
  ny = vol->ny;
  nz = vol->nz;

  BITMAP4 *bits; // {unsigned char r,g,b,a}
  bits = (BITMAP4 *)calloc(nx, ny * sizeof(BITMAP4));
  char thisof[FNAMELEN+1];
  float fval;
  unsigned char grey;
  FILE *fout;
  int i, j, k;
  int idx;

  for (k = 0; k < nz; k++) {
    for (i = 0; i < nx; i++) {
      for (j = 0; j < ny; j++) {
	fval = vol->data[i][j][k];
	grey = fval * 255.0;
	idx = i + j * nx;
	bits[idx].r = bits[idx].g = bits[idx].b = grey;
	bits[idx].a = 255;
      }
    }
    sprintf(thisof, ofname, k);
    if (!(k % 10)) {
      fprintf(stderr, "Writing file %s (%d/%d)\n", thisof, k, nz);
    }
    fout = fopen(thisof, "wb");
    Write_Bitmap(fout, bits, nx, ny, 1); // format == 1 => TGA see bitmaplib.c
    fclose(fout);
  }

  return 0;
}

