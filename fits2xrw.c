/* fits2xrw.c
 *
 * Copyright 2017 David G. Barnes - based on code by Dany Vohl in encube
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
 * $Id: nifti2xrw.c 224 2014-04-23 22:54:03Z barnesd $
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "fitsio.h" 
#include "libxrw.c"



typedef struct {
   fitsfile *fptr;                      /* Pointer to FITS file */
   long naxes[3];
   float crv[3], crp[3], cde[3];
   float dmin, dmax;
   float low, high;
   float rdmin, rdmax;
   float obsfreq;
   char **label;

   float ***array;

   XYZ min, max;
   XYZ rflag;
   XYZ range, mp;
   long axmax;
   XYZ vp;
} FITSCube;

int errorFITS(int status);
FITSCube readFITScubeHeader(char *fname, int debug);
void readFITScube(FITSCube *cube, int debug);
float ***initVolume(int nx, int ny, int nz, float val);





void usage(char *exename) {
    fprintf(stderr, "usage: %s [options] -f fitsfilename -o xrwfilename\n", exename);
    fprintf(stderr, "options:\n");
    fprintf(stderr, "\t -N\t\t do not normalise\n");
    fprintf(stderr, "\t -R min max\t\t normalise to range min-max\n");
}

int main(int argc, char *argv[]) {

  #define FNAMELEN 400
  char ifname[FNAMELEN+1], ofname[FNAMELEN+1];
  int haveifname = 0, haveofname = 0;
  int normalise = 1;
  int range_normalise = 0;
  float range_min = 0, range_max = 0;

  if (argc < 5) {
    usage(argv[0]);
    exit(-1);
  }

  int ic = 1;
  while (ic < argc) {
    if (!strcmp(argv[ic], "-f")) {
      strncpy(ifname, argv[++ic], FNAMELEN);
      haveifname = 1;
    } else if (!strcmp(argv[ic], "-N")) {
      normalise = 0;
    } else if (!strcmp(argv[ic], "-o")) {
      strncpy(ofname, argv[++ic], FNAMELEN);
      haveofname = 1;
    } else if (!strcmp(argv[ic], "-R")) {
      range_normalise = 1;
      sscanf(argv[++ic], "%f", &range_min);
      sscanf(argv[++ic], "%f", &range_max);
    }
    ic++;
  }
  
  if (!haveifname) {
    usage(argv[0]);
    return -1;
  }

  FILE *fp = fopen(ifname, "r");
  if (!fp) {
    fprintf(stderr, "Couldn't open FITS image (%s)\n", ifname);
    exit(-1);
  }
  fclose(fp);

  FITSCube c = readFITScubeHeader(ifname, 1);
  
  readFITScube(&c, 1);

  int pdid_Nx = c.naxes[0];
  int pdid_Ny = c.naxes[1];
  int pdid_Nz = c.naxes[2];
  
  // normalise cube to be [0,1]
  //float *** ndat = (float ***)malloc(pdid_Nx * sizeof(float **));
  int i, j, k;
  float dmin=9e99, dmax=-9e99;
  
  //float dx, dy, dz;
  for (i = 0; i < pdid_Nx; i++) {
    for (j = 0; j < pdid_Ny; j++) {
      for (k = 0; k < pdid_Nz; k++) {
	if (c.array[i][j][k] < dmin) {
	  dmin = c.array[i][j][k];
	}
	if (c.array[i][j][k] > dmax) {
	  dmax = c.array[i][j][k];
	}
      }
    }
  }
  
  fprintf(stderr, "Normalising data to [0,1] range.\n");
  
  float scale = 1.0 / fabs(dmax-dmin);
  
  for (i=0; i<pdid_Nx; i++) {
    for (j=0; j<pdid_Ny; j++) {
      for (k=0; k<pdid_Nz; k++){
	// in place modification... may be an issue?
	c.array[i][j][k] = (c.array[i][j][k] - dmin) * scale;
	
      }
    }
  }
  
  fprintf(stderr, "HERE\n");
  
  VOL_STRUCT fitsv;
  strcpy(fitsv.filename, "nnnn.nnn");
  
  // input crop
  int blc[3], trc[3];
  blc[0] = blc[1] = blc[2] = 0;
  trc[0] = pdid_Nx - 1;
  trc[1] = pdid_Ny - 1;
  trc[2] = pdid_Nz - 1;

  int stride[3] = {1,1,1};

  int w, h;
  w = pdid_Nx;
  h = pdid_Ny;

  // fixed for input crop
  fitsv.nx = (trc[0] - blc[0] + 1) / stride[0];
  fitsv.ny = (trc[1] - blc[1] + 1) / stride[1];
  fitsv.nz = (trc[2] - blc[2] + 1) / stride[2];

  fitsv.wdx = 1.0 * stride[0];
  fitsv.wdy = 1.0 * stride[1];
  fitsv.wdz = 1.0 * stride[2];
  for (i = 0; i < 256; i++) {
    fitsv.red[i] = fitsv.green[i] = fitsv.blue[i] = powf((float)i / (float)255, 0.0);
  }
  // copy data with stride
  fitsv.data = (float ***)malloc(fitsv.nx * sizeof(float **));
  if (!fitsv.data) {
    fprintf(stderr,"Failed to allocate float volume.\n");
    //return NULL;
    exit(-1);
  }
  int ii, ij, ik, ax, ay, az;
  float sum;

  // fill y and z planes in reverse directions so that home view in 
  // s2plot (and pdf) matches orientation xrw was exported from in 
  // OsiriX - use jtarg, ktarg to control
  int jtarg, ktarg;

  fprintf(stderr, "GOLDFISH\n");
  
  for (i = 0; i < fitsv.nx; i++) {
    ii = blc[0] + i * stride[0];
    // copy this frame in to this slice
    fitsv.data[i] = (float **)malloc(fitsv.ny * sizeof(float *));
    if (!fitsv.data[i]) {
      fprintf(stderr, "Failed to allocate row in float volume.\n");
      exit(-1);
    }
    for (j = 0; j < fitsv.ny; j++) {
      ij = blc[1] + j * stride[1];

      //jtarg = fitsv.ny-1-j;
      jtarg = j;

      fitsv.data[i][jtarg] = (float *)malloc(fitsv.nz * sizeof(float));
      if (!fitsv.data[i][jtarg]) {
        fprintf(stderr, "Failed to allocate column in float volume.\n");
	exit(-1);
      }
      for (k = 0; k < fitsv.nz; k++) {
        ik = blc[2] + k * stride[2];

	//ktarg = fitsv.nz-1-k;
	ktarg = k;

        // averaging
        sum = 0.0;
	/*
	short *deref_s;
	int *deref_i;
	unsigned short *deref_us;
	unsigned char *deref_uc;
        float *deref;
	double *deref64;
	for (ax = 0; ax < stride[0]; ax++) {
          for (ay = 0; ay < stride[1]; ay++) {
            for (az = 0; az < stride[2]; az++) {
	      switch (nim->datatype) {
	      case 2: // unsigned char
		deref_uc = (unsigned char *)nim->data;
		sum += (float)(deref_uc[(ik+az)*w*h + (ij+ay)*w + (ii+ax)]) / 255.0;
		break;
	      case 4: // signed short 
		deref_s = (short *)nim->data;
		sum += (float)(deref_s[(ik+az)*w*h + (ij+ay)*w + (ii+ax)]);
		break;
	      case 8: // integer
		deref_i = (int *)nim->data;
		sum += (float)(deref_i[(ik+az)*w*h + (ij+ay)*w + (ii+ax)]);
		break;
	      case 512: // unsigned short
		deref_us = (unsigned short *)nim->data;
		sum += (float)(deref_us[(ik+az)*w*h + (ij+ay)*w + (ii+ax)]);
		break;
	      case 16: // float
		deref = nim->data;
		sum += (float)(deref[(ik+az)*w*h + (ij+ay)*w + (ii+ax)]);
		break;
	      case 64: // float64
		deref64 = nim->data;
		sum += (double)(deref64[(ik+az)*w*h + (ij+ay)*w + (ii+ax)]);
		break;
		
	      }
            }
          }
        }
	*/

	//sum += (double)c.array[(ik+az)][(ij+ay)][(ii+ax)];
	sum += (double)c.array[(ii+ax)][(ij+ay)][(ik+az)];
	
        fitsv.data[i][jtarg][ktarg] = sum / (float)(stride[0] * stride[1] * stride[2]);

      }
    }
  }

  // normaliseXraw into [0,1]
  if (range_normalise) {
    fprintf(stderr, "RANGE Normalising...\n");
    rangeNormaliseXvol(&fitsv, range_min, range_max);
  } else if (normalise) {
    fprintf(stderr, "Normalising ...\n");
    normaliseXvol(&fitsv);
  }

  XRAW_STRUCT *xraw = Xvol2Xraw(&fitsv);
  if (haveofname) {
    sprintf(xraw->filename, ofname);
  } else {
    sprintf(xraw->filename, "%s.xrw", ifname);
  }
  showXraw(xraw);
  saveXraw(xraw);

  exit(0);

   
  return 1;
}





int errorFITS(int status)
/* Report an error from dealing with FITS files */
{
   if (status) {
       fits_report_error(stderr, status);
       exit(status);
   }
   return 0;
}




FITSCube readFITScubeHeader(char *fname, int debug)
/* Read in a FITS file called fname */
{
   FITSCube cube;                       /* The FITS cube and metadata */
   int status, nfound;
   long i;
   long naxis;

/* NOTE: Should set defaults for all values - assumes FITS metadata exists */
   cube.rdmin   = 0;
   cube.rdmax   = 0;
   cube.low     = 0;
   cube.high    = 0;
   cube.obsfreq = 0;

   status = 0;                          /* Error condition */

   if (fits_open_file(&cube.fptr, fname, READONLY, &status))
      status = errorFITS(status);

   if (fits_read_keys_lng(cube.fptr, "NAXIS", 0, 1, &naxis, &nfound, &status))
      status = errorFITS(status);

#ifdef NEVER
/* CJF: FORCE FOR MANGA */
   if (naxis < 3) {
      fprintf(stderr,"NAXIS = %ld is not a cube!\n", naxis);
      fprintf(stderr,"Check for image extension\n");
      int hdupos,hdutype;
      fits_movabs_hdu(cube.fptr, 2, NULL, &status);
      fits_get_hdu_num(cube.fptr, &hdupos);
      fits_get_hdu_type(cube.fptr, &hdutype, &status); 
      fprintf(stderr,"hdupos = %d %d %d\n",hdupos,hdutype, IMAGE_HDU);
      fprintf(stderr,"status = %d\n",status);
      if (fits_read_keys_lng(cube.fptr, "NAXIS", 0, 1, &naxis, &nfound, &status))
         status = errorFITS(status);
      fprintf(stderr,"Naxis = %ld\n",naxis);
/*
      if (status == END_OF_FILE) status = 0;
      else {
         fits_get_hdu_num(cube.fptr, &hdupos);
         fits_get_hdu_type(cube.fptr, &hdutype, &status); 
         fprintf(stderr,"hdupos = %d %d %d\n",hdupos,hdutype, IMAGE_HDU);
      }
*/
      fits_close_file(cube.fptr, &status);
      exit(1);
   }
#endif

   if (fits_read_keys_lng(cube.fptr, "NAXIS", 1, 3, cube.naxes, &nfound, &status))
      status = errorFITS(status);

   if (fits_read_keys_flt(cube.fptr, "CRVAL1", 0, 1, &cube.crv[0], &nfound, &status))
      status = errorFITS(status);
   if (fits_read_keys_flt(cube.fptr, "CRVAL2", 0, 1, &cube.crv[1], &nfound, &status))
      status = errorFITS(status);
   if (fits_read_keys_flt(cube.fptr, "CRVAL3", 0, 1, &cube.crv[2], &nfound, &status))
      status = errorFITS(status);

/* Coordinates of reference pixel */
   if (fits_read_keys_flt(cube.fptr, "CRPIX1", 0, 1, &cube.crp[0], &nfound, &status))
      status = errorFITS(status);
   if (fits_read_keys_flt(cube.fptr, "CRPIX2", 0, 1, &cube.crp[1], &nfound, &status))
      status = errorFITS(status);
   if (fits_read_keys_flt(cube.fptr, "CRPIX3", 0, 1, &cube.crp[2], &nfound, &status))
      status = errorFITS(status);

/* Gradient at reference pixel */
   if (fits_read_keys_flt(cube.fptr, "CDELT1", 0, 1, &cube.cde[0], &nfound, &status))
      status = errorFITS(status);
   if (fits_read_keys_flt(cube.fptr, "CDELT2", 0, 1, &cube.cde[1], &nfound, &status))
      status = errorFITS(status);
   if (fits_read_keys_flt(cube.fptr, "CDELT3", 0, 1, &cube.cde[2], &nfound, &status))
      status = errorFITS(status);

   cube.label = (char **)calloc(3, sizeof(char *));
   for (i=0;i<3;i++) {
      cube.label[i] = (char *)calloc(256, sizeof(char));
   }

   if (fits_read_keys_str(cube.fptr, "CTYPE1", 0, 1, &cube.label[0], &nfound, &status))
      status = errorFITS(status);
   if (fits_read_keys_str(cube.fptr, "CTYPE2", 0, 1, &cube.label[1], &nfound, &status))
      status = errorFITS(status);
   if (fits_read_keys_str(cube.fptr, "CTYPE3", 0, 1, &cube.label[2], &nfound, &status))
      status = errorFITS(status);

   double MYEPS=1E-7;
   for (i=0;i<3;i++) {
      if (fabs(cube.cde[i]) < MYEPS) cube.cde[i] = 1;
   }


/* Data range */
/* NOTE: Check the FITS reference to see if this is the correct way to calculate these ranges */
   cube.min.x = (1-cube.crp[0])*cube.cde[0] + cube.crv[0];
   cube.max.x = (cube.naxes[0]-cube.crp[0])*cube.cde[0] + cube.crv[0];
   cube.min.y = (1-cube.crp[1])*cube.cde[1] + cube.crv[1];
   cube.max.y = (cube.naxes[1]-cube.crp[1])*cube.cde[1] + cube.crv[1];
   cube.min.z = (1-cube.crp[2])*cube.cde[2] + cube.crv[2];
   cube.max.z = (cube.naxes[2]-cube.crp[2])*cube.cde[2] + cube.crv[2];

   cube.rflag.x = (cube.max.x < cube.min.x) ? -1 : +1;
   cube.rflag.y = (cube.max.y < cube.min.y) ? -1 : +1;
   cube.rflag.z = (cube.max.z < cube.min.z) ? -1 : +1;

   double swap;
   if (cube.rflag.x < 0) {
      swap = cube.min.x; cube.min.x = cube.max.x; cube.max.x = swap;
   }
   if (cube.rflag.y < 0) {
      swap = cube.min.y; cube.min.y = cube.max.y; cube.max.y = swap;
   }
   if (cube.rflag.z < 0) {
      swap = cube.min.z; cube.min.z = cube.max.z; cube.max.z = swap;
   }

/* Data minimum and maximum */
   if (fits_read_keys_flt(cube.fptr, "DATAMIN", 0, 1, &cube.rdmin, &nfound, &status))
      status = errorFITS(status);
   if (fits_read_keys_flt(cube.fptr, "DATAMAX", 0, 1, &cube.rdmax, &nfound, &status))
      status = errorFITS(status);

/* Observing frequency */
/* NOTE: Is this sensible for optical data cubes? Not currently used */
   if (fits_read_keys_flt(cube.fptr, "OBSFREQ", 0, 1, &cube.obsfreq, &nfound, &status))
      status = errorFITS(status);


/* Calculate derived values used for s2plotting */


/* Data range */
   cube.range.x = fabs(cube.max.x - cube.min.x);
   cube.range.y = fabs(cube.max.y - cube.min.y);
   cube.range.z = fabs(cube.max.z - cube.min.z);

/* Mid-point */
   cube.mp.x = (cube.max.x + cube.min.x)*0.5;
   cube.mp.y = (cube.max.y + cube.min.y)*0.5;
   cube.mp.z = (cube.max.z + cube.min.z)*0.5;

/* Biggest axis */
   cube.axmax = cube.naxes[0];
   if (cube.naxes[1] > cube.axmax) cube.axmax = cube.naxes[1];
   if (cube.naxes[2] > cube.axmax) cube.axmax = cube.naxes[2];

/* Viewport range scaled by maximum axis */
   cube.vp.x = (float)cube.naxes[0]/(float)cube.axmax;
   cube.vp.y = (float)cube.naxes[1]/(float)cube.axmax;
   cube.vp.z = (float)cube.naxes[2]/(float)cube.axmax;


   return cube;
}



void readFITScube(FITSCube *cube, int debug)
{
   long fpixel = 1, i,j,k;
   float nullval = 0; /* Don't check for null values in the cube */
   int status = 0, anynull;

   long Mem = (cube->naxes[0]*cube->naxes[1]*cube->naxes[2]*sizeof(float))/(long)1000000;

   if (debug) fprintf(stderr,"Allocating FITScube array: %ld Mbytes\n",Mem);
   float ***array = initVolume(cube->naxes[0],cube->naxes[1],cube->naxes[2],0.0);
   if (array == NULL) {
      fprintf(stderr,"Could not allocate %ld Mbytes\n",Mem);
      exit(1);
   }

   float *dummy;
   Mem = (cube->naxes[0]*sizeof(float))/(long)1000;
   if (debug) fprintf(stderr,"Allocating dummy array: %ld Kbytes\n",Mem);
   dummy = (float *)calloc(cube->naxes[0], sizeof(float));
   if (dummy == NULL) {
      fprintf(stderr,"Could not allocate %ld Kbytes\n",Mem);
      exit(1);
   }

   for (k=0;k<cube->naxes[2];k++) {
      for (j=0;j<cube->naxes[1];j++) {
         if (fits_read_img(cube->fptr, TFLOAT, fpixel,
                cube->naxes[0], &nullval, dummy, &anynull, &status)) {
            status = errorFITS(status);
         }

	/* Copy from dummy array into the FITS cube array */
         for (i=0;i<cube->naxes[0];i++) {
            array[i][j][k] = dummy[i];
         }
         fpixel += cube->naxes[0];	/* Advance the pixel pointer */
      }
   }

/* Confirm min/max from the array */
/* NOTE: Calculated but not used beyond providing a warning */
   cube->dmin = array[0][0][0];
   cube->dmax = array[0][0][0];
   for (i=0;i<cube->naxes[0];i++) {
      for (j=0;j<cube->naxes[1];j++) {
         for (k=0;k<cube->naxes[2];k++) {
            if (array[i][j][k] < cube->dmin) cube->dmin = array[i][j][k];
            if (array[i][j][k] > cube->dmax) cube->dmax = array[i][j][k];
         }
      }
   }
 
   float TEPS = 1.0E-5;
   if (fabs(cube->rdmin-cube->dmin) > TEPS) {
      fprintf(stderr,"Data minimum discrepancy: %12.3f %12.3f\n",
		cube->rdmin,cube->dmin);
   }
   if (fabs(cube->rdmax-cube->dmax) > TEPS) {
      fprintf(stderr,"Data minimum discrepancy: %12.3f %12.3f\n",
		cube->rdmax,cube->dmax);
   } 

   cube->low  = cube->dmin;
   cube->high = cube->dmax;

/* Transfer pointer to array into the FITS cube */
   cube->array = array;
   if (fits_close_file(cube->fptr, &status))
      status = errorFITS(status);

/* Clean up allocated memory that is not longe required */
   if (dummy != NULL) { free(dummy); dummy = NULL; }

   return;
}

float ***initVolume(int nx, int ny, int nz, float val)
/* Allocate memory and initialise a data cube */
{
   float ***volume;
   int i, j;

   float *zero;
   zero = (float *)malloc(nz *sizeof(float));
   for (i=0;i<nz;i++) {
      zero[i] = val;
   }

   volume = (float ***)malloc(nx * sizeof(float **));
   if (volume == NULL) {
      fprintf(stderr,"Failed to allocate %ld bytes\n",nx*sizeof(float **));
      exit(-1);
   }
   for (i=0;i<nx;i++) {
      volume[i] = (float **)malloc(ny * sizeof(float *));
      if (volume[i] == NULL) {
         fprintf(stderr,"Failed to allocate %ld bytes\n",nx*sizeof(float *));
         exit(-1);
      }
      for (j=0;j<ny;j++) {
         volume[i][j] = (float *)malloc(nz * sizeof(float));
         if (volume[i][j] == NULL) {
            fprintf(stderr,"Failed to allocate %ld bytes\n",nx*sizeof(float));
            exit(-1);
         }
         memcpy(volume[i][j], zero, nz*sizeof(float));
                /* Fill with zeroes! */
      }
   }

   if (zero != NULL) { free(zero); zero = NULL; }

   return volume;
}


