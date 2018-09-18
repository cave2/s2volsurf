/* libobj.h
 *
 * Copyright 2012 David G. Barnes
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

#define OBJ_STRUCT_LABEL_LEN 80
typedef struct {
  char label[OBJ_STRUCT_LABEL_LEN+1];
  COLOUR col;
  float alpha;
  int nverts;
  XYZ *verts;
  int nfacets;
  int *facets; // 3 * int per facet = 3 * vertex indices
  int *facets_tcs; // 3 * int per facet = 3 * vertex texture coordinate indices
  XYZ minP, maxP; // min,max coordinates
  XYZ meanP; // mean (~central) coordinate
  int nnorms;
  XYZ *norms; // require as many as verts
  int nvtcs; // vertex texture coordinates
  XYZ *vtcs; // use only X,Y (as u,v)
} OBJ_STRUCT;

typedef struct {
  int nverts;
  int *vert_indices;
  float *vert_weights;
  float minW, maxW;
} OBJ_WGT_STRUCT;

OBJ_STRUCT *loadObj(char *filename, char *label, float r, float g, float b, float a);
OBJ_STRUCT *loadObjFromFS(char *filename, char *label, float r, float g, float b, float a);
OBJ_WGT_STRUCT *loadObjWgtFromFS(char *filename);
OBJ_STRUCT *loadObjFromSTL(char *filename, char *label, float r, float g, float b, float a);
void calcObjMinMaxMean(OBJ_STRUCT *obj);
void translateObj(OBJ_STRUCT *obj, XYZ vec);
void transformObj(OBJ_STRUCT *obj, float *m); // m[4c,3r] matrix
void drawObj(OBJ_STRUCT *obj); //, COLOUR c);
void drawWireObj(OBJ_STRUCT *obj);
void drawObjWgts(OBJ_STRUCT *obj, OBJ_WGT_STRUCT *wgt, int wgted_vx_only,
		 int solid_col_idx);
void saveWgtedObj(OBJ_STRUCT *obj, OBJ_STRUCT *obj_norms, OBJ_WGT_STRUCT *wgt, char *fname);


