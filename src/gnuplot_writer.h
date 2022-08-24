/*
 This file is part of FRAME3DD:
 Static and dynamic structural analysis of 2D and 3D frames and trusses with
 elastic and geometric stiffness.
 ---------------------------------------------------------------------------
 http://frame3dd.sourceforge.net/
 ---------------------------------------------------------------------------
 Copyright (C) 1992-2014  Henri P. Gavin

 FRAME3DD is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 FRAME3DD is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with FRAME3DD.  If not, see <http://www.gnu.org/licenses/>.
*//**
	@file
	Gnuplot file writing functions
*/

#include <math.h>

#include "microstran/vec3.h"
#include "types.h"

/*
 * FORCE_BENT_BEAM
 * 	reads internal frame element forces and deflections
 * 	from the internal force and deflection data file.  
 *	Saves deflected shapes to a file.  These bent shapes are exact. 
 */
void force_bent_beam(
	FILE *fpm,	/**< deformed mesh data file pointer	*/
	FILE *fpif,	/**< internal force data file pointer	*/
	char fnif[],	/**< internal force data file name	*/
	int nx,		/**< number of x-axis increments	*/
	int n1, int n2,	/**< node 1 and node 2 of the frame element */
	vec3 *xyz,	/**< node coordinates			*/
	double L,	/**< frame element lengths		*/
	float p,	/**< frame element local rotations	*/
	double *D,	/**< node displacements		*/
	double exagg	/**< mesh exaggeration factor		*/
);


/*
 * CUBIC_BENT_BEAM
 *	computes cubic deflection functions from end deflections
 *	and end rotations.  Saves deflected shapes to a file.
 *	These bent shapes are exact for mode-shapes, and for frames
 *	loaded at their nodes.
 */
void cubic_bent_beam(
	FILE *fpm,	/**< deformed mesh data file pointer	*/
	int n1, int n2,	/**< node 1 and node 2 of the frame element */
	vec3 *xyz,	/**< node coordinates			*/
	double L,	/**< frame element lengths		*/
	float p,	/**< frame element local rotations	*/
	double *D,	/**< node displacements		*/
	double exagg	/**< mesh exaggeration factor		*/
);


/*
 * STATIC_MESH
 *	create mesh data of deformed and undeformed mesh, use gnuplot	22feb99
 *	useful gnuplot options: set noxtics noytics noztics noborder view nokey
 */
void static_mesh(
	char OUT_file[],
	char infcpath[], char meshpath[], char plotpath[],
	char *title, int nN, int nE, int nL, int lc, int DoF,
	vec3 *xyz, double *L,
	int *N1, int *N2, float *p, double *D,
	double exagg_static, int D3_flag, int anlyz, float dx, float scale,
    LoadcaseData *load_cases
);


/*
 * MODAL_MESH
 *	create mesh data of the mode-shape meshes, use gnuplot	19oct98
 *	useful gnuplot options: set noxtics noytics noztics noborder view nokey
 */
void modal_mesh(
	char OUT_file[], char meshpath[], char modepath[],
	char plotpath[], char *title,
	int nN, int nE, int DoF, int nM,
	vec3 *xyz, double *L,
	int *N1, int *N2, float *p,
	double **M, double *f, double **V,
	double exagg_modal, int D3_flag, int anlyz
);


/*
 * ANIMATE
 *	create mesh data of animated mode-shape meshes, use gnuplot	16dec98
 *	useful gnuplot options: set noxtics noytics noztics noborder view nokey
 *	mpeg movie example:   % convert mesh_file-03-f-*.ps mode-03.mpeg
 *	... requires ImageMagick and mpeg2vidcodec packages
 */
void animate(
	char OUT_file[], char meshpath[], char modepath[], char plotpath[],
	char *title,
	int anim[],
	int nN, int nE, int DoF, int nM,
	vec3 *xyz, double *L, float *p,
	int *N1, int *N2, double *f, double **V,
	double exagg_modal, int D3_flag, float pan, float scale
);


