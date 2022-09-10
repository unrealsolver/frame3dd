/*
 This file is part of FRAME3DD:
 Static and dynamic structural analysis of 2D and 3D frames and trusses with
 elastic and geometric stiffness.
 ---------------------------------------------------------------------------
 http://frame3dd.sourceforge.net/
 ---------------------------------------------------------------------------
 Copyright (C) 1992-2015  Henri P. Gavin

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
*/


#include <stdlib.h>

#include "HPGmatrix.h"
#include "HPGutil.h"
#include "NRutil.h"
#include "common.h"
#include "coordtrans.h"
#include "gnuplot_writer.h"

/*
 * FORCE_BENT_BEAM  -  reads internal frame element forces and deflections
 * from the internal force and deflection data file.
 * Saves deflected shapes to a file.  These bent shapes are exact.
 * Note: It would not be difficult to adapt this function to plot
 * internal axial force, shear force, torques, or bending moments.
 * 9 Jan 2010
 */
void force_bent_beam(
	FILE *fpm, FILE *fpif, char fnif[], int nx, int n1, int n2, vec3 *xyz,
	double L, float p, double *D, double exagg
){
	double	t1, t2, t3, t4, t5, t6, t7, t8, t9; 	/* coord xfmn	*/
	double	xi, dX, dY, dZ;
	float	x, Nx, Vy, Vz, Tx, My, Mz, Dx, Dy, Dz, Rx;
	double	Lx, Ly, Lz;
	int	n;
	int	sfrv=0;		/* *scanf return value	*/

	Lx = xyz[n2].x - xyz[n1].x;
	Ly = xyz[n2].y - xyz[n1].y;
	Lz = xyz[n2].z - xyz[n1].z;

	coord_trans ( xyz, L, n1, n2,
			&t1, &t2, &t3, &t4, &t5, &t6, &t7, &t8, &t9, p );

	x = -1.0;
	n = 0;
	for ( xi = 0; xi <= 1.01*L && n < nx; xi += 0.10*L ) {

		while ( x < xi && n < nx ) {
		    /* read the deformed shape in local coordinates */
		    sfrv=fscanf(fpif,"%f %f %f %f %f %f %f %f %f %f %f",
			&x, &Nx, &Vy, &Vz, &Tx, &My, &Mz, &Dx, &Dy, &Dz, &Rx );
//		    printf("x = %12.4f\n", x );		/* debug */
		    if (sfrv != 11) sferr(fnif);
		    ++n;
		}

		/* exaggerated deformed shape in global coordinates */
		dX = exagg * ( t1*Dx + t4*Dy + t7*Dz );
		dY = exagg * ( t2*Dx + t5*Dy + t8*Dz );
		dZ = exagg * ( t3*Dx + t6*Dy + t9*Dz );

		fprintf (fpm," %12.4e %12.4e %12.4e\n",
					xyz[n1].x + (x/L)*Lx + dX ,
					xyz[n1].y + (x/L)*Ly + dY ,
					xyz[n1].z + (x/L)*Lz + dZ );

//		printf("...  x = %7.3f  n = %3d  Dx = %10.3e   Dy = %10.3e   Dz = %10.3e \n", x,n,Dx,Dy,Dz ); /* debug */
//		printf("                           dX = %10.3e   dY = %10.3e   dZ = %10.3e \n", dX,dY,dZ ); /* debug */

	}

	fprintf(fpm,"\n\n");

	return;
}

/*
 * CUBIC_BENT_BEAM  -  computes cubic deflection functions from end deflections
 * and end rotations.  Saves deflected shapes to a file.  These bent shapes
 * are exact for mode-shapes, and for frames loaded at their nodes.
 * 15 May 2009
 */
void cubic_bent_beam(
	FILE *fpm, int n1, int n2, vec3 *xyz,
	double L, float p, double *D, double exagg
){
	double	t1, t2, t3, t4, t5, t6, t7, t8, t9, 	/* coord xfmn	*/
		u1, u2, u3, u4, u5, u6, u7, u8, u9, u10, u11, u12,
		*a, *b, **A,
		s, v, w, dX, dY, dZ;
	int	i1, i2, pd;
	char	errMsg[MAXL];

	A = dmatrix(1,4,1,4);
	a = dvector(1,4);
	b = dvector(1,4);

	coord_trans ( xyz, L, n1, n2,
			&t1, &t2, &t3, &t4, &t5, &t6, &t7, &t8, &t9, p );

	i1 = 6*(n1-1);	i2 = 6*(n2-1);

		/* compute end deflections in local coordinates */

	u1  = exagg*(t1*D[i1+1] + t2*D[i1+2] + t3*D[i1+3]);
	u2  = exagg*(t4*D[i1+1] + t5*D[i1+2] + t6*D[i1+3]);
	u3  = exagg*(t7*D[i1+1] + t8*D[i1+2] + t9*D[i1+3]);

	u4  = exagg*(t1*D[i1+4] + t2*D[i1+5] + t3*D[i1+6]);
	u5  = exagg*(t4*D[i1+4] + t5*D[i1+5] + t6*D[i1+6]);
	u6  = exagg*(t7*D[i1+4] + t8*D[i1+5] + t9*D[i1+6]);

	u7  = exagg*(t1*D[i2+1] + t2*D[i2+2] + t3*D[i2+3]);
	u8  = exagg*(t4*D[i2+1] + t5*D[i2+2] + t6*D[i2+3]);
	u9  = exagg*(t7*D[i2+1] + t8*D[i2+2] + t9*D[i2+3]);

	u10 = exagg*(t1*D[i2+4] + t2*D[i2+5] + t3*D[i2+6]);
	u11 = exagg*(t4*D[i2+4] + t5*D[i2+5] + t6*D[i2+6]);
	u12 = exagg*(t7*D[i2+4] + t8*D[i2+5] + t9*D[i2+6]);

		/* curve-fitting problem for a cubic polynomial */

	a[1] =  u2;		b[1] =  u3;
	a[2] =  u8;   		b[2] =  u9;
	a[3] =  u6;		b[3] = -u5;
	a[4] =  u12;		b[4] = -u11;

	u7 += L;
	A[1][1] = 1.0;   A[1][2] = u1;   A[1][3] = u1*u1;   A[1][4] = u1*u1*u1;
	A[2][1] = 1.0;   A[2][2] = u7;   A[2][3] = u7*u7;   A[2][4] = u7*u7*u7;
	A[3][1] = 0.0;   A[3][2] = 1.;   A[3][3] = 2.*u1;   A[3][4] = 3.*u1*u1;
	A[4][1] = 0.0;   A[4][2] = 1.;   A[4][3] = 2.*u7;   A[4][4] = 3.*u7*u7;
	u7 -= L;

	lu_dcmp ( A, 4, a, 1, 1, &pd );		/* solve for cubic coef's */

	if (!pd) {
	 sprintf(errMsg," n1 = %d  n2 = %d  L = %e  u7 = %e \n", n1,n2,L,u7);
	 errorMsg(errMsg);
	 exit(30);
	}

	lu_dcmp ( A, 4, b, 0, 1, &pd );		/* solve for cubic coef's */

	// debug ... if deformed mesh exageration is too big, some elements
	// may not be plotted.
	//fprintf( fpm, "# u1=%e  L+u7=%e, dx = %e \n",
	//				u1, fabs(L+u7), fabs(L+u7-u1)/10.0);
	for ( s = u1; fabs(s) <= 1.01*fabs(L+u7); s += fabs(L+u7-u1) / 10.0 ) {

			/* deformed shape in local coordinates */
		v = a[1] + a[2]*s + a[3]*s*s + a[4]*s*s*s;
		w = b[1] + b[2]*s + b[3]*s*s + b[4]*s*s*s;

			/* deformed shape in global coordinates */
		dX = t1*s + t4*v + t7*w;
		dY = t2*s + t5*v + t8*w;
		dZ = t3*s + t6*v + t9*w;

		fprintf (fpm," %12.4e %12.4e %12.4e\n",
			xyz[n1].x + dX , xyz[n1].y + dY , xyz[n1].z + dZ );
	}
	fprintf(fpm,"\n\n");

	free_dmatrix(A,1,4,1,4);
	free_dvector(a,1,4);
	free_dvector(b,1,4);

	return;
}


/* Part of static_mesh */
inline void write_header(
	FILE *fpm, char *title,
	vec3 *xyz, int *N1, int *N2,
	float scale, char D3, time_t now, const Frame *frame
) {
	int j, m, n1, n2 = 0;
	double	mx, my, mz; /* coordinates of the frame element number labels */
	fprintf(fpm,"# FRAME3DD ANALYSIS RESULTS  http://frame3dd.sf.net/");
	fprintf(fpm," VERSION %s \n", VERSION);
	fprintf(fpm,"# %s\n", title );
	fprintf(fpm,"# %s", ctime(&now) );
	fprintf(fpm,"# G N U P L O T   S C R I P T   F I L E \n");
	/* fprintf(fpm,"#  X=%d , Y=%d , Z=%d, D3=%d  \n", X,Y,Z,D3_flag); */

	fprintf(fpm,"set autoscale\n");
	fprintf(fpm,"unset border\n");
	fprintf(fpm,"set pointsize 1.0\n");
	fprintf(fpm,"set xtics; set ytics; set ztics; \n");
	fprintf(fpm,"unset zeroaxis\n");
	fprintf(fpm,"unset key\n");
	fprintf(fpm,"unset label\n");
	fprintf(fpm,"set size ratio -1    # 1:1 2D axis scaling \n");
	fprintf(fpm,"# set view equal xyz # 1:1 3D axis scaling \n");

	fprintf(fpm,"# NODE NUMBER LABELS\n");
	for (j = 0; j < frame->nodes.size; j++) {
		Node *node = &frame->nodes.data[j];
		fprintf(
			fpm, "set label ' %d' at %12.4e, %12.4e, %12.4e\n",
			j + 1, node->position.x, node->position.y, node->position.z
		);
	}

	fprintf(fpm,"# ELEMENT NUMBER LABELS\n");
	for (m=1; m <= frame->edges.size; m++) {
		n1 = N1[m];	n2 = N2[m];
		mx = 0.5 * ( xyz[n1].x + xyz[n2].x );
		my = 0.5 * ( xyz[n1].y + xyz[n2].y );
		mz = 0.5 * ( xyz[n1].z + xyz[n2].z );
		fprintf(
			fpm,"set label ' %d' at %12.4e, %12.4e, %12.4e\n",
			m, mx, my, mz );
	}

	// 3D plot setup commands
	fprintf(fpm,"%c set parametric\n", D3 );
	fprintf(fpm,"%c set view 60, 70, %5.2f \n", D3, scale );
	fprintf(fpm,"%c set view equal xyz # 1:1 3D axis scaling \n", D3 );
	fprintf(fpm,"%c unset key\n", D3 );
	fprintf(fpm,"%c set xlabel 'x'\n", D3 );
	fprintf(fpm,"%c set ylabel 'y'\n", D3 );
	fprintf(fpm,"%c set zlabel 'z'\n", D3 );
	//fprintf(fpm,"%c unset label\n", D3 );
}

/* Part of static_mesh */
inline void writer_undeformed_mesh(
	char meshpath[], FILE *fpm, char *title,
	vec3 *xyz, int *N1, int *N2,
	time_t now, char *errMsg, Frame *frame
) {
	int m, n;
	// open the undeformed mesh data file for writing
	if ((fpm = fopen (meshpath, "w")) == NULL) {
		sprintf(errMsg,"\n  error: cannot open gnuplot undeformed mesh data file: %s\n", meshpath );
		errorMsg(errMsg);
		exit(21);
	}

	fprintf(fpm,"# FRAME3DD ANALYSIS RESULTS  http://frame3dd.sf.net/");
	fprintf(fpm," VERSION %s \n", VERSION);
	fprintf(fpm,"# %s\n", title );
	fprintf(fpm,"# %s", ctime(&now) );
	fprintf(fpm,"# U N D E F O R M E D   M E S H   D A T A   (global coordinates)\n");
	fprintf(fpm,"# Node        X            Y            Z \n");

	for (m=1; m <= frame->edges.size; m++) {
		n = N1[m];	// i = 6*(n-1);
		fprintf (fpm,"%5d %12.4e %12.4e %12.4e \n",
					n , xyz[n].x , xyz[n].y , xyz[n].z );
		n = N2[m];	// i = 6*(n-1);
		fprintf (fpm,"%5d %12.4e %12.4e %12.4e",
					n , xyz[n].x , xyz[n].y , xyz[n].z );
		fprintf (fpm,"\n\n\n");
	}
	fclose(fpm);
}

/*
 * STATIC_MESH  - create mesh data of deformed and undeformed mesh  22 Feb 1999
 * use gnuplot
 * useful gnuplot options: unset xtics ytics ztics border view key
 * This function illustrates how to read the internal force output data file.
 * The internal force output data file contains all the information required
 * to plot deformed meshes, internal axial force, internal shear force, internal
 * torsion, and internal bending moment diagrams.
 */
void static_mesh(
		char IN_file[],
		char infcpath[], char meshpath[], char plotpath[],
		char *title, int lc, int DoF,
		vec3 *xyz, double *L,
		int *N1, int *N2, float *p, double *D,
		double exagg_static, int D3_flag, int anlyz, float dx, float scale,
		LoadCases *load_cases,
		Frame *frame
){
	// Hookup the old variable names
	const int nN = frame->nodes.size;
	const int nE = frame->edges.size;

	FILE	*fpif=NULL, *fpm=NULL;
	char	fnif[FILENMAX], meshfl[FILENMAX],
		D2='#', D3='#',	/* indicates plotting in 2D or 3D	*/
		errMsg[MAXL],
		ch = 'a';
	int	sfrv=0,		/* *scanf return value			*/
		frel, nx,	/* frame element number, number of increments */
		n1, n2;		/* node numbers			*/
	float	x1, y1, z1,	/* coordinates of node n1		*/
		x2, y2, z2;	/* coordinates of node n2		*/
	int	j=0, m=0, n=0,
		X=0, Y=0, Z=0,
		lw = 1;		/* line width of deformed mesh		*/
	time_t  now;		/* modern time variable type		*/

	(void) time(&now);

	// write gnuplot plotting script commands

	for ( j=1; j<=nN; j++ ) { // check for three-dimensional frame
		if (xyz[j].x != 0.0) X=1;
		if (xyz[j].y != 0.0) Y=1;
		if (xyz[j].z != 0.0) Z=1;
	}
	if ( (X && Y && Z) || D3_flag ) {
		D3 = ' '; D2 = '#';
	} else {
		D3 = '#'; D2 = ' ';
	}

	if (lc <= 1) {	// open plotting script file for writing
		if ((fpm = fopen (plotpath, "w")) == NULL) {
			sprintf (errMsg,"\n  error: cannot open gnuplot script file: %s \n", plotpath);
			errorMsg(errMsg);
			exit(23);
		}
	} else {	// open plotting script file for appending
		if ((fpm = fopen (plotpath, "a")) == NULL) {
			sprintf (errMsg,"\n  error: cannot open gnuplot script file: %s \n", plotpath);
			errorMsg(errMsg);
			exit(24);
		}
	}

	// file name for deformed mesh data for load case "lc"
	if ( lc >= 1 && anlyz )	{
		sprintf( meshfl, "%sf.%03d", meshpath, lc );
	}

	// write header, plot-setup cmds, node label, and element label data

	if (lc <= 1) {	// header & node number & element number labels
		write_header(fpm, title, xyz, N1, N2, scale, D3, now, frame);
	}

	// different plot title for each load case
	fprintf(fpm,"set title \"%s\\n", title );
	fprintf(fpm,"analysis file: %s ", IN_file );
	if ( anlyz ) {
		fprintf(fpm,"  deflection exaggeration: %.1f ", exagg_static );
		fprintf(fpm,"  load case %d of %d \"\n", lc, load_cases->size );
	} else {
		fprintf(fpm,"  data check only \"\n");
	}
	fprintf(fpm,"unset clip; \nset clip one; set clip two\n");
	fprintf(fpm,"set xyplane 0 \n"); // requires Gnuplot >= 4.6
	fprintf(fpm, "unset arrow\n");

	if (lc == 0) {
		fprintf(stderr, "Panic. load_case index is negative\n");
		// TODO code?
		exit(255);
	}

	const LoadCase load_case = load_cases->data[lc - 1];

	// Draw static force arrows
	for (unsigned i = 0; i < load_case.loads.point.size; i++) {
		PointLoad point_load = load_case.loads.point.data[i];
		Node node = frame->nodes.data[point_load.node_id];
		fprintf(
			fpm, "set arrow %d from %f, %f, %f rto %f, %f, %f lw 2\n",
			i + 1,
			node.position.x, node.position.y, node.position.z,
			point_load.force.x, point_load.force.y, point_load.force.z
		);
	}

	// 2D plot command
	fprintf(fpm,"%c plot '%s' u 2:3 t 'undeformed mesh' w lp ", D2, meshpath);
	if (!anlyz) fprintf(fpm,"lw %d lt 1 pt 6 \n", lw );
	else fprintf(
		fpm,"lw 1 lt 5 pt 6, '%s' u 1:2 t 'load case %d of %d' w l lw %d lt 3\n",
		meshfl, lc, load_cases->size, lw
	);

	// 3D plot command
	fprintf(fpm,"%c splot '%s' u 2:3:4 t 'load case %d of %d' w lp ", D3, meshpath, lc, load_cases->size);
	if (!anlyz) fprintf(fpm," lw %d lt 1 pt 6 \n", lw );
	else fprintf(
		fpm, " lw 1 lt 5 pt 6, '%s' u 1:2:3 t 'load case %d of %d' w l lw %d lt 3\n",
		meshfl, lc, load_cases->size, lw
	);

	if ( lc < load_cases->size && anlyz ) fprintf(fpm, "pause -1\n");

	fclose(fpm);

	// write undeformed mesh data
	if (lc <= 1) {
		writer_undeformed_mesh(meshpath, fpm, title, xyz, N1, N2, now, errMsg, frame);
	}

	if (!anlyz) return; // no deformed mesh

	// write deformed mesh data

	// open the deformed mesh data file for writing
	if ((fpm = fopen (meshfl, "w")) == NULL) {
		sprintf (errMsg,"\n  error: cannot open gnuplot deformed mesh data file %s \n", meshfl );
		errorMsg(errMsg);
		exit(22);
	}

	fprintf(fpm,"# FRAME3DD ANALYSIS RESULTS  http://frame3dd.sf.net/");
	fprintf(fpm," VERSION %s \n", VERSION);
	fprintf(fpm,"# %s\n", title );
	fprintf(fpm,"# L O A D  C A S E   %d  of   %d \n", lc, load_cases->size);
	fprintf(fpm,"# %s", ctime(&now) );
	fprintf(fpm,"# D E F O R M E D   M E S H   D A T A ");
	fprintf(fpm,"  deflection exaggeration: %.1f\n", exagg_static );
	fprintf(fpm,"#       X-dsp        Y-dsp        Z-dsp\n");

	// open the interior force data file for reading
	if ( dx > 0.0 && anlyz ) {
		// file name for internal force data for load case "lc"
		sprintf( fnif, "%s%02d", infcpath, lc );
		if ((fpif = fopen (fnif, "r")) == NULL) {
			sprintf (errMsg,"\n  error: cannot open interior force data file: %s \n",fnif);
			errorMsg(errMsg);
			exit(20);
		}
	}

	for (m=1; m<=nE; m++) { // write deformed shape data for each element
		ch = 'a';

		fprintf( fpm, "\n# element %5d \n", m );
		if ( dx < 0.0 && anlyz ) {
			cubic_bent_beam ( fpm,
				N1[m],N2[m], xyz, L[m],p[m], D, exagg_static );
		}
		if ( dx > 0.0 && anlyz ) {
			while ( ch != '@' )	ch = getc(fpif);
			sfrv=fscanf(fpif,"%d %d %d %f %f %f %f %f %f %d",
			 &frel, &n1, &n2, &x1, &y1, &z1, &x2, &y2, &z2, &nx);
			if (sfrv != 10) sferr(fnif);
			if ( frel != m || N1[m] != n1 || N2[m] != n2 ) {
			 fprintf(stderr," error in static_mesh parsing\n");
			 fprintf(stderr,"  frel = %d; m = %d; nx = %d \n", frel,m,nx );
			}
			/* debugging ... check mesh data
			printf("  frel = %3d; m = %3d; n1 =%4d; n2 = %4d; nx = %3d L = %f \n", frel,m,n1,n2,nx,L[m] );
			*/
			while ( ch != '~' )	ch = getc(fpif);
			force_bent_beam ( fpm, fpif, fnif, nx,
				N1[m],N2[m], xyz, L[m],p[m], D, exagg_static );
		}

	}

	if ( dx > 0.0 && anlyz ) fclose(fpif);

	fclose(fpm);

	return;
}


/*
 * MODAL_MESH  -  create mesh data of the mode-shape meshes, use gnuplot	19oct98
 * useful gnuplot options: unset xtics ytics ztics border view key
 */
void modal_mesh(
		char IN_file[], char meshpath[], char modepath[],
		char plotpath[], char *title,
		int DoF, int nM,
		vec3 *xyz, double *L,
		int *J1, int *J2, float *p,
		double **M, double *f, double **V,
		double exagg_modal, int D3_flag, int anlyz, Frame *frame
){
	// Hookup the old variable names
	const int nN = frame->nodes.size;
	const int nE = frame->edges.size;

	FILE	*fpm;
	double mpfX, mpfY, mpfZ;	/* mode participation factors	*/
	double *msX, *msY, *msZ;
	double *v;		/* a mode-shape vector */

	int	i, j, m,n, X=0, Y=0, Z=0;
	int	lw = 1;		/*  line thickness of deformed mesh	*/
	char	D2='#', D3 = '#',	/* indicate 2D or 3D frame	*/
		modefl[FILENMAX],
		errMsg[MAXL];


	msX = dvector(1,DoF);
	msY = dvector(1,DoF);
	msZ = dvector(1,DoF);
	v   = dvector(1,DoF);

	for (i=1; i<=DoF; i++) {	/* modal participation factors */
		msX[i] = msY[i] = msZ[i] = 0.0;
		for (j=1; j<=DoF; j+=6) msX[i] += M[i][j];
		for (j=2; j<=DoF; j+=6) msY[i] += M[i][j];
		for (j=3; j<=DoF; j+=6) msZ[i] += M[i][j];
	}

	if (!anlyz) exagg_modal = 0.0;

	for (m=1; m<=nM; m++) {

		sprintf( modefl,"%s-%02d-", modepath, m );

		if ((fpm = fopen (modefl, "w")) == NULL) {
			sprintf (errMsg,"\n  error: cannot open gnuplot modal mesh file: %s \n", modefl);
			errorMsg(errMsg);
			exit(27);
		}

		fprintf(fpm,"# FRAME3DD ANALYSIS RESULTS  http://frame3dd.sf.net/");
		fprintf(fpm," VERSION %s \n", VERSION);
		fprintf(fpm,"# %s\n", title );
		fprintf(fpm,"# M O D E   S H A P E   D A T A   F O R   M O D E");
		fprintf(fpm,"   %d\t(global coordinates)\n", m );
		fprintf(fpm,"# deflection exaggeration: %.1f\n\n", exagg_modal );
		mpfX = 0.0;	for (i=1; i<=DoF; i++)    mpfX += V[i][m]*msX[i];
		mpfY = 0.0;	for (i=1; i<=DoF; i++)    mpfY += V[i][m]*msY[i];
		mpfZ = 0.0;	for (i=1; i<=DoF; i++)    mpfZ += V[i][m]*msZ[i];
		fprintf(fpm,"# MODE %5d:   f= %lf Hz, T= %lf sec\n", m,f[m],1./f[m]);
		fprintf(fpm,"#\t\tX- modal participation factor = %12.4e \n", mpfX);
		fprintf(fpm,"#\t\tY- modal participation factor = %12.4e \n", mpfY);
		fprintf(fpm,"#\t\tZ- modal participation factor = %12.4e \n", mpfZ);

		for(i=1; i<=DoF; i++)	v[i] = V[i][m];

		fprintf(fpm,"#      X-dsp       Y-dsp       Z-dsp\n\n");

		for(n=1; n<=nE; n++) {
			fprintf( fpm, "\n# element %5d \n", n );
			cubic_bent_beam ( fpm, J1[n], J2[n], xyz, L[n], p[n], v, exagg_modal );
		}

		fclose(fpm);

		for ( j=1; j<=nN; j++ ) { // check for three-dimensional frame
			if (xyz[j].x != 0.0) X=1;
			if (xyz[j].y != 0.0) Y=1;
			if (xyz[j].z != 0.0) Z=1;
		}

		if ( (X && Y && Z) || D3_flag ) {
			D3 = ' '; D2 = '#';
		} else {
			D3 = '#'; D2 = ' ';
		}


		if ((fpm = fopen (plotpath, "a")) == NULL) {
			sprintf (errMsg,"\n  error: cannot append gnuplot script file: %s \n",plotpath);
			errorMsg(errMsg);
			exit(25);
		}

		fprintf(fpm,"pause -1\n");

		if (m==1) {
			fprintf(fpm,"unset label\n");
			fprintf(fpm,"%c unset key\n", D3 );
		}

		fprintf(fpm,"set title '%s     mode %d     %lf Hz'\n",IN_file,m,f[m]);
		fprintf(fpm, "unset arrow\n");

		// 2D plot command

		fprintf(fpm,"%c plot '%s' u 2:3 t 'undeformed mesh' w l ", D2, meshpath );
		if (!anlyz) fprintf(fpm," lw %d lt 1 \n", lw );
		else fprintf(fpm," lw 1 lt 5 , '%s' u 1:2 t 'mode-shape %d' w l lw %d lt 3\n", modefl, m, lw );

		// 3D plot command

		fprintf(fpm,"%c splot '%s' u 2:3:4 t 'undeformed mesh' w l ",
								D3, meshpath);
		if (!anlyz) fprintf(fpm," lw %d lt 1 \n", lw );
		else fprintf(fpm," lw 1 lt 5 , '%s' u 1:2:3 t 'mode-shape %d' w l lw %d lt 3\n", modefl, m, lw );

		fclose(fpm);

	}

	free_dvector(msX,1,DoF);
	free_dvector(msY,1,DoF);
	free_dvector(msZ,1,DoF);
	free_dvector(v,1,DoF);
}


/*
 * ANIMATE -  create mesh data of animated mode-shape meshes, use gnuplot	16dec98
 * useful gnuplot options: unset xtics ytics ztics border view key
 * mpeg movie example:   % convert mesh_file-03-f-*.ps mode-03.mpeg
 * ... requires ImageMagick and mpeg2vidcodec packages
 */
void animate(
	char IN_file[], char meshpath[], char modepath[], char plotpath[],
	char *title,
	int anim[],
	int DoF, int nM,
	vec3 *xyz, double *L, float *p,
	int *J1, int *J2, double *f, double **V,
	double exagg_modal, int D3_flag,
	float pan,		/* pan rate for animation	     */
	float scale,		/* inital zoom scale in 3D animation */
	Frame *frame
){
	// Hookup the old variable names
	const int nN = frame->nodes.size;
	const int nE = frame->edges.size;

	FILE	*fpm;

	float	x_min = 0.0, x_max = 0.0,
		y_min = 0.0, y_max = 0.0,
		z_min = 0.0, z_max = 0.0,
		Dxyz = 0.0, 		/* "diameter" of the structure	*/
		rot_x_init  =  70.0,	/* inital x-rotation in 3D animation */
		rot_x_final =  60.0,	/* final  x-rotation in 3D animation */
		rot_z_init  = 100.0,	/* inital z-rotation in 3D animation */
		rot_z_final = 120.0,	/* final  z-rotation in 3D animation */
		zoom_init  = 1.0*scale,	/* init.  zoom scale in 3D animation */
		zoom_final = 1.1*scale, /* final  zoom scale in 3D animation */
		frames = 25;		/* number of frames in animation */

	double	ex=10,		/* an exageration factor, for animation */
		*v;

	int	fr, i,j, m,n, X=0, Y=0, Z=0, c, CYCLES=3,
		frame_number = 0,
		lw = 1,		/*  line thickness of deformed mesh	*/
		total_frames;	/* total number of frames in animation	*/

	char	D2 = '#', D3 = '#',	/* indicate 2D or 3D frame	*/
		Movie = '#',	/* use '#' for no-movie  -OR-  ' ' for movie */
		modefl[FILENMAX], framefl[FILENMAX];
	char	errMsg[MAXL];

	for (j=1; j<=nN; j++) {		// check for three-dimensional frame
		if (xyz[j].x != 0.0) X=1;
		if (xyz[j].y != 0.0) Y=1;
		if (xyz[j].z != 0.0) Z=1;
		if (j==1) {
			x_min = x_max = xyz[j].x;
			y_min = y_max = xyz[j].y;
			z_min = z_max = xyz[j].z;
		}
		if (xyz[j].x < x_min ) x_min = xyz[j].x;
		if (xyz[j].y < y_min ) y_min = xyz[j].y;
		if (xyz[j].z < z_min ) z_min = xyz[j].z;
		if ( x_max < xyz[j].x ) x_max = xyz[j].x;
		if ( y_max < xyz[j].y ) y_max = xyz[j].y;
		if ( z_max < xyz[j].z ) z_max = xyz[j].z;
	}
	if ( (X && Y && Z) || D3_flag ) {
		D3 = ' '; D2 = '#';
	} else {
		D3 = '#'; D2 = ' ';
	}

	Dxyz = sqrt( (x_max-x_min)*(x_max-x_min) + (y_max-y_min)*(y_max-y_min) + (z_max-z_min)*(z_max-z_min) );


	if ((fpm = fopen (plotpath, "a")) == NULL) {
		sprintf (errMsg,"\n  error: cannot append gnuplot script file: %s \n",plotpath);
		errorMsg(errMsg);
		exit(26);
	}
	i = 1;
	while ( (m = anim[i]) != 0 && i < 100) {
	 if ( i==1 ) {

	   fprintf(fpm,"\n# --- M O D E   S H A P E   A N I M A T I O N ---\n");
	   fprintf(fpm,"# rot_x_init  = %7.2f\n", rot_x_init );
	   fprintf(fpm,"# rot_x_final = %7.2f\n", rot_x_final );
	   fprintf(fpm,"# rot_z_init  = %7.2f\n", rot_z_init );
	   fprintf(fpm,"# rot_z_final = %7.2f\n", rot_z_final );
	   fprintf(fpm,"# zoom_init   = %7.2f\n", zoom_init );
	   fprintf(fpm,"# zoom_final  = %7.2f\n", zoom_init );
	   fprintf(fpm,"# pan rate    = %7.2f \n", pan );
	   fprintf(fpm,"set autoscale\n");
	   fprintf(fpm,"unset border\n");
	   fprintf(fpm,"%c unset xlabel \n", D3 );
	   fprintf(fpm,"%c unset ylabel \n", D3 );
	   fprintf(fpm,"%c unset zlabel \n", D3 );
	   fprintf(fpm,"%c unset label \n", D3 );
	   fprintf(fpm,"unset key\n");
	   fprintf(fpm,"%c set parametric\n", D3 );

	   fprintf(fpm,"# x_min = %12.5e     x_max = %12.5e \n", x_min, x_max);
	   fprintf(fpm,"# y_min = %12.5e     y_max = %12.5e \n", y_min, y_max);
	   fprintf(fpm,"# z_min = %12.5e     z_max = %12.5e \n", z_min, z_max);
	   fprintf(fpm,"# Dxyz = %12.5e \n", Dxyz );
	   fprintf(fpm,"set xrange [ %lf : %lf ] \n",
			x_min-0.2*Dxyz, x_max+0.1*Dxyz );
	   fprintf(fpm,"set yrange [ %lf : %lf ] \n",
			y_min-0.2*Dxyz, y_max+0.1*Dxyz );
	   fprintf(fpm,"set zrange [ %lf : %lf ] \n",
			z_min-0.2*Dxyz, z_max+0.1*Dxyz );

/*
 *	   if ( x_min != x_max )
 *		fprintf(fpm,"set xrange [ %lf : %lf ] \n",
 *	 		x_min-0.2*(x_max-x_min), x_max+0.2*(x_max-x_min) );
 *	   else fprintf(fpm,"set xrange [ %lf : %lf ] \n",
 *			x_min-exagg_modal, x_max+exagg_modal );
 *	   if (y_min != y_max)
 *		fprintf(fpm,"set yrange [ %lf : %lf ] \n",
 *	 		y_min-0.2*(y_max-y_min), y_max+0.2*(y_max-y_min) );
 *	   else fprintf(fpm,"set yrange [ %lf : %lf ] \n",
 *			y_min-exagg_modal, y_max+exagg_modal );
 *	   if (z_min != z_max)
 *	   	fprintf(fpm,"set zrange [ %lf : %lf ] \n",
 *			z_min-0.2*(z_max-z_min), z_max+0.2*(z_max-z_min) );
 *	   else fprintf(fpm,"set zrange [ %lf : %lf ] \n",
 *			z_min-exagg_modal, z_max+exagg_modal );
 */

	   fprintf(fpm,"unset xzeroaxis; unset yzeroaxis; unset zzeroaxis\n");
	   fprintf(fpm,"unset xtics; unset ytics; unset ztics; \n");
	   fprintf(fpm,"%c set view 60, 70, %5.2f \n", D3, scale );
	   fprintf(fpm,"set size ratio -1    # 1:1 2D axis scaling \n");
	   fprintf(fpm,"%c set view equal xyz # 1:1 3D axis scaling \n", D3 );

	 }

	 fprintf(fpm,"pause -1 \n");
	 fprintf(fpm,"set title '%s     mode %d      %lf Hz'\n",IN_file,m,f[m]);
	fprintf(fpm, "unset arrow\n");

	 frame_number = 0;
	 total_frames = 2*CYCLES*frames;
	 for ( c=1; c <= CYCLES; c++ ) {

	  for ( fr=0; fr<=frames; fr++ ) {

	    ++frame_number;

	    sprintf(modefl,"%s-%02d.%03d", modepath, m, fr  );
	    sprintf(framefl,"%s-%02d-f-%03d.ps", modepath, m, fr  );

	    fprintf(fpm,"%c plot '%s' u 2:3 w l lw 1 lt 5, ", D2,meshpath );
	    fprintf(fpm," '%s' u 1:2 w l lw %d lt 3 ; \n", modefl, lw );
	    if ( pan != 0.0 )
	     fprintf(fpm,"%c set view %7.2f, %7.2f, %5.3f # pan = %f\n", D3,
		rot_x_init + pan*(rot_x_final-rot_x_init)*frame_number/total_frames,
		rot_z_init + pan*(rot_z_final-rot_z_init)*frame_number/total_frames,
		zoom_init + pan*(zoom_final-zoom_init)*frame_number/total_frames, pan );
	    fprintf(fpm,"%c splot '%s' u 2:3:4 w l lw 1 lt 5, ",D3,meshpath);
            fprintf(fpm," '%s' u 1:2:3 w l lw %d lt 3;", modefl, lw );

	    if ( fr==0 && c==1 )	fprintf(fpm,"  pause 1.5 \n");
	    else			fprintf(fpm,"  pause 0.05 \n");
	    fprintf(fpm,"%c  load 'saveplot';\n",Movie);
	    fprintf(fpm,"%c  !mv my-plot.ps %s\n", Movie, framefl );
	  }
	  for ( fr = frames-1; fr > 0; fr-- ) {

	    ++frame_number;

	    sprintf(modefl,"%s-%02d.%03d", modepath, m, fr  );
	    sprintf(framefl,"%s-%02d-f-%03d.ps", modepath, m, fr  );

	    fprintf(fpm,"%c plot '%s' u 2:3 w l lw 1 lt 5, ", D2,meshpath );
	    fprintf(fpm," '%s' u 1:2 w l lw %d lt 3; \n", modefl, lw );
	    if ( pan != 0.0 )
	     fprintf(fpm,"%c set view %7.2f, %7.2f, %5.3f # pan = %f\n", D3,
		rot_x_init + pan*(rot_x_final-rot_x_init)*frame_number/total_frames,
		rot_z_init + pan*(rot_z_final-rot_z_init)*frame_number/total_frames,
		zoom_init + pan*(zoom_final-zoom_init)*frame_number/total_frames, pan );
	    fprintf(fpm,"%c splot '%s' u 2:3:4 w l lw 1 lt 5, ",D3,meshpath);
	    fprintf(fpm," '%s' u 1:2:3 w l lw %d lt 3;", modefl, lw );
	    fprintf(fpm,"  pause 0.05 \n");
	    fprintf(fpm,"%c  load 'saveplot';\n",Movie);
	    fprintf(fpm,"%c  !mv my-plot.ps %s\n", Movie, framefl );
	  }

	 }
	 fr = 0;

	 sprintf(modefl,"%s-%02d.%03d", modepath, m, fr  );

	 fprintf(fpm,"%c plot '%s' u 2:3 w l lw %d lt 5, ", D2, meshpath, lw );
	 fprintf(fpm," '%s' u 1:2 w l lw 3 lt 3 \n", modefl );
	 fprintf(fpm,"%c splot '%s' u 2:3:4 w l lw %d lt 5, ",D3,meshpath, lw );
	 fprintf(fpm," '%s' u 1:2:3 w l lw 3 lt 3 \n", modefl );

	 i++;
	}
	fclose(fpm);

	v = dvector(1,DoF);

	i = 1;
	while ( (m = anim[i]) != 0 ) {
	  for ( fr=0; fr<=frames; fr++ ) {

	    sprintf(modefl,"%s-%02d.%03d", modepath, m, fr  );

	    if ((fpm = fopen (modefl, "w")) == NULL) {
		sprintf (errMsg,"\n  error: cannot open gnuplot modal mesh data file: %s \n", modefl);
		errorMsg(errMsg);
		exit(28);
	    }

	    ex = exagg_modal*cos( PI*fr/frames );

	    fprintf(fpm,"# FRAME3DD ANALYSIS RESULTS  http://frame3dd.sf.net/");
	    fprintf(fpm," VERSION %s \n", VERSION);
	    fprintf(fpm,"# %s\n", title );
	    fprintf(fpm,"# A N I M A T E D   M O D E   S H A P E   D A T A \n");
	    fprintf(fpm,"# deflection exaggeration: %.1f\n", ex );
	    fprintf(fpm,"# MODE %5d: f= %lf Hz  T= %lf sec\n\n",m,f[m],1./f[m]);

	    for (j=1; j<=DoF; j++)	v[j] = V[j][m];		/* mode "m" */

	    fprintf(fpm,"#      X-dsp       Y-dsp       Z-dsp\n\n");

	    for (n=1; n<=nE; n++) {
		fprintf( fpm, "\n# element %5d \n", n );
		cubic_bent_beam ( fpm, J1[n], J2[n], xyz, L[n], p[n], v, ex );
	    }

	    fclose(fpm);
	  }
	  i++;
	}

	free_dvector(v,1,DoF);

	return;
}
