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
*//** @file
	Main FRAME3DD program driver
*//** @mainpage
FRAME3DD: a program for static and dynamic structural analysis of 2D and 3D
frames and trusses with elastic and geometric stiffness.

Also included is a system for parsing Microstran .arc 'Archive' files and
for parsing calculated force and displacement output files (.p1 format) from
Microstran. See @ref mstranp. It is intended that ultimately the .arc format
be an alternative method of inputting data to the FRAME3DD program,
but currently these two parts of the code are distinct.

For more information go to http://frame3dd.sourceforge.net/

The input file format for FRAME is defined in doc/user_manual.html

Henri P. Gavin hpgavin@duke.edu (main FRAME3DD code)
John Pye john.pye@anu.edu.au (Microstran parser and viewer)

For compilation/installation, see README.txt.

*/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "common.h"
#include "compat_types.h"
#include "frame3dd.h"
#include "frame3dd_io.h"
#include "eig.h"
#include "gnuplot_writer.h"
#include "struct_writer.h"
#include "HPGmatrix.h"
#include "HPGutil.h"
#include "NRutil.h"

// compile the Frame3DD analysis into another code, such as a GUI
#ifdef WITH_GLOBALS
 int run_kernel ( int argc, char *argv[] ) {
#endif

// compile Frame3DD to run as a stand-alone code through the terminal
#ifndef WITH_GLOBALS
 int main ( int argc, char *argv[] ) {
#endif

	char	IN_file[FILENMAX],	// the input  data filename
		OUT_file[FILENMAX],	// the output data filename
		title[MAXL],		// the title of the analysis
		errMsg[MAXL],		// the text of an error message
		meshpath[FRAME3DD_PATHMAX] = "EMPTY_MESH", // mesh data path
		plotpath[FRAME3DD_PATHMAX] = "EMPTY_PLOT", // plot file path
		infcpath[FRAME3DD_PATHMAX] = "EMPTY_INFC", // int  file path
		modepath[FRAME3DD_PATHMAX] = "EMPTY_MODE", // mode data path
		strippedInputFile[FRAME3DD_PATHMAX] = "EMPTY_TEMP"; // temp data path

	FILE	*fp;		// input and output file pointer

	float
		scale=1.0,	// zoom scale for 3D plotting in Gnuplot
		dx=1.0;		// x-increment for internal force data

	double	**K=NULL,	// equilibrium stiffness matrix
		// **Ks=NULL,	// Broyden secant stiffness matrix
		traceK = 0.0,	// trace of the global stiffness matrix
		**M = NULL,	// global mass matrix
		traceM = 0.0,	// trace of the global mass matrix
		*R  = NULL,	// total reaction force vector
		*dR = NULL,	// incremental reaction force vector
		*D  = NULL,	// displacement vector
		*dD = NULL,	// incremental displacement vector
		//dDdD = 0.0,	// dD' * dD
		*dF = NULL,	// equilibrium error in nonlinear anlys
		*f  = NULL,	// resonant frequencies
		**V = NULL,	// resonant mode-shapes
		rms_resid=1.0,	// root mean square of residual displ. error
		error = 1.0,	// rms equilibrium error and reactions
		Cfreq = 0.0,	// frequency used for Guyan condensation
		**Kc, **Mc,	// condensed stiffness and mass matrices
		exagg_static=10;// exaggerate static displ. in mesh data

		// peak internal forces, moments, and displacments
		// in each frame element and each load case
	double	**pkNx, **pkVy, **pkVz, **pkTx, **pkMy, **pkMz,
		**pkDx, **pkDy, **pkDz, **pkRx, **pkSy, **pkSz;

	int	nN=0,		// number of Nodes
		nE=0,		// number of frame Elements
		nL=0, lc=0,	// number of Load cases
		DoF=0, i, j,	// number of Degrees of Freedom
		shear=0,	// indicates shear deformation
		geom=0,		// indicates  geometric nonlinearity
		anlyz=1,	// 1: stiffness analysis, 0: data check
		nM_calc,	// number of modes to calculate
		iter=0,		// number of iterations
		ok=1,		// number of (-ve) diag. terms of L D L'
		filetype=0,	// 1 if .CSV, 2 if file is Matlab
		axial_strain_warning = 0, // 0: "ok", 1: strain > 0.001
		ExitCode = 0;	// error code returned by Frame3DD

	int	sfrv=0;		// *scanf return value for err checking
	char	extension[16];	// Input Data file name extension

	InputScope scope; // All input data combined; Compatibility adaptor

	Frame *frame = (Frame *) malloc(sizeof(Frame));
	LoadCases *load_cases = (LoadCases *) malloc(sizeof(LoadCases));
	Results *results = (Results *) malloc(sizeof(Results));
	RunOptions *run_options = (RunOptions *) malloc(sizeof(RunOptions));
	run_options->simulation.dx = dx;
	run_options->simulation.nonlinear = geom;
	run_options->simulation.shear = shear;
	run_options->visual.pan = 1.0;
	run_options->visual.scale = scale;
	run_options->visual.exagg_static = exagg_static;
	run_options->visual.exagg_modal = exagg_static;

	const RuntimeArgs args = parse_options ( argc, argv, IN_file, OUT_file );
	// FIXME tmp aliases
	const int8_t debug = args.debug;
	const int8_t verbose = args.verbose;

	if ( verbose ) { /*  display program name, version and license type */
		textColor('w','b','b','x');
		fprintf(stdout,"\n FRAME3DD version: %s\n", VERSION);
		fprintf(stdout," Analysis of 2D and 3D structural frames with elastic and geometric stiffness.\n");
		fprintf(stdout," http://frame3dd.sf.net\n");
		fprintf(stdout," GPL Copyright (C) 1992-2014, Henri P. Gavin\n");
		fprintf(stdout," This is free software with absolutely no warranty.\n");
		fprintf(stdout," For details, see the GPL license file, LICENSE.txt\n");
		color(0); fprintf(stdout,"\n");
	}

	/* open the input data file */

	if ((fp = fopen (IN_file, "r")) == NULL) { /* open input data file */
		sprintf (errMsg,"\n ERROR: cannot open input data file '%s'\n", IN_file);
		errorMsg(errMsg);
		display_help();
		if ( argc == 1 ) {
			fprintf(stderr," Press the 'Enter' key to close.\n");
			(void) getchar();	// clear the buffer ??
			while( !getchar() ) ;	// wait for the Enter key
		}
		exit(11);
	}

	filetype = get_file_ext( IN_file, extension ); /* .CSV or .FMM or other? */

//	temp_file_location("frame3dd.3dd",strippedInputFile,FRAME3DD_PATHMAX);
	output_path("frame3dd.3dd", strippedInputFile, FRAME3DD_PATHMAX, NULL);

	parse_input(fp, strippedInputFile);	/* strip comments from input data */
	fclose(fp);

	if ((fp = fopen (strippedInputFile, "r")) == NULL) { /* open stripped input file */
		sprintf(errMsg,"\n ERROR: cannot open stripped input data file '%s'\n", strippedInputFile);
		errorMsg(errMsg);
		exit(13);
	}

	frame3dd_getline(fp, title, MAXL);
	if ( verbose ) {	/*  display analysis title */
		textColor('w','g','b','x');
		fprintf(stdout,"\n");
		fprintf(stdout," ** %s ** \n", title );
		color(0);
		fprintf(stdout,"\n");
	}

	/* Start actually parsing file */
	// Read nN
	read_node_number(fp, &nN, verbose);
	scope.nN = nN;
	frame->nodes.size = nN;
	// FIXME deallocte
	frame->nodes.data = (Node *) calloc(nN, sizeof(Node));

	read_node_data(fp, frame, &scope);
	if ( verbose )	printf(" ... complete\n");

	DoF = 6*nN;		/* total number of degrees of freedom	*/
	scope.DoF = 6 * nN;

	// TODO read to Frame
	read_reaction_data (fp, verbose, frame, &scope);
	if ( verbose )	fprintf(stdout," ... complete\n");

	read_element_number(fp, &nE, verbose);
	scope.nE = nE;

	if ( nN > nE + 1) {	/* not enough elements */
		fprintf(stderr,"\n  warning: %d nodes and %d members...", nN, nE );
		fprintf(stderr," not enough elements to connect all nodes.\n");
	}

	frame->edges.size = nE;
	// FIXME deallocte
	frame->edges.data = (Edge *) calloc(nE, sizeof(Edge));

	read_frame_element_data(fp, frame, &scope);
	if ( verbose) 	fprintf(stdout," ... complete\n");

	read_run_data ( fp, OUT_file, &shear, &geom,
			meshpath, plotpath, infcpath,
			&exagg_static, &scale, &dx,
			&anlyz, run_options, args);

	sfrv=fscanf(fp, "%d", &nL );	/* number of load cases		*/
	load_cases->size = nL;
	results->size = nL;
	scope.nL = nL;

	if (sfrv != 1)	sferr("nL value for number of load cases");
	if ( verbose ) {	/* display nL */
		fprintf(stdout," number of load cases ");
		dots(stdout,31); fprintf(stdout," nL = %3d \n", load_cases->size);
	}

	if (load_cases->size < 1) {	/* not enough load cases */
		errorMsg("\n ERROR: the number of load cases must be at least 1\n");
		exit(101);
	}
	if (load_cases->size >= _NL_) { /* too many load cases */
		sprintf(errMsg,"\n ERROR: maximum of %d load cases allowed\n", _NL_-1);
		errorMsg(errMsg);
		exit(102);
	}

	// FIXME Deallocate
	load_cases->data = (LoadCase *) calloc(nL, sizeof(LoadCase));
	// FIXME deallocte
	results->data = (LoadCaseResult*) calloc(nL, sizeof(LoadCaseResult));

	for (unsigned i = 0; i < nL; i++) {
		load_cases->data[i].loads.point.size = 0;
		load_cases->data[i].loads.uniform.size = 0;
		// FIXME deallocte
		results->data[i].displacements = calloc(nN, sizeof(NodeDisplacement));
		// FIXME deallocte
		results->data[i].edges = calloc(nE, sizeof(EdgeResult));
	}


	dF	= dvector(1,DoF);	/* equilibrium error {F} - [K]{D} */

	K   = dmatrix(1,DoF,1,DoF);	/* global stiffness matrix	*/

	D   = dvector(1,DoF);	/* displacments of each node		*/
	dD  = dvector(1,DoF);	/* incremental displ. of each node	*/
	R   = dvector(1,DoF);	/* reaction forces			*/
	dR  = dvector(1,DoF);	/* incremental reaction forces		*/


	// peak axial forces, shears, torques, and moments along each element
	pkNx = dmatrix(1,nL,1,nE);
	pkVy = dmatrix(1,nL,1,nE);
	pkVz = dmatrix(1,nL,1,nE);
	pkTx = dmatrix(1,nL,1,nE);
	pkMy = dmatrix(1,nL,1,nE);
	pkMz = dmatrix(1,nL,1,nE);

	// peak displacements and slopes along each element
	pkDx = dmatrix(1,nL,1,nE);
	pkDy = dmatrix(1,nL,1,nE);
	pkDz = dmatrix(1,nL,1,nE);
	pkRx = dmatrix(1,nL,1,nE);
	pkSy = dmatrix(1,nL,1,nE);
	pkSz = dmatrix(1,nL,1,nE);

	read_and_assemble_loads(fp, shear, verbose, frame, load_cases, &scope);

	if ( verbose ) {	/* display load data complete */
		fprintf(stdout,"                                                     ");
		fprintf(stdout," load data ... complete\n");
	}

	read_mass_data(fp, IN_file, modepath, args, &scope);
	if ( verbose ) {	/* display mass data complete */
		fprintf(stdout,"                                                     ");
		fprintf(stdout," mass data ... complete\n");
	}

	read_condensation_data(fp, verbose, args, &scope);

	if( scope.nC>0 && verbose ) {	/*  display condensation data complete */
		fprintf(stdout,"                                      ");
		fprintf(stdout," matrix condensation data ... complete\n");
	}

	fclose(fp);		/* close the input data file */

	fp = fopen(OUT_file, "a"); /* open the output data file for appending */

	if(fp==NULL) {	/* unable to append to output data file */
		fprintf(stderr,"Unable to append to output data file '%s'!\n",
								OUT_file);
		exit(14);
	}

	write_input_data(
		fp, title, scope.nD, scope.nR, scope.nF, scope.nU, scope.nW, scope.nP, scope.nT,
		scope.p, scope.d,
		scope.F_temp, scope.F_mech, scope.Dp, scope.r, scope.U, scope.W, scope.P, scope.T,
		shear, anlyz, geom, frame, load_cases
	);


	if ( anlyz ) {			/* solve the problem	*/
		srand(time(NULL));
		for (lc=1; lc<=nL; lc++) {	/* begin load case analysis loop */
			if ( verbose ) {	/* display the load case number  */
				fprintf(stdout,"\n");
				textColor('y','g','b','x');
				fprintf(stdout," Load Case %d of %d ... ", lc,nL );
				fprintf(stdout,"                                          ");

				fflush(stdout);
				color(0);
				fprintf(stdout,"\n");
			}

			LoadCaseResult *lc_result = &results->data[lc - 1];

			/*  initialize displacements and displ. increment to {0}  */
			/*  initialize reactions     and react. increment to {0}  */
			for (i=1; i<=DoF; i++) D[i] = dD[i] = R[i] = dR[i] = 0.0;

			// FIXME Possible duplicates initialization in the read_and_assemble_loads
			/*  initialize internal element end forces Q = {0}	*/
			for (i=1; i<=nE; i++) for (j=1;j<=12;j++) scope.Q[i][j] = 0.0;

			/*  elastic stiffness matrix  [K({D}^(i))], {D}^(0)={0} (i=0) */
			assemble_K(
				K, DoF, nE, scope.xyz, scope.rj, scope.L, scope.Le, scope.N1, scope.N2,
				scope.Ax, scope.Asy, scope.Asz, scope.Jx,scope.Iy,scope.Iz, scope.E, scope.G, scope.p,
				shear, geom, scope.Q, debug
			);

#ifdef MATRIX_DEBUG
			save_dmatrix ( "Ku", K, 1,DoF, 1,DoF, 0, "w" ); // unloaded stiffness matrix
#endif

			/* first apply temperature loads only, if there are any ... */
			if (scope.nT[lc] > 0) {
				if (verbose)
					fprintf(stdout," Linear Elastic Analysis ... Temperature Loads\n");

				/*  solve {F_t} = [K({D=0})] * {D_t} */
				solve_system(K,dD,scope.F_temp[lc],dR,DoF,scope.q,scope.r,&ok,verbose,&rms_resid);

				/* increment {D_t} = {0} + {D_t} temp.-induced displ */
				for (i=1; i<=DoF; i++)	if (scope.q[i]) D[i] += dD[i];
				/* increment {R_t} = {0} + {R_t} temp.-induced react */
				for (i=1; i<=DoF; i++)	if (scope.r[i]) R[i] += dR[i];

				if (geom) {	/* assemble K = Ke + Kg */
				 /* compute   {Q}={Q_t} ... temp.-induced forces     */
					element_end_forces(
						scope.Q, nE, scope.xyz, scope.L, scope.Le, scope.N1, scope.N2,
						scope.Ax, scope.Asy,scope.Asz, scope.Jx,scope.Iy,scope.Iz, scope.E,scope.G, scope.p,
						scope.eqF_temp[lc], scope.eqF_mech[lc], D, shear, geom,
						&axial_strain_warning
					);
					/* assemble temp.-stressed stiffness [K({D_t})]     */
					assemble_K(
						K, DoF, nE, scope.xyz, scope.rj,
						scope.L, scope.Le, scope.N1, scope.N2,
						scope.Ax,scope.Asy,scope.Asz, scope.Jx,scope.Iy,scope.Iz, scope.E, scope.G, scope.p,
						shear,geom, scope.Q, debug
					);
				}
			}

			/* ... then apply mechanical loads only, if there are any ... */
			if (
				scope.nF[lc]>0 || scope.nU[lc]>0 || scope.nW[lc]>0 || scope.nP[lc]>0 || scope.nD[lc]>0 ||
				scope.gX[lc] != 0 || scope.gY[lc] != 0 || scope.gZ[lc] != 0
			) {
				if ( verbose )
					fprintf(stdout," Linear Elastic Analysis ... Mechanical Loads\n");
				/* incremental displ at react'ns = prescribed displ */
				for (i=1; i<=DoF; i++)	if (scope.r[i]) dD[i] = scope.Dp[lc][i];

				/*  solve {F_m} = [K({D_t})] * {D_m}	*/
				solve_system(K,dD,scope.F_mech[lc],dR,DoF,scope.q,scope.r,&ok,verbose,&rms_resid);

				/* combine {D} = {D_t} + {D_m}	*/
				for (i=1; i<=DoF; i++) {
					if (scope.q[i])	D[i] += dD[i];
					else {		D[i]  = scope.Dp[lc][i]; dD[i] = 0.0; }
				}
				/* combine {R} = {R_t} + {R_m} --- for linear systems */
				for (i=1; i<=DoF; i++) if (scope.r[i]) R[i] += dR[i];
			}


			/*  combine {F} = {F_t} + {F_m} */
			for (i=1; i<=DoF; i++)	scope.F[i] = scope.F_temp[lc][i] + scope.F_mech[lc][i];

			/*  element forces {Q} for displacements {D}	*/
			element_end_forces ( scope.Q, nE, scope.xyz, scope.L, scope.Le, scope.N1, scope.N2,
					scope.Ax, scope.Asy,scope.Asz, scope.Jx,scope.Iy,scope.Iz, scope.E,scope.G, scope.p,
					scope.eqF_temp[lc], scope.eqF_mech[lc], D, shear, geom,
					&axial_strain_warning );

			/*  check the equilibrium error	*/
			error = equilibrium_error ( dF, scope.F, K, D, DoF, scope.q,scope.r );

			if ( geom && verbose )
				fprintf(stdout,"\n Non-Linear Elastic Analysis ...\n");

/*
 *	 		if ( geom ) { // initialize Broyden secant stiffness matrix, Ks
 *				Ks  = dmatrix( 1, DoF, 1, DoF );
 *				for (i=1;i<=DoF;i++) {
 *					for(j=i;j<=DoF;j++) {
 *						Ks[i][j]=Ks[j][i]=K[i][j];
 *					}
 *				}
 *			}
 */

			/* quasi Newton-Raphson iteration for geometric nonlinearity  */
			if (geom) { error = 1.0; ok = 0; iter = 0; } /* re-initialize */
			while ( geom && error > scope.tol && iter < 500 && ok >= 0) {
				++iter;

				/*  assemble stiffness matrix [K({D}^(i))]	      */
				assemble_K ( K, DoF, nE, scope.xyz, scope.rj, scope.L, scope.Le, scope.N1, scope.N2,
					scope.Ax,scope.Asy,scope.Asz, scope.Jx,scope.Iy,scope.Iz, scope.E, scope.G, scope.p,
					shear,geom, scope.Q, debug );

				/*  compute equilibrium error, {dF}, at iteration i   */
				/*  {dF}^(i) = {F} - [K({D}^(i))]*{D}^(i)	      */
				/*  convergence criteria = || {dF}^(i) ||  /  || F || */
				error = equilibrium_error ( dF, scope.F, K, D, DoF, scope.q,scope.r );

				/*  Powell-Symmetric-Broyden secant stiffness update  */
				// PSB_update ( Ks, dF, dD, DoF );  /* not helpful?   */

				/*  solve {dF}^(i) = [K({D}^(i))] * {dD}^(i)	      */
				solve_system(K,dD,dF,dR,DoF,scope.q,scope.r,&ok,verbose,&rms_resid);

				if ( ok < 0 ) {	/*  K is not positive definite	      */
					fprintf(stderr,"   The stiffness matrix is not pos-def. \n");
					fprintf(stderr,"   Reduce loads and re-run the analysis.\n");
					ExitCode = 181;
					break;
				}

				/*  increment {D}^(i+1) = {D}^(i) + {dD}^(i)	      */
				for (i=1; i<=DoF; i++)	if (scope.q[i])	D[i] += dD[i];

				/*  element forces {Q} for displacements {D}^(i)      */
				element_end_forces ( scope.Q, nE, scope.xyz, scope.L, scope.Le, scope.N1, scope.N2,
					scope.Ax, scope.Asy,scope.Asz, scope.Jx,scope.Iy,scope.Iz, scope.E,scope.G, scope.p,
					scope.eqF_temp[lc], scope.eqF_mech[lc], D, shear, geom,
					&axial_strain_warning );

				if ( verbose ) { /*  display equilibrium error        */
				 fprintf(stdout,"   NR iteration %3d ---", iter);
				 fprintf(stdout," RMS relative equilibrium error = %8.2e \n",error);
				}
			}			/* end quasi Newton-Raphson iteration */


			/*   strain limit failure ... */
			if (axial_strain_warning > 0 && ExitCode == 0)   ExitCode = 182;
			/*   strain limit _and_ buckling failure ... */
			if (axial_strain_warning > 0 && ExitCode == 181) ExitCode = 183;

			if ( geom )	compute_reaction_forces( R,scope.F,K, D, DoF, scope.r );

			/*  dealocate Broyden secant stiffness matrix, Ks */
			// if ( geom )	free_dmatrix(Ks, 1, DoF, 1, DoF );

			if ( args.write_matrix )	/* write static stiffness matrix */
				save_ut_dmatrix ( "Ks", K, DoF, "w" );

			/*  display RMS equilibrium error */
			if ( verbose && ok >= 0 ) evaluate ( error, rms_resid, scope.tol );

			write_static_struct(lc_result, frame, D, scope.Q);

			write_static_results ( fp, lc, DoF, scope.N1, scope.N2,
					scope.F,D,R, scope.r,scope.Q, rms_resid, ok, args.axial_sign, frame, load_cases );

			if ( filetype == 1 ) {		// .CSV format output
				write_static_csv(OUT_file, title,
				    lc, DoF, scope.N1, scope.N2, scope.F,D,R, scope.r,scope.Q, error, ok, frame, load_cases );
			}

			if ( filetype == 2 ) {		// .m matlab format output
				write_static_mfile (OUT_file, title, lc, DoF,
						scope.N1,scope.N2, scope.F,D,R, scope.r,scope.Q, error, ok, frame, load_cases);
			}

/*
 * 			if ( verbose )
 * 			 printf("\n   If the program pauses here for very long,"
 * 			 " hit CTRL-C to stop execution, \n"
 * 			 "    reduce exagg_static in the Input Data,"
 * 			 " and re-run the analysis. \n");
 */

			write_internal_forces ( OUT_file, fp, infcpath, lc, title, dx, scope.xyz,
						scope.Q, scope.L, scope.N1, scope.N2,
						scope.p,
						scope.d, scope.gX[lc], scope.gY[lc], scope.gZ[lc],
						scope.nU[lc],scope.U[lc],scope.nW[lc],scope.W[lc],scope.nP[lc],scope.P[lc],
						D, shear, error, frame, load_cases );

			static_mesh ( IN_file, infcpath, meshpath, plotpath, title,
						lc, DoF,
						scope.xyz, scope.L, scope.N1, scope.N2, scope.p, D,
						exagg_static, anlyz,
						dx, scale, load_cases, frame, args );

		} /* end load case loop */
	} else {		/*  data check only  */
		if ( verbose ) {	/* display data check only */
			fprintf(stdout,"\n * %s *\n", title );
			fprintf(stdout,"  DATA CHECK ONLY.\n");
		}

		static_mesh(
			IN_file, infcpath, meshpath, plotpath, title,
			lc, DoF,
			scope.xyz, scope.L, scope.N1, scope.N2, scope.p, D,
			exagg_static, anlyz,
			dx, scale, load_cases, frame, args
		);
	}


	if ( scope.nM > 0 ) { /* carry out modal analysis */

		if(verbose & anlyz) fprintf(stdout,"\n\n Modal Analysis ...\n");

		nM_calc = scope.nM + 8 < scope.nM * 2
			? scope.nM + 8
			: scope.nM * 2;		/* Bathe */

		M   = dmatrix(1,DoF,1,DoF);
		f   = dvector(1,nM_calc);
		V   = dmatrix(1,DoF,1,nM_calc);

		assemble_M ( M, DoF, nN, nE, scope.xyz, scope.rj, scope.L, scope.N1, scope.N2,
				scope.Ax, scope.Jx,scope.Iy,scope.Iz, scope.p, scope.d, scope.EMs, scope.NMs, scope.NMx, scope.NMy, scope.NMz,
				scope.lump, debug );

#ifdef MATRIX_DEBUG
		save_dmatrix ( "Mf", M, 1,DoF, 1,DoF, 0, "w" ); /* free mass matrix */
#endif

		for (j=1; j<=DoF; j++) { /*  compute traceK and traceM */
			if ( !scope.r[j] ) {
				traceK += K[j][j];
				traceM += M[j][j];
			}
		}
		for (i=1; i<=DoF; i++) { /*  modify K and M for reactions    */
			if ( scope.r[i] ) {	/* apply reactions to upper triangle */
				K[i][i] = traceK * 1e4;
				M[i][i] = traceM;
				for (j=i+1; j<=DoF; j++)
					K[j][i]=K[i][j]=M[j][i]=M[i][j] = 0.0;
		    }
		}

		if ( args.write_matrix ) {	/* write Kd and Md matrices */
			save_ut_dmatrix ( "Kd", K, DoF, "w" );/* dynamic stff matx */
			save_ut_dmatrix ( "Md", M, DoF, "w" );/* dynamic mass matx */
		}

		if ( anlyz ) {	/* subspace or stodola methods */
			if( scope.Mmethod == 1 )
				subspace( K, M, DoF, nM_calc, f, V, scope.tol,scope.shift,&iter,&ok, verbose );
			if( scope.Mmethod == 2 )
				stodola ( K, M, DoF, nM_calc, f, V, scope.tol,scope.shift,&iter,&ok, verbose );

			for (j=1; j<=nM_calc; j++) f[j] = sqrt(f[j])/(2.0*PI);

			write_modal_results ( fp, scope.nI, DoF, M,f,V,
					scope.total_mass, scope.struct_mass,
					iter, scope.sumR, scope.nM, scope.shift, scope.lump, scope.tol, ok, frame );
		}
	}

	fprintf(fp,"\n");
	fclose (fp);

	if ( scope.nM > 0 && anlyz ) {	/* write modal analysis results */
		modal_mesh ( IN_file, meshpath, modepath, plotpath, title,
				DoF, scope.nM, scope.xyz, scope.L, scope.N1, scope.N2, scope.p,
				M, f, V, scope.exagg_modal, anlyz, frame, args );

		animate ( IN_file, meshpath, modepath, plotpath, title,scope.anim,
				DoF, scope.nM, scope.xyz, scope.L, scope.p, scope.N1, scope.N2, f,
				V, scope.exagg_modal, scope.pan, scale, frame, args );
	}

	if ( scope.nC > 0 ) {		/* matrix condensation of stiffness and mass */

		if ( verbose ) fprintf(stdout,"\n Matrix Condensation ...\n");

		if(scope.Cdof > scope.nM && scope.Cmethod == 3){
			fprintf(stderr,"  Cdof > nM ... Cdof = %d  nM = %d \n",
				 scope.Cdof, scope.nM );
			fprintf(stderr,"  The number of condensed degrees of freedom");
			fprintf(stderr," may not exceed the number of computed modes");
			fprintf(stderr," when using dynamic condensation.\n");
			exit(94);
		}

		Kc = dmatrix(1,scope.Cdof,1,scope.Cdof);
		Mc = dmatrix(1,scope.Cdof,1,scope.Cdof);

		if ( scope.m[1] > 0 && scope.nM > 0 )	Cfreq = f[scope.m[1]];

		if ( scope.Cmethod == 1 && anlyz) {	/* static condensation only */
			static_condensation(K, DoF, scope.c, scope.Cdof, Kc, 0 );
			if ( verbose )
				fprintf(stdout,"   static condensation of K complete\n");
		}
		if ( scope.Cmethod == 2 && anlyz ) {  /*  dynamic condensation  */
			paz_condensation(M, K, DoF, scope.c, scope.Cdof, Mc,Kc, Cfreq, 0 );
			if ( verbose ) {
				fprintf(stdout,"   Paz condensation of K and M complete");
				fprintf(stdout," ... dynamics matched at %f Hz.\n", Cfreq );
			}
		}
		if ( scope.Cmethod == 3 && scope.nM > 0 && anlyz ) {
			modal_condensation(M,K, DoF, scope.r, scope.c, scope.Cdof, Mc,Kc, V,f, scope.m, 0 );
			if ( verbose )
				fprintf(stdout,"   modal condensation of K and M complete\n");
		}
		save_dmatrix("Kc", Kc, 1,scope.Cdof, 1,scope.Cdof, 0, "w" );
		save_dmatrix("Mc", Mc, 1,scope.Cdof, 1,scope.Cdof, 0, "w" );

		free_dmatrix(Kc, 1,scope.Cdof,1,scope.Cdof );
		free_dmatrix(Mc, 1,scope.Cdof,1,scope.Cdof );
	}

	/* deallocate memory used for each frame analysis variable */
	deallocate ( nN, nE, nL, scope.nF, scope.nU, scope.nW, scope.nP, scope.nT, DoF, scope.nM,
			scope.xyz, scope.rj, scope.L, scope.Le, scope.N1, scope.N2, scope.q,scope.r,
			scope.Ax, scope.Asy, scope.Asz, scope.Jx, scope.Iy, scope.Iz, scope.E, scope.G, scope.p,
			scope.U,scope.W,scope.P,scope.T, scope.Dp, scope.F_mech, scope.F_temp,
			scope.eqF_mech, scope.eqF_temp, scope.F, dF,
			K, scope.Q, D, dD, R, dR,
			scope.d,scope.EMs,scope.NMs,scope.NMx,scope.NMy,scope.NMz, M,f,V, scope.c, scope.m,
			pkNx, pkVy, pkVz, pkTx, pkMy, pkMz,
			pkDx, pkDy, pkDz, pkRx, pkSy, pkSz
	);

	free(frame);
	free(load_cases);
	free(results);

	if ( verbose ) fprintf(stdout,"\n");

	if ( argc == 1 ) { /* wait for keyboard entry to close the terminal */
	   fprintf(stderr," The Output Data was appended to %s \n", OUT_file );
	   fprintf(stderr," A Gnuplot script was written to %s \n", plotpath );
	   fprintf(stderr," Press the 'Enter' key to close.\n");
	   (void) getchar();	// clear the buffer ??
	   while( !getchar() ) ;	// wait for the Enter key to be hit
	}
	color(0);

	return( ExitCode );
}
