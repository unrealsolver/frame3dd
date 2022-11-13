#include "core.h"

#include "NRutil.h"
#include "frame3dd_io.h"

/**
 * Run simulation for given load case
 */
uint8_t solve(
	const InputScope scope, const RuntimeArgs args, Results *results,
	SolverContext ctx, ResultScope rs, const int lc
) {
	// Iterators
	uint16_t i, j;

	// Aliases
	const uint16_t DoF = scope.DoF;
	const uint16_t nL = scope.nL;
	const uint16_t nE = scope.nE;

	uint8_t ExitCode = 0;
	int axial_strain_warning = 0; // 0: "ok", 1: strain > 0.001

	double *dF  = dvector(1,DoF);	/* equilibrium error {F} - [K]{D} */
	double *dD  = dvector(1,DoF);	/* incremental displ. of each node	*/
	double *dR  = dvector(1,DoF);	/* incremental reaction forces		*/

	if ( args.verbose ) {	/* display the load case number  */
		fprintf(stdout,"\n");
		textColor('y','g','b','x');
		fprintf(stdout," Load Case %d of %d ... ", lc, scope.nL );
		fprintf(stdout,"                                          ");

		fflush(stdout);
		color(0);
		fprintf(stdout,"\n");
	}

	/*  initialize displacements and displ. increment to {0}  */
	/*  initialize reactions     and react. increment to {0}  */
	for (i=1; i<=  DoF; i++) rs.D[i] = dD[i] = rs.R[i] = dR[i] = 0.0;

	/*  initialize internal element end forces Q = {0}	*/
	for (i=1; i<= nE; i++) for (j=1;j<=12;j++) rs.Q[i][j] = 0.0;

	/*  elastic stiffness matrix  [K({D}^(i))], {D}^(0)={0} (i=0) */
	assemble_K(
		rs.K, DoF, nE, scope.xyz, scope.rj, scope.L, scope.Le, scope.N1, scope.N2,
		scope.Ax, scope.Asy, scope.Asz, scope.Jx,scope.Iy,scope.Iz, scope.E, scope.G, scope.p,
		scope.shear, scope.geom, rs.Q, args.debug
	);

#ifdef MATRIX_DEBUG
	save_dmatrix ( "Ku", K, 1,DoF, 1,DoF, 0, "w" ); // unloaded stiffness matrix
#endif

	/* first apply temperature loads only, if there are any ... */
	if (scope.nT[lc] > 0) {
		if (args.verbose)
			fprintf(stdout," Linear Elastic Analysis ... Temperature Loads\n");

		/*  solve {F_t} = [K({D=0})] * {D_t} */
		solve_system(rs.K,dD,scope.F_temp[lc],dR,DoF,scope.q,scope.r, &ctx.ok, args.verbose, &ctx.rms_resid);

		/* increment {D_t} = {0} + {D_t} temp.-induced displ */
		for (i=1; i<=DoF; i++)	if (scope.q[i]) rs.D[i] += dD[i];
		/* increment {R_t} = {0} + {R_t} temp.-induced react */
		for (i=1; i<=DoF; i++)	if (scope.r[i]) rs.R[i] += dR[i];

		if (scope.geom) {	/* assemble K = Ke + Kg */
		 /* compute   {Q}={Q_t} ... temp.-induced forces     */
			element_end_forces(
				rs.Q, nE, scope.xyz, scope.L, scope.Le, scope.N1, scope.N2,
				scope.Ax, scope.Asy,scope.Asz, scope.Jx,scope.Iy,scope.Iz, scope.E,scope.G, scope.p,
				scope.eqF_temp[lc], scope.eqF_mech[lc], rs.D, scope.shear, scope.geom,
				&axial_strain_warning
			);
			/* assemble temp.-stressed stiffness [K({D_t})]     */
			assemble_K(
				rs.K, DoF, nE, scope.xyz, scope.rj,
				scope.L, scope.Le, scope.N1, scope.N2,
				scope.Ax,scope.Asy,scope.Asz, scope.Jx,scope.Iy,scope.Iz, scope.E, scope.G, scope.p,
				scope.shear, scope.geom, rs.Q, args.debug
			);
		}
	}

	/* ... then apply mechanical loads only, if there are any ... */
	if (
		scope.nF[lc]>0 || scope.nU[lc]>0 || scope.nW[lc]>0 || scope.nP[lc]>0 || scope.nD[lc]>0 ||
		scope.gX[lc] != 0 || scope.gY[lc] != 0 || scope.gZ[lc] != 0
	) {
		if ( args.verbose )
			fprintf(stdout," Linear Elastic Analysis ... Mechanical Loads\n");
		/* incremental displ at react'ns = prescribed displ */
		for (i=1; i<=DoF; i++)	if (scope.r[i]) dD[i] = scope.Dp[lc][i];

		/*  solve {F_m} = [K({D_t})] * {D_m}	*/
		solve_system(rs.K,dD,scope.F_mech[lc],dR,DoF,scope.q,scope.r, &ctx.ok, args.verbose, &ctx.rms_resid);

		/* combine {D} = {D_t} + {D_m}	*/
		for (i=1; i<=DoF; i++) {
			if (scope.q[i])	rs.D[i] += dD[i];
			else {		rs.D[i]  = scope.Dp[lc][i]; dD[i] = 0.0; }
		}
		/* combine {R} = {R_t} + {R_m} --- for linear systems */
		for (i=1; i<=DoF; i++) if (scope.r[i]) rs.R[i] += dR[i];
	}


	/*  combine {F} = {F_t} + {F_m} */
	for (i=1; i<=DoF; i++)	scope.F[i] = scope.F_temp[lc][i] + scope.F_mech[lc][i];

	/*  element forces {Q} for displacements {D}	*/
	element_end_forces (rs.Q, nE, scope.xyz, scope.L, scope.Le, scope.N1, scope.N2,
			scope.Ax, scope.Asy,scope.Asz, scope.Jx,scope.Iy,scope.Iz, scope.E,scope.G, scope.p,
			scope.eqF_temp[lc], scope.eqF_mech[lc], rs.D, scope.shear, scope.geom,
			&axial_strain_warning );

	/*  check the equilibrium error	*/
	ctx.error = equilibrium_error ( dF, scope.F, rs.K, rs.D, DoF, scope.q,scope.r );

	if ( scope.geom && args.verbose )
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
	if (scope.geom) { ctx.error = 1.0; ctx.ok = 0; ctx.iter = 0; } /* re-initialize */
	while ( scope.geom && ctx.error > scope.tol && ctx.iter < 500 && ctx.ok >= 0) {
		++ctx.iter;

		/*  assemble stiffness matrix [K({D}^(i))]	      */
		assemble_K(rs.K, DoF, nE, scope.xyz, scope.rj, scope.L, scope.Le, scope.N1, scope.N2,
			scope.Ax,scope.Asy,scope.Asz, scope.Jx,scope.Iy,scope.Iz, scope.E, scope.G, scope.p,
			scope.shear,scope.geom, rs.Q, args.debug );

		/*  compute equilibrium error, {dF}, at iteration i   */
		/*  {dF}^(i) = {F} - [K({D}^(i))]*{D}^(i)	      */
		/*  convergence criteria = || {dF}^(i) ||  /  || F || */
		ctx.error = equilibrium_error(dF, scope.F, rs.K, rs.D, DoF, scope.q,scope.r );

		/*  Powell-Symmetric-Broyden secant stiffness update  */
		// PSB_update ( Ks, dF, dD, DoF );  /* not helpful?   */

		/*  solve {dF}^(i) = [K({D}^(i))] * {dD}^(i)	      */
		solve_system(rs.K,dD,dF,dR,DoF,scope.q,scope.r,&ctx.ok,args.verbose, &ctx.rms_resid);

		if ( ctx.ok < 0 ) {	/*  K is not positive definite	      */
			fprintf(stderr,"   The stiffness matrix is not pos-def. \n");
			fprintf(stderr,"   Reduce loads and re-run the analysis.\n");
			ExitCode = 181;
			break;
		}

		/*  increment {D}^(i+1) = {D}^(i) + {dD}^(i)	      */
		for (i=1; i<=DoF; i++)	if (scope.q[i])	rs.D[i] += dD[i];

		/*  element forces {Q} for displacements {D}^(i)      */
		element_end_forces(rs.Q, nE, scope.xyz, scope.L, scope.Le, scope.N1, scope.N2,
			scope.Ax, scope.Asy,scope.Asz, scope.Jx,scope.Iy,scope.Iz, scope.E,scope.G, scope.p,
			scope.eqF_temp[lc], scope.eqF_mech[lc], rs.D, scope.shear, scope.geom,
			&axial_strain_warning );

		if ( args.verbose ) { /*  display equilibrium error        */
		 fprintf(stdout,"   NR iteration %3d ---", ctx.iter);
		 fprintf(stdout," RMS relative equilibrium error = %8.2e \n", ctx.error);
		}
	}			/* end quasi Newton-Raphson iteration */


	/*   strain limit failure ... */
	if (axial_strain_warning > 0 && ExitCode == 0)   ExitCode = 182;
	/*   strain limit _and_ buckling failure ... */
	if (axial_strain_warning > 0 && ExitCode == 181) ExitCode = 183;

	if ( scope.geom ) compute_reaction_forces(rs.R, scope.F, rs.K, rs.D, DoF, scope.r );

	/*  dealocate Broyden secant stiffness matrix, Ks */
	// if ( geom )	free_dmatrix(Ks, 1, DoF, 1, DoF );

	if ( args.write_matrix )	/* write static stiffness matrix */
		save_ut_dmatrix("Ks", rs.K, DoF, "w" );

	/*  display RMS equilibrium error */
	if ( args.verbose && ctx.ok >= 0 ) evaluate ( ctx.error, ctx.rms_resid, scope.tol );

	return ExitCode;
}
