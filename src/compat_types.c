#include "compat_types.h"
#include "NRutil.h"
#include "HPGutil.h"
#include "coordtrans.h"
#include <string.h>

Error *Error_new(const uint8_t code, const char *message) {
	Error *err = calloc(1, sizeof(Error));
	err->code = code;
	err->message = message;
	return err;
};

void Error_destroy(Error *self) {
	free(self);
};

void Error_handle(Error *self) {
	if (self != NULL) {
		const char *msg = self->message;
		const uint8_t code = self->code;
		Error_destroy(self);
		errorMsg(msg);
		exit(code);
	}
}


void RS_init_for_IS(ResultScope *self, const InputScope *is) {
	const uint16_t DoF = is->DoF;
	const uint16_t nE = is->nE;
	self->K = dmatrix(1,DoF,1,DoF);	/* global stiffness matrix	*/
	self->Q = dmatrix(1,nE,1,12);	/* end forces for each member	*/
	self->D = dvector(1,DoF);	/* displacments of each node	*/
	self->R = dvector(1,DoF);	/* reaction forces		*/
};

void IS_set_nN(InputScope *self, const uint16_t nN) {
	self->nN = nN;
	self->DoF = nN * 6;
	self->rj = vector(1, nN);	/* rigid radius around each node */
	self->xyz = (vec3 *) calloc(nN + 1, sizeof(vec3)); /* node coordinates */
}

void IS_set_nE(InputScope *self, const uint16_t nE) {
	self->nE  = nE;
	self->L   = dvector(1,nE);	/* length of each element		*/
	self->Le  = dvector(1,nE);	/* effective length of each element	*/

	self->N1  = ivector(1,nE);	/* node #1 of each element		*/
	self->N2  = ivector(1,nE);	/* node #2 of each element		*/

	self->Ax  =  vector(1,nE);	/* cross section area of each element	*/
	self->Asy =  vector(1,nE);	/* shear area in local y direction 	*/
	self->Asz =  vector(1,nE);	/* shear area in local z direction	*/
	self->Jx  =  vector(1,nE);	/* torsional moment of inertia 		*/
	self->Iy  =  vector(1,nE);	/* bending moment of inertia about y-axis */
	self->Iz  =  vector(1,nE);	/* bending moment of inertia about z-axis */

	self->E   =  vector(1,nE);	/* frame element Young's modulus	*/
	self->G   =  vector(1,nE);	/* frame element shear modulus		*/
	self->p   =  vector(1,nE);	/* element rotation angle about local x axis */
	self->d   =  vector(1,nE);	/* element mass density			*/
}

void IS_set_nL(InputScope *self, const uint8_t nL) {
	const uint16_t nN  = self->nN;
	const uint16_t nE  = self->nE;
	const uint16_t DoF = self->DoF;

	self->nL  = nL;
	self->U   =  D3matrix(1,nL,1,nE,1,4);    /* uniform load on each member */
	self->W   =  D3matrix(1,nL,1,10*nE,1,13);/* trapezoidal load on each member */
	self->P   =  D3matrix(1,nL,1,10*nE,1,5); /* internal point load each member */
	self->T   =  D3matrix(1,nL,1,nE,1,8);    /* internal temp change each member*/
	self->Dp  =  matrix(1,nL,1,DoF); /* prescribed displacement of each node */

	self->F_mech  = dmatrix(1,nL,1,DoF);	/* mechanical load vector	*/
	self->F_temp  = dmatrix(1,nL,1,DoF);	/* temperature load vector	*/
	self->F       = dvector(1,DoF);	/* external load vector	*/

	self->eqF_mech =  D3dmatrix(1,nL,1,nE,1,12); /* eqF due to mech loads */
	self->eqF_temp =  D3dmatrix(1,nL,1,nE,1,12); /* eqF due to temp loads */

	/* initialize load data vectors and matrices to zero */
	for (uint16_t j=1; j<=DoF; j++)	self->F[j] = 0.0;
	for (uint16_t j=1; j<=DoF; j++)
		for (uint16_t lc=1; lc <= nL; lc++)
			self->F_temp[lc][j] = self->F_mech[lc][j] = 0.0;
	for (uint16_t i=1; i<=12; i++)
		for (uint16_t n=1; n<=nE; n++)
			for (uint16_t lc=1; lc <= nL; lc++)
				self->eqF_mech[lc][n][i] = self->eqF_temp[lc][n][i] = 0.0;

	for (uint16_t i=1; i<=DoF; i++)
		for (uint16_t lc=1; lc<=nL; lc++)
			self->Dp[lc][i] = 0.0;

}

void IS_init_eqF_mech(InputScope *self, const uint8_t lc) {
	const float gX = self->gX[lc];
	const float gY = self->gY[lc];
	const float gZ = self->gZ[lc];
	const double *L = self->L;
	const float *p = self->p;
	const float *d = self->d;

	double t1, t2, t3, t4, t5, t6, t7, t8, t9;	/* 3D coord Xfrm coeffs */

	for (uint16_t n = 1; n <= self->nE; n++) {
		const float Ax = self->Ax[n];
		const uint16_t n1 = self->N1[n];
		const uint16_t n2 = self->N2[n];

		coord_trans(self->xyz, L[n], n1, n2,
			&t1, &t2, &t3, &t4, &t5, &t6, &t7, &t8, &t9, p[n]);

		// CONSIDER REDEFINITIONS
		self->eqF_mech[lc][n][1]  = d[n] * Ax * L[n] * gX / 2.0;
		self->eqF_mech[lc][n][2]  = d[n] * Ax * L[n] * gY / 2.0;
		self->eqF_mech[lc][n][3]  = d[n] * Ax * L[n] * gZ / 2.0;

		self->eqF_mech[lc][n][4]  = d[n] * Ax * L[n]*L[n] / 12.0 *
			( (-t4*t8+t5*t7)*gY + (-t4*t9+t6*t7)*gZ );
		self->eqF_mech[lc][n][5]  = d[n] * Ax *L[n]*L[n] / 12.0 *
			( (-t5*t7+t4*t8)*gX + (-t5*t9+t6*t8)*gZ );
		self->eqF_mech[lc][n][6]  = d[n] * Ax *L[n]*L[n] / 12.0 *
		   ( (-t6*t7+t4*t9)*gX + (-t6*t8+t5*t9)*gY );

		self->eqF_mech[lc][n][7]  = d[n] * Ax * L[n] * gX / 2.0;
		self->eqF_mech[lc][n][8]  = d[n] * Ax * L[n] * gY / 2.0;
		self->eqF_mech[lc][n][9]  = d[n] * Ax * L[n] * gZ / 2.0;

		self->eqF_mech[lc][n][10] = d[n] * Ax * L[n]*L[n] / 12.0 *
		   ( ( t4*t8-t5*t7)*gY + ( t4*t9-t6*t7)*gZ );
		self->eqF_mech[lc][n][11] = d[n] * Ax * L[n]*L[n] / 12.0 *
		   ( ( t5*t7-t4*t8)*gX + ( t5*t9-t6*t8)*gZ );
		self->eqF_mech[lc][n][12] = d[n] * Ax * L[n]*L[n] / 12.0 *
			( ( t6*t7-t4*t9)*gX + ( t6*t8-t5*t9)*gY );

		/* debugging ... check eqF data
		printf("n=%d ", n);
		for (l=1;l<=12;l++) {
			if (eqF_mech[lc][n][l] != 0)
			   printf(" eqF %d = %9.2e ", l, eqF_mech[lc][n][l] );
		}
		printf("\n");
		*/
	}					/* end gravity loads */
}

/**
 * assemble all element equivalent loads into
 * separate load vectors for mechanical and thermal loading
 */
void IS_assemble_eq_loads(InputScope *self, uint8_t lc) {
	const uint16_t nE  = self->nE;
	for (uint16_t n = 1; n <= nE; n++) {
		uint16_t n1 = self->N1[n];
		uint16_t n2 = self->N2[n];
		uint8_t i;
		for (i = 1; i <= 6;  i++) self->F_mech[lc][6*n1- 6+i] += self->eqF_mech[lc][n][i];
		for (i = 7; i <= 12; i++) self->F_mech[lc][6*n2-12+i] += self->eqF_mech[lc][n][i];
		for (i = 1; i <= 6;  i++) self->F_temp[lc][6*n1- 6+i] += self->eqF_temp[lc][n][i];
		for (i = 7; i <= 12; i++) self->F_temp[lc][6*n2-12+i] += self->eqF_temp[lc][n][i];
	}
}

Error *IS_init_reactions(InputScope *self) {
	char errMsg[MAXL];
	self->q = ivector(1, self->DoF);	/* allocate memory for reaction data ... */

	// Calculate total number of restrained nodes
	self->sumR = 0;
	for (uint8_t i = 1; i <= self->DoF; i++)
		self->sumR += self->r[i];

	if (self->sumR < 4) {
		sprintf(errMsg,
			"\n  Warning:  un-restrained structure   %d imposed reactions.\n"
			"  At least 4 reactions are required to support static loads.\n",
			self->sumR
		);
		return Error_new(84, errMsg);
	}

	if (self->sumR >= self->DoF) {
		sprintf(errMsg,
			"\n  error in reaction data:  Fully restrained structure\n"
			"   %d imposed reactions >= %d degrees of freedom\n",
			self->sumR, self->DoF
		);
		return Error_new(85, errMsg);
	}

	// Create inverted matrix for q matrix
	for (uint8_t i = 1; i <= self->DoF; i++)
		self->q[i] = self->r[i] == 0;

	return NULL;
};

Error *IS_init_elements_length(InputScope *self) {
	uint16_t i;
	char errMsg[MAXL];

	// Aliases
	const vec3 *xyz = self->xyz;
	const int *N1 = self->N1;
	const int *N2 = self->N2;
	const float *rj = self->rj;

	for (i = 1; i <= self->nE; i++) {		/* calculate frame element lengths */
		const int n1 = N1[i];
		const int n2 = N2[i];
#define SQ(X) ((X)*(X))
		self->L[i] =	SQ( xyz[n2].x - xyz[n1].x ) +
				SQ( xyz[n2].y - xyz[n1].y ) +
				SQ( xyz[n2].z - xyz[n1].z );
#undef SQ

		self->L[i] = sqrt( self->L[i] );
		self->Le[i] = self->L[i] - rj[n1] - rj[n2];
		if ( n1 == n2 || self->L[i] == 0.0 ) {
		   sprintf(errMsg,
			" Frame elements must start and stop at different nodes\n  frame element %d  N1= %d N2= %d L= %e\n   Perhaps frame element number %d has not been specified.\n  or perhaps the Input Data file is missing expected data.\n",
		   i, n1,n2, self->L[i], i);
		   return Error_new(60, errMsg);
		}
		if ( self->Le[i] <= 0.0 ) {
		   sprintf(errMsg, " Node  radii are too large.\n  frame element %d  N1= %d N2= %d L= %e \n  r1= %e r2= %e Le= %e \n",
		   i, n1,n2, self->L[i], rj[n1], rj[n2], self->Le[i] );
		   return Error_new(61, errMsg);
		}
	}
	return NULL;
};
