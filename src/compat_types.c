#include "compat_types.h"
#include "NRutil.h"
#include "HPGutil.h"
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
