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

Error *IS_set_derived_reaction_data(InputScope *self, const uint_fast8_t *reactions) {
	char errMsg[MAXL];
	self->q = ivector(1, self->DoF);	/* allocate memory for reaction data ... */

	// Calculate total number of restrained nodes
	self->sumR=0;
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
