#ifndef COMPAT_TYPES_H
#define COMPAT_TYPES_H

#include <stdlib.h>

struct Overrides {
	double	exagg,
		tol,
		shift;
	float	pan;
	int8_t	shear,
		geom,
		anlyz,
		D3,
		lump,
		modal,
		condense;
};

typedef struct RuntimeArgs {
	struct	Overrides overrides; // Command line overrides
	int8_t	write_matrix,	// write stiffness and mass matrix
		axial_sign,	// suppress 't' or 'c' in output data
		debug,		// 1: debugging screen output, 0: none
		verbose;	// 1: copious screen output, 0: none
} RuntimeArgs;
#endif /* COMPAT_TYPES_H */
