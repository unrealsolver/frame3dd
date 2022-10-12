#ifndef COMPAT_TYPES_H
#define COMPAT_TYPES_H

#include <stdlib.h>

struct Overrides {
	double	exagg	= -1.0; // over-ride input file value
	double	tol	= -1.0; // over-ride input file value
	double	shift	= -1.0; // over-ride input file value
	float	pan	= -1.0; // over-ride input file value
	int8_t	shear	= -1;	// over-ride input file value
	int8_t	geom	= -1;	// over-ride input file value
	int8_t	anlyz	= -1;	// over-ride input file value
	int8_t	D3	= 0;	// over-ride 3D plotting check
	int8_t	lump	= -1;	// over-ride input file value
	int8_t	modal	= -1;	// over-ride input file value
	int8_t	condense= -1;	// over-ride input file value
};

typedef struct RuntimeArgs {
	struct	Overrides overrides;
	int8_t	write_matrix	= -1;	// write stiffness and mass matrix
	int8_t	axial_sign	= -1;	// suppress 't' or 'c' in output data
	int8_t	filetype	= 0;	// 1 if .CSV, 2 if file is Matlab
	int8_t	debug		= 0;	// 1: debugging screen output, 0: none
	int8_t	verbose		= 1;	// 1: copious screen output, 0: none
} RuntimeArgs;
#endif /* COMPAT_TYPES_H */
