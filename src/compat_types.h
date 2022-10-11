#ifndef COMPAT_TYPES_H
#define COMPAT_TYPES_H

#include <stdlib.h>

typedef struct {
	int8_t	shear,	//   over-ride input file value
		geom,	//   over-ride input file value
		anlyz,	//   over-ride input file value
		D3,	//   over-ride 3D plotting check
		lump,	//   over-ride input file value
		modal,	//   over-ride input file value
		write,	//   write stiffness and mass matrix
		axial,	//   suppress 't' or 'c' in output data
		condense; // over-ride input file value
} OverrideFlags;
#endif /* COMPAT_TYPES_H */
