#ifndef COMPAT_TYPES_H
#define COMPAT_TYPES_H

#include <stdlib.h>
#include "frame3dd.h"

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

// TODO check which properties were never used
typedef struct {
	vec3	*xyz;		// X,Y,Z node coordinates (global)
	int
		*q,*r, sumR;	// reaction data, total no. of reactions
	int
		nN,		// number of Nodes
		nE,		// number of frame Elements
		nL,		// number of Load cases
		nD[_NL_],	// number of prescribed nodal displ'nts
		nF[_NL_],	// number of loaded nodes
		nU[_NL_],	// number of members w/ unifm dist loads
		nW[_NL_],	// number of members w/ trapz dist loads
		nP[_NL_],	// number of members w/ conc point loads
		nT[_NL_],	// number of members w/ temp. changes
		DoF,		// number of Degrees of Freedom
		nR,		// number of restrained nodes
		*N1, *N2;	// begin and end node numbers
	float
		*rj,		// node size radius, for finite sizes
		*Ax,*Asy, *Asz,	// cross section areas, incl. shear
		*Jx,*Iy,*Iz,	// section inertias
		*E, *G,		// elastic modulus and shear moduli
		*p,		// roll of each member, radians
		***U,		// uniform distributed member loads
		***W,		// trapizoidal distributed member loads
		***P,		// member concentrated loads
		***T,		// member temperature  loads
		**Dp,		// prescribed node displacements
		*d, *EMs,	// member densities and extra inertia
		*NMs,		// mass of a node
		*NMx,*NMy,*NMz,	// inertia of a node in global coord
		gX[_NL_],	// gravitational acceleration in global X
		gY[_NL_],	// gravitational acceleration in global Y
		gZ[_NL_];	// gravitational acceleration in global Z
	double
		*L,		// node-to-node length of each element
		*Le,		// effcve lngth, accounts for node size
		*F,	 	// total load vectors for a load case
		***eqF_mech,	// equivalent end forces from mech loads global
		***eqF_temp,	// equivalent end forces from temp loads global
		**F_mech,	// mechanical load vectors, all load cases
		**F_temp,	// thermal load vectors, all load cases
		**Q;		// local member node end-forces
} InputScope;

#endif /* COMPAT_TYPES_H */
