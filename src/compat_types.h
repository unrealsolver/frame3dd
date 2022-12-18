/*
 *
*//** @file
	Compatibility types
 */
#ifndef COMPAT_TYPES_H
#define COMPAT_TYPES_H

#include <stdlib.h>
#include <stdint.h>
#include "frame3dd.h"
#include "HPGutil.h"

typedef struct {
	uint8_t code;
	const char *message;
} Error;

Error *Error_new(const uint8_t code, const char *message);

void Error_handle(Error *self);

void Error_destroy(Error *self);

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

typedef struct {
	double rms_resid; // root mean square of residual displ. error
	double error;	// rms equilibrium error and reactions
	int ok;		// number of (-ve) diag. terms of L D L'
	int iter;	// number of iterations
} SolverContext;

// TODO check which properties were never used
typedef struct {
	vec3	*xyz;		// X,Y,Z node coordinates (global)
	int
		*q,*r, sumR;	// reaction data, total no. of reactions
	int
		nN,		// number of Nodes
		nE,		// number of frame Elements
		nL,		// number of Load cases
		nX,		// number of elemts w/ extra mass
		nI,		// number of nodes w/ extra inertia
		nM,		// number of desired modes
		nD[_NL_],	// number of prescribed nodal displ'nts
		nF[_NL_],	// number of loaded nodes
		nU[_NL_],	// number of members w/ unifm dist loads
		nW[_NL_],	// number of members w/ trapz dist loads
		nP[_NL_],	// number of members w/ conc point loads
		nT[_NL_],	// number of members w/ temp. changes
		DoF,		// number of Degrees of Freedom
		nR,		// number of restrained nodes
		nC,		// number of condensed nodes
		*N1, *N2,	// begin and end node numbers
		Mmethod,	// 1: Subspace Jacobi, 2: Stodola
		Cmethod,	// matrix condensation method
		Cdof,		// number of condensed degrees o freedom
		*c,		// vector of DoF's to condense
		*m,		// vector of modes to condense
		lump,		// 1: lumped, 0: consistent mass matrix
		shear,		// indicates shear deformation
		geom,		// indicates  geometric nonlinearity
		anlyz,		// 1: stiffness analysis, 0: data check
		anim[128];	// the modes to be animated
	float
		pan,		// >0: pan during animation; 0: don't
		scale,		// zoom scale for 3D plotting in Gnuplot
		dx,		// x-increment for internal force data
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
		tol,		// tolerance for modal convergence
		shift,		// shift-factor for rigid-body-modes
		struct_mass,	// mass of structural system
		total_mass,	// total structural mass and extra mass
		exagg_static,	// exaggerate static displ. in mesh data
		exagg_modal;	// exaggerate modal displ. in mesh data
} InputScope;

typedef struct {
	double
		**K,	// equilibrium stiffness matrix
		*R,	// total reaction force vector
		*D,	// displacement vector
		**Q;	// local member node end-forces
} ResultScope;

/**
 * Result Scope initialize for given InputScope
 */
void RS_init_for_IS(ResultScope *self, const InputScope *is);

/**
 * Initialize number of nodes nN and derived data (rj, xyz)
 */
void IS_set_nN(InputScope *self, const uint16_t nN);

/**
 * Initialize number of elements nE and derived data
 */
void IS_set_nE(InputScope *self, const uint16_t nE);

/**
 * Initialize number of elements nE and derived data
 */
void IS_set_nL(InputScope *self, const uint8_t nL);

/**
 * Initialize eqF_mech for given Load Case lc
 * Requires gravity loads to be set for given lc
 */
void IS_init_eqF_mech(InputScope *self, const uint8_t lc);

/**
 * assemble all element equivalent loads into
 * separate load vectors for mechanical and thermal loading
 */
void IS_assemble_eq_loads(InputScope *self, uint8_t lc);

/**
 * Initialize derived reactions data (q, sumR)
 */
Error *IS_init_reactions(
	InputScope *self
);

/**
 * Initialize L (length) and Le (effective length)
 * of all frame elements
 */
Error *IS_init_elements_length(
	InputScope *self
);

#endif /* COMPAT_TYPES_H */
