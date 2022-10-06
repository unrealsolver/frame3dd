#ifndef F3DD_TYPES_H
#define F3DD_TYPES_H

#include <stdbool.h>
#include "microstran/vec3.h"

/* ####### LOAD CASE DATA ####### */
/* Point load */
typedef struct {
	vec3 force;
	vec3 momentum;
	unsigned node_id;
} PointLoad;

typedef struct {
	PointLoad *data;
	unsigned size;
} PointLoads;

/* Uniform load */
typedef struct {
	vec3 force;
	unsigned edge_id;
} UniformLoad;

typedef struct {
	unsigned size;
	UniformLoad* data;
} UniformLoads;

typedef struct {
	PointLoads point;
	UniformLoads uniform;
} Loads;

typedef struct {
	vec3 gravity;
	Loads loads;
} LoadCase;

typedef struct {
	unsigned size;
	LoadCase* data;
} LoadCases;

/* ####### FRAME DATA ####### */
/* Nodes */
typedef struct {
	bool x;
	bool y;
	bool z;
	bool xx;
	bool yy;
	bool zz;
} DepthOfFreedom;

typedef struct {
	DepthOfFreedom dof;
	vec3 position;
	double radius;
} Node;

typedef struct {
	unsigned size;
	/* Some cells might be uninitialized for O(1) search */
	Node *data;
} Nodes;

/* Edges */
typedef struct {
	float density;
	/* Young's modulus */
	float E;
	/* shear modulus */
	float G;
} Material;

typedef struct {
	/* cross section area of each element*/
	float Ax;
	/* shear area in local y direction*/
	float Asy;
	/* shear area in local z direction*/
	float Asz;
	/* torsional moment of inertia */
	float Jx;
	/* bending moment of inertia about y-axis */
	float Iy;
	/* bending moment of inertia about z-axis */
	float Iz;
} Profile;

typedef struct {
	// TODO Add edge text label
	/* Start node id */
	unsigned start;
	/* End node id */
	unsigned end;
	/* Rotation of the element along the axis */
	float roll;
	/* Element's material */
	Material *material;
	/* Element's profile data */
	Profile *profile;
} Edge;

typedef struct {
	/* Number of elements */
	unsigned size;
	/* Some cells might be uninitialized for O(1) search */
	Edge *data;
} Edges;

typedef struct {
	Nodes nodes;
	Edges edges;
} Frame;

/* ####### SIMULATION OPTIONS ####### */
/* Options */
// FIXME Refactor
typedef struct {
	// TODO Review required
        // "-1" means no internal frame forces calcilatinons
	float dx; /**< x-increment for internal force data */
	// TODO Review required
	bool shear; /**< Enable shear deformation modeling */
	// TODO Review required
	bool nonlinear; /**< Enable geometrical nonlinearity modeling */
} SimulationOptions;

typedef struct {
	// TODO Review required
	double exagg_static; /**< Exaggerate display */
	double exagg_modal; /**< Exaggerate display */
	// TODO Review required
	float scale; /**< zoom scale for 3D plotting in Gnuplot */
	float pan; /**< >0: pan during animation; 0: don't */
} VisualizationOptions;

typedef struct {
	SimulationOptions simulation;
	VisualizationOptions visual;
} RunOptions;

#endif /* F3DD_TYPES_H */
