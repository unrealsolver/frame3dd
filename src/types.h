#ifndef F3DD_TYPES_H
#define F3DD_TYPES_H

#include <stdbool.h>
#include "microstran/vec3.h"

typedef struct {
	vec3 force;
	vec3 momentum;
	unsigned node_id;
} PointLoad;

typedef struct {
	PointLoad *data;
	unsigned size;
} PointLoads;

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
} LoadcaseData;

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
	Node *data;
} Nodes;

typedef struct {
	unsigned start_node_id;
	unsigned end_node_id;
} Edge;

typedef struct {
	unsigned size;
	Edge *data;
} Edges;

typedef struct {
	Nodes nodes;
	Edges edges;
} Frame;

#endif /* F3DD_TYPES_H */
