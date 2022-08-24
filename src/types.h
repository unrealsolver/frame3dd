#ifndef F3DD_TYPES_H
#define F3DD_TYPES_H

#include <stdbool.h>
#include "microstran/vec3.h"

typedef struct {
	unsigned int node_id;
	vec3 force;
} UniformLoad;

typedef struct {
	unsigned int size;
	UniformLoad* data;
} UniformLoads;

typedef struct {
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
	unsigned int size;
	Node *data;
} Nodes;

typedef struct {
	unsigned int start_node_id;
	unsigned int end_node_id;
} Edge;

typedef struct {
	unsigned int size;
	Edge *data;
} Edges;

typedef struct {
	Nodes nodes;
	Edges edges;
} Frame;

#endif /* F3DD_TYPES_H */
