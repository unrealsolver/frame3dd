#ifndef F3DD_TYPES_H
#define F3DD_TYPES_H

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

struct InputData { };

#endif /* F3DD_TYPES_H */
