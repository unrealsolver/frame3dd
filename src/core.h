#pragma once

#include <stdint.h>

#include "compat_types.h"
#include "types.h"

uint8_t solve(
	const InputScope scope,
	const RuntimeArgs args,
	Results *results,
	SolverContext ctx,
	ResultScope rs,
	const int lc // Load case id
);

