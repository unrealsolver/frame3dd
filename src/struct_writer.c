#include <math.h>
#include <stdint.h>
#include "struct_writer.h"

void write_static_struct(
	LoadCaseResult *lc_result, Frame *frame,
	const double *D, const double **Q
) {
	// TODO write 'ok'

	// Write Node Displacements
	double abs_disp;
	for (uint16_t j = 0; j < frame->nodes.size; j++) {
		abs_disp = 0.0;
		NodeDisplacement *disp = &lc_result->displacements[j];
		/*
		// Calculate total absolute displacement
		for (uint8_t i = 1; i <= 6; i++) abs_disp += fabs(D[6 * j + i]);
		// FIXME this will cause existance of unitialized ram areas
		// Skip iteration if there are no node displacements
		if (abs_disp == 0.0) continue; // or set NULL
		*/
		disp->pos.x = D[6 * j + 1];
		disp->pos.y = D[6 * j + 2];
		disp->pos.z = D[6 * j + 3];
		disp->rot.x = D[6 * j + 4];
		disp->rot.y = D[6 * j + 5];
		disp->rot.z = D[6 * j + 6];
	}

	// Write Element End Forces
	for (uint16_t j = 0; j < frame->edges.size; j++) {
		const Edge edge = frame->edges.data[j];
		// Convert to 'natural' ordering (starting w/ 1)
		const uint16_t _j = j + 1;
		EdgeResult *er = &lc_result->edges[j];
		// For each edge ending (start and the end)
		for (uint8_t ending = 0; ending < 2; ending++) {
			const uint8_t idx_offset = 6 * ending;
			ElementEndForces *forces = &er->forces[ending];
			forces->Nx = Q[_j][1 + idx_offset];
			forces->Nx_rot = forces->Nx >= 0.0 /* TODO && axial_sign */ ? 'c' : 't';
			forces->Vy = Q[_j][2 + idx_offset];
			forces->Vz = Q[_j][3 + idx_offset];
			forces->Txx = Q[_j][4 + idx_offset];
			forces->Myy = Q[_j][5 + idx_offset];
			forces->Mzz = Q[_j][6 + idx_offset];
		}
	}
}

void test_print(Frame *frame, LoadCases *load_cases, Results *results) {
	for (int j = 0; j < load_cases->size; j++) {
		printf("Load Case %d\n", j);
		LoadCaseResult *lc_result = &results->data[j];
		for (int i = 0; i < frame->nodes.size; i++) {
			NodeDisplacement *disp = &lc_result->displacements[i];
			printf(
				"%d:\tNodeDisplacement(x=%.2f, y=%.2f, z=%.2f, xx=%.2f, yy=%.2f, zz=%.2f, )\n",
				i, disp->pos.x, disp->pos.y, disp->pos.z, disp->rot.x, disp->rot.y, disp->rot.z
			);
		}
	}
}
