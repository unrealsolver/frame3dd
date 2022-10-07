#include <math.h>
#include <stdint.h>
#include "struct_writer.h"

void write_static_struct(
	LoadCaseResult *lc_result, Frame *frame, double *D
) {
	// TODO write 'ok'

	for (uint16_t j = 0; j < frame->nodes.size; j++) {
		NodeDisplacement *disp = &lc_result->displacements[j];
		disp->pos.x = D[6 * j + 1];
		disp->pos.y = D[6 * j + 2];
		disp->pos.z = D[6 * j + 3];
		disp->rot.x = D[6 * j + 4];
		disp->rot.y = D[6 * j + 5];
		disp->rot.z = D[6 * j + 6];
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
