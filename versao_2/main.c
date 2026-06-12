#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include "tridiagonal_system.h"
#include "../helper.h"
#include "../utils.h"

#include <likwid.h>

// #define BENCHMARKING




int main(int argc, char* argv[]) {
	
	FILE* out = stdout;
	if (argc == 3 && (strcmp(argv[1], "-o") == 0)) {
		out = fopen(argv[2], "w");
		if (!out) return 1;
	} else if (argc != 1) {
		puts("Uso inválido.");
		return 1;
	}

	int N, MAX;
	double x0, epsilon;
	
	if (scanf("%d %lf %lf %d", &N, &x0, &epsilon, &MAX) != 4) {
		fputs("Entrada inválida.\n", out);
		return 1;
	}

	double* x = malloc(N * sizeof(*x));
	TridiagonalSystem* system = create_system(N);

	if (!x || !system) {
		free(x);
		free_system(system);
		return 1;
	}


	rtime_t tempo_total, tempo_jacobiana = 0, tempo_sistema = 0;


	// preenche aproximação inicial com x0
	for (int i = 0; i < N; ++i) x[i] = x0;

	LIKWID_MARKER_INIT;
	LIKWID_MARKER_REGISTER("Total");
	LIKWID_MARKER_REGISTER("Jacobiana");
	LIKWID_MARKER_REGISTER("Sistema Linear");
	
	LIKWID_MARKER_START("Total");
	
	tempo_total = timestamp();
	for (int i = 0; i < MAX; ++i) {

		LIKWID_MARKER_START("Jacobiana");
		tempo_jacobiana = timestamp() - tempo_jacobiana;
		if (build_system_compute_norm(system, x) < epsilon) break;
		tempo_jacobiana = timestamp() - tempo_jacobiana;
		LIKWID_MARKER_STOP("Jacobiana");

		LIKWID_MARKER_START("Sistema Linear");
		tempo_sistema = timestamp() - tempo_sistema;
		double norma = solve_system_compute_norm(system, x);
		tempo_sistema = timestamp() - tempo_sistema;
		LIKWID_MARKER_STOP("Sistema Linear");

#ifndef BENCHMARKING
		if (i != 0) fputs("#\n", out);
		print_vec(x, N, out);
#endif

		if (norma < epsilon) break;
	}
	tempo_total = timestamp() - tempo_total;
	
	LIKWID_MARKER_STOP("Total");
	LIKWID_MARKER_CLOSE;

	fprintf(out, "###########\n# Tempo Total: %lf\n# Tempo Jacobiana: %lf\n# Tempo SL: %lf\n###########\n", tempo_total, tempo_jacobiana, tempo_sistema);
	
	free(x);
	free_system(system);

	return 0;
}
