#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
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
	double* delta = malloc(N * sizeof(*delta));
	TridiagonalSystem* system = create_system(N);

	if (!x || !delta || !system) {
		free(x);
		free(delta);
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
		build_system(system, x);
		tempo_jacobiana = timestamp() - tempo_jacobiana;
		LIKWID_MARKER_STOP("Jacobiana");

		if (calcula_norma(system->rhs, N) < epsilon) break;

		LIKWID_MARKER_START("Sistema Linear");
		tempo_sistema = timestamp() - tempo_sistema;
		solve_system(system, delta);
		tempo_sistema = timestamp() - tempo_sistema;
		LIKWID_MARKER_STOP("Sistema Linear");

		// atualiza x da próxima iteração
		for (int j = 0; j < N; ++j) x[j] += delta[j];

#ifndef BENCHMARKING
		if (i != 0) fputs("#\n", out);
		print_vec(x, N, out);
#endif

		if (calcula_norma(delta, N) < epsilon) break;
	}
	tempo_total = timestamp() - tempo_total;
	
	LIKWID_MARKER_STOP("Total");
	LIKWID_MARKER_CLOSE;

	fprintf(out, "###########\n# Tempo Total: %lf\n# Tempo Jacobiana: %lf\n# Tempo SL: %lf\n###########\n", tempo_total, tempo_jacobiana, tempo_sistema);
	
	free(x);
	free(delta);
	free_system(system);
	
	if (out != stdout) {
		fclose(out);
	}

	return 0;
}
