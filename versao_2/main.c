#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include "tridiagonal_system.h"
#include "helper.h"
#include "../utils.h"

#define BENCHMARKING




int main() {

	int N, MAX;
	double x0, epsilon;
	
	if (scanf("%d %lf %lf %d", &N, &x0, &epsilon, &MAX) != 4) {
		puts("Entrada inválida.");
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

	tempo_total = timestamp();
	for (int i = 0; i < MAX; ++i) {

		tempo_jacobiana = timestamp() - tempo_jacobiana;
		if (build_system_compute_norm(system, x) < epsilon) break;
		tempo_jacobiana = timestamp() - tempo_jacobiana;

		tempo_sistema = timestamp() - tempo_sistema;
		double norma = solve_system_compute_norm(system, x);
		tempo_sistema = timestamp() - tempo_sistema;

#ifndef BENCHMARKING
		puts(""); print_vec(x, N); printf("#");
#endif

		if (norma < epsilon) break;
	}
	tempo_total = timestamp() - tempo_total;

	printf("##########\n# Tempo Total: %lf\n# Tempo Jacobiana: %lf\n# Tempo SL: %lf\n###########\n", tempo_total, tempo_jacobiana, tempo_sistema);
	
	free(x);
	free_system(system);

	return 0;
}