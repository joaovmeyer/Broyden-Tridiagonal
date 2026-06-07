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

	tempo_total = timestamp();
	for (int i = 0; i < MAX; ++i) {

		tempo_jacobiana = timestamp() - tempo_jacobiana;
		build_system(system, x);
		tempo_jacobiana = timestamp() - tempo_jacobiana;

		if (calcula_norma(system->rhs, N) < epsilon) break;

		tempo_sistema = timestamp() - tempo_sistema;
		solve_system(system, delta);
		tempo_sistema = timestamp() - tempo_sistema;

		// atualiza x da próxima iteração
		for (int j = 0; j < N; ++j) x[j] += delta[j];

#ifndef BENCHMARKING
		puts(""); print_vec(x, N); printf("#");
#endif

		if (calcula_norma(delta, N) < epsilon) break;
	}
	tempo_total = timestamp() - tempo_total;

	printf("##########\n# Tempo Total: %lf\n# Tempo Jacobiana: %lf\n# Tempo SL: %lf\n###########\n", tempo_total, tempo_jacobiana, tempo_sistema);
	
	free(x);
	free(delta);
	free_system(system);

	return 0;
}