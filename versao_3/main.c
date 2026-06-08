#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include "helper.h"
#include "../utils.h"
#include <math.h>

#define BENCHMARKING




int main() {

	int N, MAX;
	double x0, epsilon;
	
	if (scanf("%d %lf %lf %d", &N, &x0, &epsilon, &MAX) != 4) {
		puts("Entrada inválida.");
		return 1;
	}

	double* x = malloc(N * sizeof(*x));
	double* main = malloc(N * sizeof(*main));
	double* rhs = malloc(N * sizeof(*rhs));

	if (!x || !main || !rhs) {
		free(x);
		free(main);
		free(rhs);
		return 1;
	}


	rtime_t tempo_total = timestamp();
	

	for (int iter = 0; iter < MAX; ++iter) {

		double norma1 = -1.0;
		double norma2 = -1.0;

		/***********************************
		*  primeira iteração fora do loop  *
		***********************************/

		// acha a solução
		double delta = (iter == 0 ? epsilon : rhs[0]); // epsilon para não convergir na primeira iteração
		x[0] = (iter == 0 ? x0 : x[0] + delta); // x[0] = x0 na primeira iteração

		// calcula a norma
		double abs = fabs(delta);
		norma1 = norma1 > abs ? norma1 : abs;

		// monta o sistema
		rhs[0] = (2.0 * x[0] - 3.0) * x[0] - 1.0; // falta somar 2 * x[1]
		main[0] = -4.0 * x[0] + 3.0;

		// triangulariza (inferior = 0)
		double inv = 1.0 / main[0];
		main[0] = -2.0 * inv;

		for (int i = 1; i < N; ++i) {

			// acha a solução
			delta = (iter == 0 ? epsilon : rhs[i] - delta * main[i]);
			x[i] = (iter == 0 ? x0 : x[i] + delta);

			// calcula a norma
			abs = fabs(delta);
			norma1 = norma1 > abs ? norma1 : abs;
			rhs[i - 1] += 2.0 * x[i]; // termina de montar a última linha
			abs = fabs(rhs[i - 1]);
			norma2 = norma2 > abs ? norma2 : abs;

			// monta o sistema
			rhs[i] = (2.0 * x[i] - 3.0) * x[i] + x[i - 1] - 1.0; // falta somar 2 * x[i + 1]
			main[i] = -4.0 * x[i] + 3.0;

			// triangulariza (inferior = 0)
			rhs[i - 1] = (rhs[i - 1] - (i > 1 ? rhs[i - 2] * -1.0 : 0.0)) * inv;
			inv = 1.0 / (main[i] - main[i - 1] * -1.0);
			main[i] = -2.0 * inv;
		}

		// termina de triangularizar (tira o -1.0 da última entrada do rhs)
		rhs[N - 1] = ((rhs[N - 1] + 1.0) - rhs[N - 2] * -1.0) * inv;
		abs = fabs(rhs[N - 1]);
		norma2 = norma2 > abs ? norma2 : abs;

		if (iter != 0) {

			++iter;

#ifndef BENCHMARKING
			puts(""); print_vec(x, N); printf("#");
#endif
		}

		if ((norma1 < norma2 ? norma1 : norma2) < epsilon || iter >= MAX) break;

		// volta finalizando essa iteração e começando a próxima
		norma1 = -1.0;
		norma2 = -1.0;

		// acha a solução
		delta = rhs[N - 1];
		x[N - 1] += delta;

		// calcula a norma
		abs = fabs(delta);
		norma1 = norma1 > abs ? norma1 : abs;

		// monta o sistema
		rhs[N - 1] = (2.0 * x[N - 1] - 3.0) * x[N - 1]; // falta somar 1 * x[N - 2]
		main[N - 1] = -4.0 * x[N - 1] + 3.0;

		// triangulariza (superior = 0)
		inv = 1.0 / main[N - 1];
		main[N - 1] = -1.0 * inv;

		for (int i = N - 2; i >= 0; --i) {

			// acha a solução
			delta = rhs[i] - delta * main[i];
			x[i] += delta;

			// calcula a norma
			abs = fabs(delta);
			norma1 = norma1 > abs ? norma1 : abs;
			rhs[i + 1] += x[i];
			abs = fabs(rhs[i + 1]);
			norma2 = norma2 > abs ? norma2 : abs;

			// monta o sistema
			rhs[i] = (2.0 * x[i] - 3.0) * x[i] + 2.0 * x[i + 1] - 1.0; // falta somar 2 * x[i + 1]
			main[i] = -4.0 * x[i] + 3.0;

			// triangulariza (superior = 0)
			rhs[i + 1] = (rhs[i + 1] - (i < N - 2 ? rhs[i + 2] * -2.0 : 0.0)) * inv;
			inv = 1.0 / (main[i] - main[i + 1] * -2.0);
			main[i] = -1.0 * inv;
		}

		// termina de triangularizar
		rhs[0] = (rhs[0] - rhs[1] * -2.0) * inv;
		abs = fabs(rhs[0]);
		norma2 = norma2 > abs ? norma2 : abs;

#ifndef BENCHMARKING
		puts(""); print_vec(x, N); printf("#");
#endif

		if ((norma1 < norma2 ? norma1 : norma2) < epsilon) break;
	}


	tempo_total = timestamp() - tempo_total;

	printf("##########\n# Tempo Total: %lf\n###########\n", tempo_total);
	
	free(x);
	free(main);
	free(rhs);

	return 0;
}
