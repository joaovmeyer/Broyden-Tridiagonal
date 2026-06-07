#include "tridiagonal_system.h"
#include "helper.h"
#include <stdlib.h>
#include <math.h>




void free_system(TridiagonalSystem* system) {
	if (!system) return;
	free(system->main);
	free(system->rhs);
	free(system);
}

TridiagonalSystem* create_system(int n) {
	TridiagonalSystem* system = malloc(sizeof(*system));
	if (!system) return NULL;

	system->n = n;
	system->main = malloc(n * sizeof(*system->main));
	system->rhs = malloc(n * sizeof(*system->rhs));

	if (!system->main || !system->rhs) {
		free_system(system);
		return NULL;
	}

	return system;
}



double build_system_compute_norm(TridiagonalSystem* system, double* x) {
	// F_1(x) = -2x_1^2 + 3x_1 - 2x_2 + 1
	// F_i(x) = -2x_i^2 + 3x_i - x_{i-1} - 2x_{i+1} + 1
	// F_n(x) = -2x_n^2 + 3x_n - x_{n-1} ( + 1?)

	if (system->n == 1) {
		system->rhs[0] = (2.0 * x[0] - 3.0) * x[0]; // com ou sem + 1 ???
		system->main[0] = -4.0 * x[0] + 3.0;
		return fabs(system->rhs[0]);
	}

	double norma = -1.0;

	// primeira linha
	system->rhs[0] = (2.0 * x[0] - 3.0) * x[0] + 2.0 * x[1] - 1.0; // rhs = -F
	system->main[0] = -4.0 * x[0] + 3.0;

	double abs = fabs(system->rhs[0]);
	norma = norma > abs ? norma : abs;

	// meio
	int i = 1;
	#pragma omp simd reduction(max:norma)
	for (i = 1; i < system->n - 1; ++i) {
		system->rhs[i] = (2.0 * x[i] - 3.0) * x[i] + x[i - 1] + 2.0 * x[i + 1] - 1.0; // rhs = -F
		system->main[i] = -4.0 * x[i] + 3.0;

		abs = fabs(system->rhs[i]);
		norma = norma > abs ? norma : abs;
	}

	// última linha (i == n - 1)
	system->rhs[i] = (2.0 * x[i] - 3.0) * x[i] + x[i - 1]; // rhs = -F
	system->main[i] = -4.0 * x[i] + 3.0;

	abs = fabs(system->rhs[i]);
	norma = norma > abs ? norma : abs;

	return norma;
}


// método de Thomas é basicamente eliminação gaussiana, mas simplificando o código
// considerando o fato de que nossa matriz tem banda 3. Sem pivoteamento, pode
// ter problemas se a diagonal for próxima de 0, mas fazer o pivoteamento pode
// destruir a estrutura tridiagonal. Por isso, só ignoramos esse problema

// https://sistemas.fciencias.unam.mx/~ediaz/Cursos/Anal-numerico/Numerical%20mathematics.pdf (p. 110)
// https://research.nvidia.com/sites/default/files/pubs/2010-01_Fast-Tridiagonal-Solvers/Zhang_Fast_2009.pdf (GPU)
/*double solve_system_compute_norm(TridiagonalSystem* system, double* ans) {

	// triangularização
	for (int i = 0; i < system->n - 1; ++i) {
		double m = -1.0 / system->main[i];

		// L[i] <- L[i] - m * L[k]
		system->main[i + 1] -= m * -2.0;
		system->rhs[i + 1] -= m * system->rhs[i];
	}

	double norma = -1.0;

	// retro-substituição
	double delta = system->rhs[system->n - 1] / system->main[system->n - 1];
	ans[system->n - 1] += delta;

	double abs = fabs(delta);
	norma = norma > abs ? norma : abs;

	for (int i = system->n - 2; i >= 0; --i) {
		delta = (system->rhs[i] - delta * -2.0) / system->main[i];
		ans[i] += delta;

		abs = fabs(delta);
		norma = norma > abs ? norma : abs;
	}

	return norma;
}*/

double solve_system_compute_norm(TridiagonalSystem* system, double* ans) {

	if (system->n == 1) {
		double delta = system->rhs[0] / system->main[0];
		ans[0] += delta;
		return fabs(delta);
	}

	double inv = 1.0 / system->main[0];
	system->main[0] = -2.0 * inv;
	system->rhs[0] *= inv;

	// triangularização
	int i = 1;
	for (i = 1; i < system->n - 1; ++i) {
		inv = 1.0 / (system->main[i] - system->main[i - 1] * -1.0);
		system->rhs[i] = (system->rhs[i] - system->rhs[i - 1] * -1.0) * inv;
		system->main[i] = -2.0 * inv;
	}

	inv = 1.0 / (system->main[i] - system->main[i - 1] * -1.0);
	system->rhs[i] = (system->rhs[i] - system->rhs[i - 1] * -1.0) * inv;


	double norma = -1.0;

	// retro-substituição
	double delta = system->rhs[system->n - 1];
	ans[system->n - 1] += delta;

	double abs = fabs(delta);
	norma = norma > abs ? norma : abs;

	for (int i = system->n - 2; i >= 0; --i) {
		delta = system->rhs[i] - delta * system->main[i];
		ans[i] += delta;

		abs = fabs(delta);
		norma = norma > abs ? norma : abs;
	}

	return norma;
}