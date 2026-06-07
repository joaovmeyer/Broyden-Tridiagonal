#include "tridiagonal_system.h"
#include <stdlib.h>




void free_system(TridiagonalSystem* system) {
	if (!system) return;
	free(system->upper);
	free(system->main);
	free(system->lower);
	free(system->rhs);
	free(system);
}

TridiagonalSystem* create_system(int n) {
	TridiagonalSystem* system = malloc(sizeof(*system));
	if (!system) return NULL;

	system->n = n;
	system->upper = malloc((n - 1) * sizeof(*system->upper));
	system->main = malloc(n * sizeof(*system->main));
	system->lower = malloc((n - 1) * sizeof(*system->lower));
	system->rhs = malloc(n * sizeof(*system->rhs));

	if (!system->upper || !system->main || !system->lower || !system->rhs) {
		free_system(system);
		return NULL;
	}

	return system;
}



void build_system(TridiagonalSystem* system, double* x) {
	// F_1(x) = -2x_1^2 + 3x_1 - 2x_2 + 1
	// F_i(x) = -2x_i^2 + 3x_i - x_{i-1} - 2x_{i+1} + 1
	// F_n(x) = -2x_n^2 + 3x_n - x_{n-1} ( + 1?)

	if (system->n == 1) {
		system->rhs[0] = (2.0 * x[0] - 3.0) * x[0]; // com ou sem + 1 ???
		system->main[0] = -4.0 * x[0] + 3.0;
		return;
	}

	// primeira linha
	system->rhs[0] = (2.0 * x[0] - 3.0) * x[0] + 2.0 * x[1] - 1.0; // rhs = -F
	system->upper[0] = -2.0;
	system->main[0] = -4.0 * x[0] + 3.0;

	// meio
	int i;
	for (i = 1; i < system->n - 1; ++i) {
		system->rhs[i] = (2.0 * x[i] - 3.0) * x[i] + x[i - 1] + 2.0 * x[i + 1] - 1.0; // rhs = -F
		system->upper[i] = -2.0;
		system->main[i] = -4.0 * x[i] + 3.0;
		system->lower[i - 1] = -1.0;
	}

	// última linha (i == n - 1)
	system->rhs[i] = (2.0 * x[i] - 3.0) * x[i] + x[i - 1]; // rhs = -F
	system->main[i] = -4.0 * x[i] + 3.0;
	system->lower[i - 1] = -1.0;

}


// método de Thomas é basicamente eliminação gaussiana, mas simplificando o código
// considerando o fato de que nossa matriz tem banda 3. Sem pivoteamento, pode
// ter problemas se a diagonal for próxima de 0, mas fazer o pivoteamento pode
// destruir a estrutura tridiagonal. Por isso, só ignoramos esse problema
void solve_system(TridiagonalSystem* system, double* ans) {

	if (system->n == 1) {
		ans[0] = system->rhs[0] / system->main[0];
		return;
	}

	double inv = 1.0 / system->main[0];
	system->upper[0] *= inv;
	system->rhs[0] *= inv;

	// triangularização
	int i = 1;
	for (i = 1; i < system->n - 1; ++i) {
		inv = 1.0 / (system->main[i] - system->upper[i - 1] * system->lower[i - 1]);
		system->rhs[i] = (system->rhs[i] - system->rhs[i - 1] * system->lower[i - 1]) * inv;
		system->upper[i] *= inv;
	}

	inv = 1.0 / (system->main[i] - system->upper[i - 1] * system->lower[i - 1]);
	system->rhs[i] = (system->rhs[i] - system->rhs[i - 1] * system->lower[i - 1]) * inv;

	// retro-substituição
	ans[system->n - 1] = system->rhs[system->n - 1];
	for (int i = system->n - 2; i >= 0; --i) {
		ans[i] = system->rhs[i] - ans[i + 1] * system->upper[i];
	}
}