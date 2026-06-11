#include "helper.h"
#include <math.h>
#include <stdio.h>


// imprime os elementos do vetor com o prefixo "xi = " separados por "\n"
void print_vec(double* x, int n) {
    for (int i = 0; i < n; ++i) {
        printf("x%i = %f\n", i + 1, x[i]);
    }
}


// calcula norma do infinito (maior valor do vetor)
double calcula_norma(double* vec, int n) {
	double norma = -1.0;

	for (int i = 0; i < n; ++i) {
		double abs = fabs(vec[i]);
		norma = norma > abs ? norma : abs;
	}

	return norma;
}
