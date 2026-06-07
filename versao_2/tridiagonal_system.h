
typedef struct {

	int n;
	double* main;
	double* rhs;

} TridiagonalSystem;


void free_system(TridiagonalSystem*);
TridiagonalSystem* create_system(int);


// preenche o sistema linear resultante da linearização do sistema
// não linear de Broyden dado o vetor de entradas x
// ao mesmo tempo, calcula a norma do rhs e retorna
double build_system_compute_norm(TridiagonalSystem*, double*);

// resolve um sistema tridiagonal e soma o resultado no vetor passado
// ao mesmo tempo, calcula a norma do deslocamento e retorna
double solve_system_compute_norm(TridiagonalSystem*, double*);