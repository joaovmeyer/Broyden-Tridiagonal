
typedef struct {

	int n;
	double* upper;
	double* main;
	double* lower;
	double* rhs;

} TridiagonalSystem;


void free_system(TridiagonalSystem*);
TridiagonalSystem* create_system(int);


// preenche o sistema linear resultante da linearização do sistema
// não linear de Broyden dado o vetor de entradas x
void build_system(TridiagonalSystem*, double*);

// resolve um sistema tridiagonal e armazena o resultado
void solve_system(TridiagonalSystem*, double*);