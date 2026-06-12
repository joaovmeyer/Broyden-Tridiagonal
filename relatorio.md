# Resolução de Sistemas Não-Lineares de Broyden

## Introdução
Sistemas não-lineares de Broyden são os sistemas da forma **F**(**X**) = **0**, com 
$$
f_1(\bold{X}) = -2x_1^2 + 3x_1 - 2x_2 + 1
$$ $$
f_i(\bold{X}) = -2x_i^2 + 3x_i - x_{i-1} - 2x_{i+1} + 1
$$ $$
f_n(\bold{X}) = -2x_n^2 + 3x_n - x_{n-1}
$$
Para a solução do sistema não-linear, foi empregado o Método de Newton. 
### Método de Newton
O Método de Newton é um método iterativo que, dada uma aproximação inicial, lineariza o problema usando sua aproximação de Taylor de primeira ordem e usa a solução do sistema linearizado como próxima aproximação. Com isso, concluímos que dado uma aproximação $\bold{X}_k$, a próxima aproximação $\bold{X}_{k+1}$ satisfaz
$$
\bold{F}(\bold{X}_{k+1}) \approx \bold{F}(\bold{X}_k) + \bold{J}_\bold{F}(\bold{X}_k) (\bold{X}_{k+1} - \bold{X}_k) = \bold{0}
$$
Ou seja, temos um sistema linear da forma
$$
\bold{J}_\bold{F}(\bold{X}_k)\bold{\Delta}_k = -\bold{F}(\bold{X}_k)
$$
Encontrado $\bold{\Delta}_k$, recuperamos a incógnita $\bold{X}_{k+1} = \bold{\Delta}_k + \bold{X}_k$.

### Cálculo da Jacobiana e Resolução dos Sistemas Lineares
Note que o cálculo da matriz jacobiana para o problema apresentado é relativamente simples. Cada função $f_i$ depende no máximo de $x_{i-1}$, $x_i$ e $x_{i+1}$. Com isso, podemos perceber que a jacobiana assume uma forma tridiagonal. Isto é, a jacobiana só assume valores diferentes de 0 na sua diagonal principal e nas entradas diretamente acima/abaixo dela, enquanto todos os outros valores são nulos. Essa estrutura nos permite resolver o sistema linear que aparece durante as iterações do Método de Newton de forma eficiente. Neste trabalho, foi utilizado o Algoritmo de Thomas, que nada mais é do que um caso simplificado da Eliminação Gaussiana que toma proveito (computacionalmente falando) da estrutura tridiagonal da matriz. Além disso, devido à natureza sequencial do algoritmo de Thomas, foi explorado o algoritmo da Redução Cíclica, que embora efetue mais operações, permite que parte disso seja feito sem sincronia, permitindo maior aproveitamento do pipelining e registradores SIMD do processador.

Em geral, o objetivo deste trabalho é não apenas desenvolver um programa que encontre uma solução para o sistema não-linear de Broyden utilizando o Método de Newton, mas também otimizar o programa de modo que seu tempo de execução seja mínimo. Para fins comparativos, será apresentado uma versão base, sem otimizações radicais do código. Em seguida, serão introduzidas outras versões, cada uma com o objetivo de superar a performance da última. Juntamente, serão apresentados resultados qualitativos de cada versão e análises de seus desempenhos e potenciais problemas. Ao final do trabalho, devemos ter um programa que resolve o sistema não-linear de Broyden o mais rápido possível, e com cada otimização justificada com base em conceitos teóricos. 

## Uso de IAs
Declaro que o uso de modelos de linguagem (Chat GPT etc.) esteve presente principalmente nas partes iniciais deste trabalho. Em uma etapa de "brainstorming", chatbots foram utilizados para direcionar ideias e revisar possibilidades. Quanto à geração de código, declaro que todo código entregue foi escrito exclusivamente por mim. Nesta etapa, chatbots foram usados exclusivamente para esclarecer dúvidas simples sobre a linguagem C (cuja qual eu não tenho costume de utilizar) e suas melhores práticas. Algumas vezes, também foram utilizados para verificar a correção do código ou detecção prematura de bugs (por mais que não muito eficazes em encontrar bugs complexos, se mostraram de grande ajuda para corrigir problemas simples de desatenção). 
Declaro ainda que para análise dos resultados e escrita do relatório, não houve uso algum de ferramentas de IA.

## Versão Base

Para versão inicial foi utilizada uma estrutura para armazenar o sistema linear resultante da linearização do Método de Newton que apenas armazena as 3 diagonais necessárias da matriz do sistema, além do vetor do sistema e seu tamanho:
```cpp
typedef struct {
	int n;
	double* upper;
	double* main;
	double* lower;
	double* rhs;
} TridiagonalSystem;
```
Esta representação é fundamental para que seja possível armazenar sistemas grandes (N = 20000) sem esgotar a memória do computador. Além disso, otimiza o uso dos caches, que não terão que armazenar valores de 0 da matriz. Para preencher o sistema linear em cada iteração do Método de Newton, foi utilizada uma função que simplesmente calcula o vetor $-\bold{F}(\bold{X}_k)$ e sua jacobiana, dado o vetor $\bold{X}_k$:
```cpp
void build_system(TridiagonalSystem* system, double* x) {

	// primeira linha
	system->rhs[0] = (2.0 * x[0] - 3.0) * x[0] + 2.0 * x[1] - 1.0; // rhs = -F
	system->upper[0] = -2.0;
	system->main[0] = -4.0 * x[0] + 3.0;

	// meio
	int i;
	for (i = 1; i < system->n - 1; ++i) {
		system->rhs[i] = (2.0 * x[i] - 3.0) * x[i] + x[i - 1] + 2.0 * x[i + 1] - 1.0;
		system->upper[i] = -2.0;
		system->main[i] = -4.0 * x[i] + 3.0;
		system->lower[i - 1] = -1.0;
	}

	// última linha (i == n - 1)
	system->rhs[i] = (2.0 * x[i] - 3.0) * x[i] + x[i - 1]; // rhs = -F
	system->main[i] = -4.0 * x[i] + 3.0;
	system->lower[i - 1] = -1.0;
}
```

Para solução dos sistemas lineares gerados, foi implementado o Algoritmo de Thomas *in-place*, que modifica o sistema passado e o triangulariza:
```cpp
void solve_system(TridiagonalSystem* system, double* ans) {

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
```
Com isso e algumas outras funções auxiliares, foi possível escrever o loop principal do método de Newton:
```cpp
	for (int i = 0; i < MAX; ++i) {

		build_system(system, x);
		if (calcula_norma(system->rhs, N) < epsilon) break;

		solve_system(system, delta);

		// atualiza x da próxima iteração
		for (int j = 0; j < N; ++j) x[j] += delta[j];

		puts(""); print_vec(x, N); printf("#");

		if (calcula_norma(delta, N) < epsilon) break;
	}
```
Aqui, `calcula_norma` retorna o maior elemento em módulo do vetor passado, de modo que os critérios de parada se alinhem com a versão do Método de Newton mostrada durante as aulas.

Abaixo, estão apresentados gráficos mensurando a performance desta versão quando considerados diferentes tamanhos do sistema linear:

Todos os gráficos (no total, 12) podem ser encontrados em outro lugar.

Embora o programa seja estruturado de uma maneira decentemente eficiente, onde cada parte lida bem com seu problema isolado, ele não considera o problema como um todo. Em especial, podemos observar que nossa solução está limitada pela velocidade da memória (o que causa o efeito de "escada" no gráfico da performance, cada vez que o problema deixa de caber em uma área mais rápida da memória). Minimizar acessos desnecessários de memória será o principal desafio da próxima versão.

## Versão 2
É definido como o objetivo da versão 2 minimizar as "passagens por vetores". Isso é uma simplificação, já que, por exemplo, ler um valor e depois sobrescreve-lo na memória é considerado igual a somente ler ou escrever, que tecnicamente faz menos acessos à memória.
Outro objetivo é reduzir o número de vetores utilizados. Cada vetor desnecessário ocupa espaço nas caches, que faz com que o padrão de "escada" observado na versão anterior tenha "degraus" menores, resultando em um mau aproveitamento das regiões mais rápidas da memória.
Primeiramente, observe que temos no total 6 vetores sendo utilizados na versão base: 3 para armazenar $\bold{J_F}(\bold{X}_k)$, 1 para armazenar $\bold{F}(\bold{X}_k)$, outro para $\bold{\Delta}_k$ e outro para $\bold{X}_k$. Em termos de passagens por esses vetores, temos um total de 15 passagens:
 - 4 para montar o sistema
 - 1 para avaliar o primeiro critério de convergência
 - 4 para triangularizar
 - 3 para retro-substituição
 - 2 para o cálculo de $\bold{X}_{k+1}$ a partir do $\bold{\Delta}_k$
 - 1 para avaliar o segundo critério de convergência

A otimização mais óbvia quando olhamos para o código é eliminar os vetores que armazenam a subdiagonal e superdiagonal, pois estas são constantes e iguais a -1 e -2, respectivamente. Note que na versão anterior, passávamos por essas diagonais uma vez cada para montar o sistema linear, uma vez cada para fazer a triangularização, e uma passagem adicional da superdiagonal durante a retro-substituição. No entanto, como a triangulação aplicada modifica a superdiagonal, não conseguimos eliminar esta última passagem. Para evitar a necessidade de um vetor adicionar para armazenar a superdiagonal modificada, usamos o espaço "sobrando" do vetor da diagonal principal, uma vez que essa fica constante igual a 1 após a triangularização. Após eliminarmos essas passagens desnecessárias, ficamos com 11 passagens por vetores e um total de 4 vetores.
Podemos eliminar um vetor adicional alterando a funcionalidade da função que resolve o sistema linear: invés de sobrescrever o resultado, somamos ele no vetor de entrada. Com essa simples mudança de funcionalidade, podemos passar o vetor de $\bold{X}_k$ e obteremos diretamente $\bold{X}_k + \bold{\Delta}_k = \bold{X}_{k+1}$. Com isso, eliminamos as duas passagens por vetores para o cálculo explícito de $\bold{X}_{k+1}$ e podemos também eliminar o vetor para $\bold{\Delta}_k$, nos deixando com 9 passagens por vetores e apenas 3 vetores armazenados, metade do que era na versão base.

Além disso, podemos fundir a avaliação do primeiro critério de convergência ao mesmo tempo que montamos o sistema, economizando uma passagem pelo vetor do sistema. Numa ideia similar, podemos fundir a avaliação do segundo critério de convergência com o cálculo de $\bold{X}_k$, economizando outra passagem e totalizando apenas 3 vetores e 7 passagens por vetores, metade do que era necessário na versão base. Note que, portanto, o programa não ficou 2 vezes mais rápido, uma vez que essa métrica de desempenho é uma aproximação, como dito anteriormente. Mesmo assim, o ganho de performance é significativo, como pode se observar no gráfico abaixo: 

### Nota sobre vetorização
Na versão anterior, a construção do sistema linear era facilmente vetorizada pelo compilador, calculando diversas linhas simultaneamente. Com a adição do cálculo da norma, o compilador já não conseguia vetorizar tudo tão facilmente, devido ao cálculo do valor máximo do vetor do sistema. Para ajudar o compilador, foi explicitamente declarado que a norma era uma [redução ao valor máximo](https://bisqwit.iki.fi/story/howto/openmp/#TheReductionClause_2) (do valor absoluto) usando OpenMP:

```cpp
double build_system_compute_norm(TridiagonalSystem* system, double* x) {
	// ...
	#pragma omp simd reduction(max:norma)
	for (i = 1; i < system->n - 1; ++i) {
		system->rhs[i] = (2.0 * x[i] - 3.0) * x[i] + x[i - 1] + 2.0 * x[i + 1] - 1.0;
		system->main[i] = -4.0 * x[i] + 3.0;

		abs = fabs(system->rhs[i]);
		norma = norma > abs ? norma : abs;
	}
	// ...
}
```
Desta forma, o compilador foi capaz de vetorizar a construção do sistema linear (verificado através dos logs de vetorização ao compilar o programa). No entando, a resolução do sistema linear utilizando o algoritmo de Thomas é inerentemente sequencial, impossibilitando a vetorização. 

## Versão 3
Nesta versão, a estrutura organizada do código será completamente descartada em troca da performance. Inicialmente, deixa-se de calcular explicitamente o sistema linear em cada iteração, e cada elemento é calculado na hora em que é preciso, diretamente na triangularização. Isso evita ter que preencher os vetores do sistema linear para logo depois percorrê-los novamente. Com isso, já se torna difícil distinguir entre que tempo foi gasto construindo o sistema linear e que tempo foi gasto resolvendo o mesmo, pois os dois passos estão sendo feitos ao mesmo tempo.

Feito isso, tem pouco o que podemos fazer para acelerar ainda mais a resolução de um sistema linear, pelo menos sem utilizar [múltiplos cores ou uma GPU](https://research.nvidia.com/sites/default/files/pubs/2010-01_Fast-Tridiagonal-Solvers/Zhang_Fast_2009.pdf), que foge do escopo desta atividade. No entanto, é essencial notar que não estamos resolvendo apenas um sistema linear, mas  vários, um após o outro. Este detalhe faz com que ainda tenha espaço para mais otimizações. Perceba que enquanto efetuamos a retro-substituição no algoritmo de Thomas na iteração k, já temos parte do vetor $\bold{X}_k$, mesmo que incompleto. Como o sistema linear da próxima iteração é tridiagonal, não precisamos de $\bold{X}_k$ completo para calcular parcialmente $\text{diag}(\bold{J}_\bold{F}(\bold{X}_k))$ ou $-\bold{F}(\bold{X}_k)$. De fator, para a linha $i$ do sistema, precisamos apenas das linhas $i - 1$, $i$ e $i + 1$. Assim, enquanto efetuamos a retro-substituição na iteração k, podemos simultaneamente calcular o sistema linear para a iteração k + 1. Consequentemente, podemos também triangularizar o sistema da próxima iteração, tudo enquanto ainda estamos fazendo a retro-substituição da iteração atual. Nesse caso, a triangularização transformará a matriz em triangular inferior ao invés de triangular superior, uma vez que estamos triangularizando de cima para baixo, mas isso não é um problema. Dessa forma, quando terminarnos a iteração k, já teremos metade do trabalho para a iteração k + 1 feito, e assim sucessivamente. Caso em alguma iteração verificarmos que o método já convergiu, podemos apenas encerrar o procedimento com medate da próxima iteração feita. Isso representa um pequeno desperdício, mas o ganho em eficiência no resto das iterações deve ser suficiente para compensar por essa pequena perda.

Em uma visão geral, cada iteração agora:
 - Encontra a solução (usando o sistema já tridiagonalizado na iteração anterior)
 - Calcula a norma do deslocamento
 - Monta o sistema da próxima iteração
 - Triangulariza o sistema da próxima iteração

O que resulta no código (com algumas simplificações)
```cpp
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
```

Essa otimização, mais uma vez, ignora completamente uma estrutura organizada do código. Agora, uma iteração começa antes mesmo da última terminar, não há separação entre que parte do programa monta o sistema linear, que parte resolve, que parte verifica a convergência. Será tudo apenas um emaranhado de código que serve ao objetivo único resolver esse problema específico, sendo muito difícil reaproveitar qualquer parte. Este é, portanto, um sacrifício necessário para atingirmos o máximo de performance na resolução do problema proposto nessa atividade.

## Apêndice

### Arquitetura do processador utilizado nos testes
```
--------------------------------------------------------------------------------
CPU name:	Intel(R) Core(TM) i5-10400 CPU @ 2.90GHz
CPU type:	Intel Cometlake processor
CPU stepping:	3
********************************************************************************
Hardware Thread Topology
********************************************************************************
Sockets:		1
CPU dies:		1
Cores per socket:	6
Threads per core:	2
--------------------------------------------------------------------------------
HWThread        Thread        Core        Die        Socket        Available
0               0             0           0          0             *                
1               0             1           0          0             *                
2               0             2           0          0             *                
3               0             3           0          0             *                
4               0             4           0          0             *                
5               0             5           0          0             *                
6               1             0           0          0             *                
7               1             1           0          0             *                
8               1             2           0          0             *                
9               1             3           0          0             *                
10              1             4           0          0             *                
11              1             5           0          0             *                
--------------------------------------------------------------------------------
Socket 0:		( 0 6 1 7 2 8 3 9 4 10 5 11 )
--------------------------------------------------------------------------------
********************************************************************************
Cache Topology
********************************************************************************
Level:			1
Size:			32 kB
Type:			Data cache
Associativity:		8
Number of sets:		64
Cache line size:	64
Cache type:		Non Inclusive
Shared by threads:	2
Cache groups:		( 0 6 ) ( 1 7 ) ( 2 8 ) ( 3 9 ) ( 4 10 ) ( 5 11 )
--------------------------------------------------------------------------------
Level:			2
Size:			256 kB
Type:			Unified cache
Associativity:		4
Number of sets:		1024
Cache line size:	64
Cache type:		Non Inclusive
Shared by threads:	2
Cache groups:		( 0 6 ) ( 1 7 ) ( 2 8 ) ( 3 9 ) ( 4 10 ) ( 5 11 )
--------------------------------------------------------------------------------
Level:			3
Size:			12 MB
Type:			Unified cache
Associativity:		16
Number of sets:		12288
Cache line size:	64
Cache type:		Inclusive
Shared by threads:	12
Cache groups:		( 0 6 1 7 2 8 3 9 4 10 5 11 )
--------------------------------------------------------------------------------
********************************************************************************
NUMA Topology
********************************************************************************
NUMA domains:		1
--------------------------------------------------------------------------------
Domain:			0
Processors:		( 0 6 1 7 2 8 3 9 4 10 5 11 )
Distances:		10
Free memory:		425.742 MB
Total memory:		7764.93 MB
--------------------------------------------------------------------------------


********************************************************************************
Graphical Topology
********************************************************************************
Socket 0:
+-------------------------------------------------------------------+
| +--------+ +--------+ +--------+ +--------+ +--------+ +--------+ |
| |  0 6   | |  1 7   | |  2 8   | |  3 9   | |  4 10  | |  5 11  | |
| +--------+ +--------+ +--------+ +--------+ +--------+ +--------+ |
| +--------+ +--------+ +--------+ +--------+ +--------+ +--------+ |
| |  32 kB | |  32 kB | |  32 kB | |  32 kB | |  32 kB | |  32 kB | |
| +--------+ +--------+ +--------+ +--------+ +--------+ +--------+ |
| +--------+ +--------+ +--------+ +--------+ +--------+ +--------+ |
| | 256 kB | | 256 kB | | 256 kB | | 256 kB | | 256 kB | | 256 kB | |
| +--------+ +--------+ +--------+ +--------+ +--------+ +--------+ |
| +---------------------------------------------------------------+ |
| |                             12 MB                             | |
| +---------------------------------------------------------------+ |
+-------------------------------------------------------------------+
```
