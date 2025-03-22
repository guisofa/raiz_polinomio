#include <stdio.h>
#include <math.h>
#include <float.h>
#include <fenv.h>

#include "utils.h"
#include "ZeroFuncao.h"

int main ()
{
	fesetround(FE_DOWNWARD);

	real_t a, b;
  	Polinomio pol;

  	scanf("%d", &pol.grau);
	pol.p = malloc(sizeof(real_t) * (pol.grau+1));

  	for (int i=pol.grau; i >=0; --i)
    	scanf("%lf", &pol.p[i]);

  	scanf("%lf %lf", &a, &b); // intervalo onde está uma das raizes.

	real_t raiz;
	int iteracoes = 0;
	double tempo;

	printf("RAPIDO\n\n");
	for (int i = 1; i <= 3; i++) {
		tempo = timestamp();
		real_t valor_crit = bisseccao(pol, a, b, i, &iteracoes, &raiz, calcPolinomio_rapido);
		tempo = timestamp() - tempo;
		printf("bissec %.15e %.15e %d %.8e\n", raiz, valor_crit, iteracoes, tempo);
	}
	for (int i = 1; i <= 3; i++) {
		tempo = timestamp();
		real_t valor_crit = newtonRaphson(pol, (a + b)/2, i, &iteracoes, &raiz, calcPolinomio_rapido);
		tempo = timestamp() - tempo;
		printf("newton %.15e %.15e %d %.8e\n", raiz, valor_crit, iteracoes, tempo);
	}
	printf("\nLENTO\n\n");
	for (int i = 1; i <= 3; i++) {
		tempo = timestamp();
		real_t valor_crit = bisseccao(pol, a, b, i, &iteracoes, &raiz, calcPolinomio_lento);
		tempo = timestamp() - tempo;
		printf("bissec %.15e %.15e %d %.8e\n", raiz, valor_crit, iteracoes, tempo);
	}
	for (int i = 1; i <= 3; i++) {
		tempo = timestamp();
		real_t valor_crit = newtonRaphson(pol, (a + b)/2, i, &iteracoes, &raiz, calcPolinomio_lento);
		tempo = timestamp() - tempo;
		printf("newton %.15e %.15e %d %.8e\n", raiz, valor_crit, iteracoes, tempo);
	}

  return 0;
}

