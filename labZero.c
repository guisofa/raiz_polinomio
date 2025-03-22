#include <stdio.h>
#include <math.h>
#include <float.h>

#include "utils.h"
#include "ZeroFuncao.h"

int main ()
{

	real_t a, b;
  	Polinomio pol;

  	scanf("%d", &pol.grau);
	pol.p = malloc(sizeof(real_t) * (pol.grau+1));

  	for (int i=pol.grau; i >=0; --i)
    	scanf("%lf", &pol.p[i]);

  	scanf("%lf %lf", &a, &b); // intervalo onde está uma das raizes.

	real_t raiz;
	int iteracoes = 0;
	real_t erro = newtonRaphson(pol, (a + b)/2, 0, &iteracoes, &raiz, calcPolinomio_rapido);

	printf("raiz = %.15e\niteracoes = %d\nerro = %.15e\n", raiz, iteracoes, erro);
  return 0;
}

