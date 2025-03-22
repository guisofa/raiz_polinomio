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
	pol.p = malloc(sizeof(real_t) * pol.grau);

  	for (int i=pol.grau; i >=0; --i)
    	scanf("%lf", &pol.p[i]);

  	scanf("%lf %lf", &a, &b); // intervalo onde está uma das raizes.


	real_t v = 0, dv = 0;
	calcPolinomio_rapido(pol, 2, &v, &dv);
	printf("v = %lf\n", v);
	printf("dv = %lf\n", dv);

  return 0;
}

