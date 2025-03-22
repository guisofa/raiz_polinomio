#include <stdio.h>
#include <math.h>

#include "utils.h"
#include "ZeroFuncao.h"

// Retorna valor do erro quando método finalizou. Este valor depende de tipoErro
real_t newtonRaphson (Polinomio p, real_t x0, int criterioParada, int *it, real_t *raiz)
{

}


// Retorna valor do erro quando método finalizou. Este valor depende de tipoErro
real_t bisseccao (Polinomio p, real_t a, real_t b, int criterioParada, int *it, real_t *raiz)
{

}


void calcPolinomio_rapido(Polinomio p, real_t x, real_t *px, real_t *dpx)
{
    real_t res_fun = 0, res_der = 0;

    for (int i = p.grau; i > 0; i--) {
        res_fun = res_fun * x + p.p[i];
        res_der = res_der * x + p.p[i] * i;
    }
    res_fun = res_fun * x + p.p[0];

    *px = res_fun;
    *dpx = res_der;
}


void calcPolinomio_lento(Polinomio p, real_t x, real_t *px, real_t *dpx)
{
    real_t res_fun = 0, res_der = 0;

    for (int i = p.grau; i > 0; i--) {
        res_fun += p.p[i] * pow(x, i);
        res_der += p.p[i] * i * pow(x, i-1);
    }
    res_fun += p.p[0];

    *px = res_fun;
    *dpx = res_der;
}
