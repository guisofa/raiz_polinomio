#include <stdio.h>
#include <math.h>

#include "utils.h"
#include "ZeroFuncao.h"

// Retorna valor do erro quando método finalizou. Este valor depende de tipoErro
real_t newtonRaphson (Polinomio p, real_t x0, int criterioParada, int *it, real_t *raiz,
                      void calcPolinomio(Polinomio, real_t, real_t*, real_t*))
{
    real_t valorParada;
    if (criterioParada == 1) valorParada =  EPS;
    else if (criterioParada == 2) valorParada = ZERO;
    else valorParada = ULPS;
    /* 1 = erro  2 = f(x)  3 = ulps*/
    real_t valoresCriticos[4];

    *it = 0;
    Double_t x; x.f = x0;
    Double_t x_next; x_next.f = x0;
    real_t func_val = 0, der_val = 0;

    calcPolinomio(p, x.f, &func_val, &der_val);
    do {
        x_next.f = x.f - func_val/der_val;

        (*it)++;
        valoresCriticos[1] = fabs(x_next.f - x.f);
        valoresCriticos[3] = labs(x.i - x_next.i) -1;
        x.f = x_next.f;
        calcPolinomio(p, x.f, &func_val, &der_val);
        valoresCriticos[2] = fabs(func_val);
    } while (*it < 500 && valoresCriticos[criterioParada] > valorParada);

    *raiz = x.f;
    if (criterioParada == 3) return (valoresCriticos[3] < 0) ? 0 : valoresCriticos[3];
    else return fabs(valoresCriticos[criterioParada]);
}


// Retorna valor do erro quando método finalizou. Este valor depende de tipoErro
real_t bisseccao (Polinomio p, real_t a, real_t b, int criterioParada, int *it, real_t *raiz,
                  void calcPolinomio(Polinomio, real_t, real_t*, real_t*))
{
    real_t valorParada;
    if (criterioParada == 1) valorParada =  EPS;
    else if (criterioParada == 2) valorParada = ZERO;
    else valorParada = ULPS;
    /* 1 = erro  2 = f(x)  3 = ulps*/
    real_t valoresCriticos[4];

    *it = 0;
    Double_t xm, x_next;
    real_t func_val_a = 0, func_val_m = 0;

    xm.f = a/2 + b/2;
    calcPolinomio(p, xm.f, &func_val_m, NULL);
    do {
        calcPolinomio(p, a, &func_val_a, NULL);

        if (func_val_a*func_val_m <= 0)
            b = xm.f;
        else
            a = xm.f;

        x_next.f = a/2 + b/2;
        (*it)++;
        valoresCriticos[1] = fabs(x_next.f - xm.f);
        valoresCriticos[3] = labs(xm.i - x_next.i) -1;
        xm.f = x_next.f;
        calcPolinomio(p, xm.f, &func_val_m, NULL);
        valoresCriticos[2] = fabs(func_val_m);
        /*printf("%d: %.16e, %.16e\n", *it, xm.f ,func_val_m);*/
    } while (*it < 500 && valoresCriticos[criterioParada] > valorParada);

    *raiz = xm.f;
    
    if (criterioParada == 3) return (valoresCriticos[3] < 0) ? 0 : valoresCriticos[3];
    else return fabs(valoresCriticos[criterioParada]);
}

void calcPolinomio_rapido(Polinomio p, real_t x, real_t *px, real_t *dpx)
{
    real_t res_fun = 0, res_der = 0;

    for (int i = p.grau; i > 0; i--) {
        res_fun = res_fun * x + p.p[i];
        res_der = res_der * x + res_fun;
    }
    res_fun = res_fun * x + p.p[0];

    *px = res_fun;
    if (dpx)
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
    if (dpx)
        *dpx = res_der;
}
