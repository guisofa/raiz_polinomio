#include <stdio.h>
#include <math.h>

#include "utils.h"
#include "ZeroFuncao.h"

// Retorna valor do erro quando método finalizou. Este valor depende de tipoErro
real_t newtonRaphson (Polinomio p, real_t x0, int criterioParada, int *it, real_t *raiz,
                      void calcPolinomio(Polinomio, real_t, real_t*, real_t*))
{
    *it = 0;
    Double_t x; x.f = x0;
    Double_t x_next; x_next.f = x0;
    real_t func_val = 0, der_val = 0;
    real_t erro = 0;
    long int ulps = 0;


    calcPolinomio(p, x.f, &func_val, &der_val);
    do {
        x_next.f = x.f - func_val/der_val;

        (*it)++;
        erro = fabs(x_next.f - x.f);
        ulps = labs(x.i - x_next.i) -1;
        x.f = x_next.f;
        calcPolinomio(p, x.f, &func_val, &der_val);
    } while (*it < 500 && (criterioParada == 1 ? (erro > EPS) : 1) && (criterioParada == 3 ? (ulps > ULPS) : 1) && (criterioParada == 2 ? (fabs(func_val) > ZERO) : 1));

    *raiz = x.f;
    
    if (criterioParada == 1) return erro;
    else if (criterioParada == 2) return func_val;
    else return ulps;
}


// Retorna valor do erro quando método finalizou. Este valor depende de tipoErro
real_t bisseccao (Polinomio p, real_t a, real_t b, int criterioParada, int *it, real_t *raiz,
                  void calcPolinomio(Polinomio, real_t, real_t*, real_t*))
{
    *it = 0;
    Double_t xm, x_next;
    real_t func_val_a = 0, func_val_m = 0;
    real_t erro = 0;
    long int ulps = 0;

    xm.f = (a + b)/2;
    calcPolinomio(p, xm.f, &func_val_m, NULL);
    do {
        calcPolinomio(p, a, &func_val_a, NULL);

        if (func_val_a*func_val_m <= 0)
            b = xm.f;
        else
            a = xm.f;

        x_next.f = (a + b)/2;
        (*it)++;
        erro = fabs(x_next.f - xm.f);
        ulps = labs(xm.i - x_next.i) -1;
        xm.f = x_next.f;
        calcPolinomio(p, xm.f, &func_val_m, NULL);
    } while (*it < 500 && (criterioParada == 1 ? (erro > EPS) : 1) && (criterioParada == 3 ? (ulps > ULPS) : 1) && (criterioParada == 2 ? (fabs(func_val_m) > ZERO) : 1));

    *raiz = xm.f;
    
    if (criterioParada == 1) return erro;
    else if (criterioParada == 2) return func_val_m;
    else return ulps;
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
