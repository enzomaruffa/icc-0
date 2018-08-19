#include <stdio.h>
#include <math.h>

#include "utils.h"
#include "ZeroFuncao.h"

double calcErroRelativo(double xmold, double xmnew)
{
    return (fabs((xmnew - xmold) / xmnew) * 100);
}

double bisseccao (Polinomio p, double a, double b, double eps,
	       int *it, double *raiz)
{
    /* xl, xu, xm, f(xl), f(xu), f(xm), f'(xl), f'(xu), f'(xm)*/
    int iter = 1;
    double xl, xu, xm, fxl, fxu, fxm, dfxl, dfxu, dfxm, xmold, sinal, erro;
    xl = a;
    xu = b;
    xm = 0;
    erro = 0;

    *it = 0;

    /* Calcula f(xl) e f(xu) */
    calcPolinomio_lento(p, xl, &fxl, &dfxl);
    calcPolinomio_lento(p, xu, &fxu, &dfxu);

    /* Passo 1 */
    if ((fxl * fxu) > 0) {
        printf("Saiu em fxl*fxu > 0, pois fxl valia %f para xl = %f e fxu valia %f para xu = %f\n", fxl, xl, fxu, xu);
        exit(2);
    }

    /* Passo 2 */

    /* Incrementa a iteração */
    (*it) += 1;

    /* Estimando xm como ponto médio */
    xmold = xm;
    xm = (xl + xu) / 2;

    /* Passo 3 */
    while (*it < MAXIT) {
        /* Calcula f(xm) */
        calcPolinomio_lento(p, xm, &fxm, &dfxm);

        sinal = fxl * fxm;

        if (sinal < 0) { /* Ou seja, a raiz está entre xl e xm */
            /* xl = xl */

            /* Atribui valores pro xu para não recalcular tudo */
            xu = xm;
            fxu = fxm;
        } else if (sinal > 0) { /* Ou seja, a raiz está entre xm e xu */
            /* xu = xu */

            /* Atribui valores pro xl para não recalcular tudo */
            xl = xm;
            fxl = fxm;
        } else { /* Ou seja, encontramos a raiz */
            *raiz = xm;
            return erro;
        }

        /* Passo 4 - achar novo xm */

        /* Estimando xm como ponto médio */
        xmold = xm;
        xm = (xl + xu) / 2;

        /* Passo 5 */
        erro = calcErroRelativo(xmold, xm);

        if (erro <= eps) {   
            *raiz = xm;
            return erro;
        }

        printf("it: %d, xm: %f, erro: %f \n", *it, xm, erro);
        (*it) += 1;
    }

    exit(2);
    return erro;
}


double newtonRaphson (Polinomio p, double x0, double eps,
		   int *it, double *raiz)
{
    return 0;
}


double secante (Polinomio p, double x0, double x1, double eps,
	     int *it, double *raiz)
{
    return 0;
}

void calcPolinomio_rapido(Polinomio p, double x, double *px, double *dpx)
{
    *px = 0;
    *dpx = 0;
    double x_mult = 1;

    /* Primeira iteração exclusiva do px onde x tem grau 0. Na derivada esse termo não existe */
    *px += p.p[0] /* * x_mult */ ;

    /* A partir do grau 1, os termos são "iguais". O que muda é que na mesma iteração, px deve ser calculado normalmente enquando dpx deve ter o x elevado até 1 grau a menos e o grau da iteração atual deve multiplicar o coeficiente */
    for (int i = 1; i <= p.grau; i++) {
        *dpx += x_mult * (p.p[i]*i));

        x_mult *= x;

        *px += x_mult * p.p[i]);
    }
}

void calcPolinomio_lento(Polinomio p, double x, double *px, double *dpx)
{
    *px = 0;
    *dpx = 0;

    /* Primeira iteração exclusiva do px onde x tem grau 0. Na derivada esse termo não existe */
    *px += pow(x, 0) * p.p[0];

    /* A partir do grau 1, os termos são "iguais". O que muda é que na mesma iteração, px deve ser calculado normalmente enquando dpx deve ter o x elevado até 1 grau a menos e o grau da iteração atual deve multiplicar o coeficiente */
    for (int i = 1; i <= p.grau; i++) {
        *px += (pow(x, i) * p.p[i]);
        *dpx += (pow(x, i-1) * (p.p[i]*i));
    }
}