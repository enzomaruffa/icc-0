#include <stdio.h>
#include <math.h>

#include "utils.h"
#include "ZeroFuncao.h"

/* Ponteiro de função utilizado para determinar o tipo de cálculo de polinômio que deve ser utilizado */
void (*calcPolinomio_funct)(Polinomio, double, double *, double *) = calcPolinomio_rapido;

double calcErroRelativo(double xmold, double xmnew)
{
    return (fabs((xmnew - xmold) / xmnew) * 100);
}

double bisseccao (Polinomio p, double a, double b, double eps,
	       int *it, double *raiz)
{
    /* xl, xu, xm, f(xl), f(xu), f(xm), f'(xl), f'(xu), f'(xm)*/
    double xl, xu, xm, fxl, fxu, fxm, dfxl, dfxu, dfxm, xmold, sinal, erro;
    xl = a;
    xu = b;
    xm = 0;
    erro = 0;

    *it = 1;

    /* Calcula f(xl) e f(xu) */
    calcPolinomio_funct(p, xl, &fxl, &dfxl);
    calcPolinomio_funct(p, xu, &fxu, &dfxu);

    /* Passo 1 */
    if ((fxl * fxu) > 0) {
        printf("[BISSECCAO] Saiu em fxl*fxu > 0, pois fxl valia %f para xl = %f e fxu valia %f para xu = %f\n", fxl, xl, fxu, xu);
        return -1;
    }

    /* Passo 2 */

    /* Estimando xm como ponto médio */
    xmold = xm;
    xm = (xl + xu) / 2;

    /* Passo 3 */
    while (*it < MAXIT) {
        /* Calcula f(xm) */
        calcPolinomio_funct(p, xm, &fxm, &dfxm);

        sinal = fxl * fxm;
        /* printf("%f * %f tem sinal %f\n", fxl, fxm, sinal); */

        if (sinal < ZERO) { /* Ou seja, a raiz está entre xl e xm */
            /* xl = xl */

            /* Atribui valores pro xu para não recalcular tudo */
            xu = xm;
            fxu = fxm;
        } else if (sinal > ZERO) { /* Ou seja, a raiz está entre xm e xu */
            /* xu = xu */   

            /* Atribui valores pro xl para não recalcular tudo */
            xl = xm;
            fxl = fxm;
        } else { /* Ou seja, encontramos a raiz */
            *raiz = xm;
            /* printf("Saiu porque encontrou a raiz\n"); */
            return erro;
        }

        /* Passo 4 - achar novo xm */

        /* Estimando xm como ponto médio */
        xmold = xm;
        xm = (xl + xu) / 2;

        /* Passo 5 */
        erro = calcErroRelativo(xmold, xm);

        if (erro <= eps) {
            /* printf("Saiu porque erro é menor\n"); */
            *raiz = xm;
            return erro;
        }

        /* printf("it: %d, xm: %f, erro: %f \n\n", *it, xm, erro); */
        (*it) += 1;
    }

    return erro;
}


double newtonRaphson (Polinomio p, double x0, double eps,
		   int *it, double *raiz)
{
    /* xn denota o x que será utilizado nas iterações, enquanto x0 será o x "old" */
    double xn, fxn, dfxn, erro;

    *it = 1;
    erro = 0;

    /* Passo 1 */
    xn = x0;

    /* Passo 2 */
    while (*it < MAXIT) {
        calcPolinomio_funct(p, xn, &fxn, &dfxn);

        /* Calcula xi+1 */
        x0 = xn;
        xn = xn - (fxn / dfxn);

        /* Passo 3 */

        /* Calcula o erro */
        erro = calcErroRelativo(x0, xn);

        /* Passo 4 */
        if (erro <= eps) {
            /* printf("Saiu porque erro é menor\n"); */
            *raiz = xn;
            return erro;
        }

        /* printf("it: %d, xm: %f, erro: %f \n\n", *it, xm, erro); */
        (*it) += 1;
    }

    return 0;
}


double secante (Polinomio p, double x0, double x1, double eps,
	     int *it, double *raiz)
{
    double fx0, dfx0, fx1, dfx1, aux, erro;

    *it = 1;
    erro = 0;

    calcPolinomio_funct(p, x0, &fx0, &dfx0);

    /* Passo 1 */
    while (*it < MAXIT) {
        /* Calculando a derivada de xi */
        calcPolinomio_funct(p, x1, &fx1, &dfx1);

        dfx1 = (fx1 - fx0) / (x1 - x0);

        /* Passo 2 */

        /* Calcula o erro */
        erro = calcErroRelativo(x0, x1);

        /* Passo 3 */
        if (erro <= eps) {
            /* printf("Saiu porque erro é menor\n"); */
            *raiz = x1;
            return erro;
        }

        /* Passo 4 */

        /* Calcula xi+1 e coloca no aux */
        aux = x1 - (fx1*(x1-x0) / (fx1-fx0));

        /* printf("x0: %f  /  x1: %f", x0, x1); */

        x0 = x1;
        fx0 = fx1;
        dfx0 = dfx1;

        x1 = aux;

        /* printf("it: %d, x1: %f, erro: %f \n\n", *it, x1, erro); */
        (*it) += 1;
    }

    return 0;

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
        *dpx += x_mult * (p.p[i]*i);

        x_mult *= x;

        *px += x_mult * p.p[i];
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