#include <stdio.h>
#include <math.h>

#include "utils.h"
#include "ZeroFuncao.h"

int main ()
{
    Polinomio p;
    p.grau = 5;
    p.p = malloc(sizeof(double) * (p.grau + 1));

    /* Criar polinômio com i = grau do x */
    p.p[0] = -15;
    p.p[1] = -10;
    p.p[2] = 3;
    p.p[3] = 4;
    p.p[4] = 2;

    double a = 0;
    double b = 3;

    int it = 0;
    double raiz = 0;
    double erro = -1;

    erro = bisseccao(p, a, b, EPS, &it, &raiz);

    printf("it: %d, raiz: %f, erro: %f \n", it, raiz, erro);

    return 0;
}

