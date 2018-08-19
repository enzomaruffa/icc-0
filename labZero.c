#include <stdio.h>
#include <math.h>

#include "utils.h"
#include "ZeroFuncao.h"

/* Ponteiro de função utilizado pra trocar o método de cálculo do polinômio */
extern void (*calcPolinomio_funct)(Polinomio, double, double *, double *);

int main ()
{
    Polinomio p;
    p.grau = 10;
    p.p = malloc(sizeof(double) * (p.grau));

    /* Criar polinômio com i = grau do x */
    p.p[0] = 3;
    p.p[1] = -1;
    p.p[2] = -5;
    p.p[3] = 2;
    p.p[4] = 1;
    p.p[5] = -2;
    p.p[6] = -5;
    p.p[7] = 5;
    p.p[8] = 4;
    p.p[9] = 2;

    int it;
    double raiz, erro;

    double a, b, x0, x1, aux;

    printf("Digite o a: ");
    scanf("%lf", &a);

    printf("Digite o b: ");
    scanf("%lf", &b);

    printf("Digite o x0: ");
    scanf("%lf", &x0);

    printf("Digite o x1: ");
    scanf("%lf", &x1);

    /* Troca de segurança */
    if (a > b) {
        aux = a;
        a = b;
        b = aux;
    }

    double tempoBisseccao, tempoNewtonRaphson, tempoSecante, tempoAux;

    for (int i = 0; i < 2; i++) {
        if (i == 0)
            printf("\nUsando o método de calculo do polinômio rápido\n");
        else
            printf("\nUsando o método de calculo do polinômio lento\n");

        it = 0;
        raiz = 0;
        erro = -1;

        printf("--------------------------------------------------------------------\n");
        printf("|      metodo      |    raiz   |   erro    | iteracoes |   tempo   |\n");
        printf("--------------------------------------------------------------------\n");

        tempoAux = timestamp();
        erro = bisseccao(p, a, b, EPS, &it, &raiz);
        tempoBisseccao = timestamp() - tempoAux;
        printf("|        bisseccao | %9f | %9f | %9d | %9f |\n", raiz, erro, it, tempoBisseccao);
        printf("--------------------------------------------------------------------\n");

        it = 0;
        raiz = 0;
        erro = -1;

        tempoAux = timestamp();
        erro = newtonRaphson(p, x0, EPS, &it, &raiz);
        tempoNewtonRaphson = timestamp() - tempoAux;
        printf("|   newton-raphson | %9f | %9f | %9d | %9f |\n", raiz, erro, it, tempoNewtonRaphson);
        printf("--------------------------------------------------------------------\n");

        it = 0;
        raiz = 0;
        erro = -1;

        tempoAux = timestamp();
        erro = secante(p, x0, x1, EPS, &it, &raiz);
        tempoSecante = timestamp() - tempoAux;
        printf("|          secante | %9f | %9f | %9d | %9f |\n", raiz, erro, it, tempoSecante);
        printf("--------------------------------------------------------------------\n");

        /* ALTERANDO PARA POLINOMIO LENTO*/
        calcPolinomio_funct = calcPolinomio_lento;
    }

    return 0;
}

