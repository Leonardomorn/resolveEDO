#include "EDO.h"

/**********************************
 * função que calcula e retorna o coeficiente p da EDO, onde p = x + 1;
 * x: valor de x no passo n 
************************************/
double calcula_p (double x)
{
    double p = x + 1;
    return p;
}
/**********************************
 * função que calcula e retorna o coeficiente q da EDO, onde q = -2;
 * x: valor de x no passo n 
************************************/
double calcula_q (double x)
{
    double q = -2.0;
    return q;
}

/**********************************
 * função que calcula o coeficiente r da EDO, onde r = (1−x²) * e^(−x) 
 * x: valor de x no passo n 
************************************/
double calcula_r(double x)
{
    double r;
    r = x * x;
    r = 1 - r;
    r = r * exp(-x);
    return r;
}
/**********************************
 * função que inicializa a equação diferencial ordinária  
 * edoq: equação diferencial ordinária
 * n: numero de passos 
************************************/
void inicializa_EDO( Edo* edoq, int n)
{
    edoq->n = n;
    edoq->a = 0;
    edoq->b = 1;
    edoq->ya = -1;
    edoq->yb = 0;
    edoq->p = calcula_p; //p aponta para calcula_p
    edoq->q = calcula_q; //q aponta para calcula_q
    edoq->r = calcula_r; //r aponta para calcula_r
}

/*******************************************************
 * função que imprime a EDO para teste
 * edoq: equação diferencial ordinária
********************************************************/
void imprime_EDO (Edo* edoq)
{
    printf("\n o número de passos n para esta EDO é: %d", edoq->n);
    printf("\n o intervalo inicial a para esta EDO é: %f", edoq->a);
    printf("\n o intervalo final b para esta EDO é: %f", edoq->b);
    printf("\n a condição de contorno inicial ya para esta EDO é: %f", edoq->ya);
    printf("\n a condição de contorno final yb para esta EDO é: %f\n", edoq->yb);

}