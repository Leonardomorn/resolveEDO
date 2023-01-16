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

/*****************************************************
 * Função que gera a matriz Tri-diagonal 
 * edoeq: equação diferencial ordinária
 * sl : sistema linear tri-diagonal
******************************************************/

void geraTriDiagonal (Edo *edoeq, SL_Tridiag *sl)
{
double xi, h;
h = (edoeq->b - edoeq->a) / (edoeq->n+1.0);
for (int i=0; i < edoeq->n; ++i) {
xi = edoeq->a + (i+1)*h; // ponto da malha
sl->Di[i] = 1 - h * edoeq->p(xi)/2.0; // diagonal inferior
sl->D[i] = -2 + h*h * edoeq->q(xi); // diagonal principal
sl->Ds[i] = 1 + h * edoeq->p(xi)/2.0; // diagonal superior
sl->B[i] = h*h * edoeq->r(xi); // termo independente
}
// Condições de contorno subtraídas do 1º e último termos independentes
sl->B[0] -= edoeq->ya * (1 - h*edoeq->p(edoeq->a+h)/2.0);
sl->B[edoeq->n-1] -= edoeq->yb * (1 + h*edoeq->p(edoeq->b-h)/2.0);
}

void gaussSeidel (double *D, double *Di, double *Ds,
double *B, double *x, int n, double tol)
{
    double erro = 1.0 + tol;
    while (erro > tol) {
    // 5(n-2)+6 ≈ 5n operações / iteração
    x[ 0 ] = (B[ 0 ] – Ds[ 0 ] * x[ 1 ]) / d[ 0 ];
    for (int i=1; i < n-1; ++i)
    x[ i ] = (B[ i ] – a[ i-1 ] * x[ i-1] – Ds[ i ] * x[ i+1 ]) / d[ i ];
    x[ n-1 ] = (B[ n-1 ] – a[ n-2 ] * x[ n-2 ] ) / d[ n-1 ];
    // Calcula erro
    }
}