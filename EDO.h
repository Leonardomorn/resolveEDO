#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

// Valor absoluto de um número. Alternativa ao uso da função 'fabs()'
#define ABS(num)  ((num) < 0.0 ? -(num) : (num))

/***********************************
 * estrutura que armazena os vetores do Sistema Linear Tridiagonal
 * D: vetor diagonal principal
 * Di: vetor diagonal inferior
 * Ds: vetor diagonal superior
 * B: vetor termos independentes
 * n : tamanho dos vetores
 ******************************/
typedef struct {
double *D, *Di, *Ds, *B;
int n;
} SL_Tridiag;

/*!
 \param estrutura que guarda as informações da equação diferencial ordinária específica com valores de contorno
 \param n: número de pontos internos na malha
 \param a: início do intervalo
 \param b: final do intervalo
 \param ya: condição de contorno inicial
 \param yb: condição de contorno final
 \param p: ponteiro de função que calcula o coeficiente que multiplica y'
 \param q: coeficiente que multiplica y
 \param r: ponteiro de função que calcula o coeficiente independente
*/
typedef struct {
int n;
double a, b;
double ya, yb;
double (* p)(double), (* q)(double), (* r)(double);
} Edo;

double calcula_p (double x);
double calcula_q (double x);
double calcula_r(double x);

void inicializa_EDO( Edo* edoq, int n);
void ajusta_n_EDO( Edo* edoq, int n);
void imprime_EDO (Edo* edoq);
void imprime_SL (Edo* edoq, SL_Tridiag* sl);

double somaKahan( double *dados, int tam );
double norma_L2_residuo(double *x, SL_Tridiag *sl, int n);
void aloca_tri_diagonal (Edo *edoeq, SL_Tridiag *sl);
void gera_tri_diagonal (Edo *edoeq, SL_Tridiag *sl);
void gaussSeidel (double *D, double *Di, double *Ds, 
double *B, double *x, int n, double *norma, SL_Tridiag *sl );
