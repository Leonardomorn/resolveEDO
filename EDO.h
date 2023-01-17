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

/*********************************
 * estrutura que guarda as informações da equação diferencial ordinária específica com valores de contorno
 * n: número de pontos internos na malha
 * a: início do intervalo
 * b: final do intervalo
 * ya: condição de contorno inicial
 * yb: condição de contorno final
 * p: ponteiro de função que calcula o coeficiente que multiplica y'
 * q: coeficiente que multiplica y
 * r: ponteiro de função que calcula o coeficiente independente
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

double somaKahan( double *dados, int tam );
double norma_L2(double *x, double *xAnt, int n);
void aloca_tri_diagonal (Edo *edoeq, SL_Tridiag *sl);
void gera_tri_diagonal (Edo *edoeq, SL_Tridiag *sl);
void gaussSeidel (double *D, double *Di, double *Ds, 
double *B, double *x, int n, double norma );
