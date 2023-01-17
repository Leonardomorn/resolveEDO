#include "EDO.h"
#define MAXIT 50
/**********************************
 * função que calcula e retorna o coeficiente p da EDO, onde p = x + 1;
 * x: valor de x no passo n 
************************************/
inline double calcula_p (double x)
{
    double p = x + 1;
    return p;
}
/**********************************
 * função que calcula e retorna o coeficiente q da EDO, onde q = -2;
 * x: valor de x no passo n 
************************************/
inline double calcula_q (double x)
{
    double q = -2.0;
    return q;
}

/**********************************
 * função que calcula o coeficiente r da EDO, onde r = (1−x²) * e^(−x) 
 * x: valor de x no passo n 
************************************/
inline double calcula_r(double x)
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

/**********************************
 * função que ajusta o número de passos do método das diferenças finitas 
 * edoq: equação diferencial ordinária
 * n: numero de passos 
************************************/
void ajusta_n_EDO( Edo* edoq, int n)
{
    edoq->n = n;
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

/*!
    \brief algoritmo que soma os numeros de um vetor de modo a evitar
        muitos cancelamentos subtrativos
    \param dados vetor
    \param tam tamanho do vetor
*/
double somaKahan( double *dados, int tam )
{
    double soma = 0.0; // Prepara o acumulador
    double compensador = 0.0;   // compensador para a perda de bits de baixa ordem

    for(int i = 0; i < tam; i++)
    {
        double y = dados[i] - compensador;
        double t = soma + y;
        compensador = (t - soma) - y;
        soma = t;
    }
    return soma;
}

/*!
  \brief Essa função calcula a norma L2 das etapas do Gauss-Seidel 

  \param x Vetor solução do sistema linear
  \param xAnt Vetor solução do sistema linear da iteração anterior
  \param n Tamanho dos vetores
*/
double norma_L2(double *x, double *xAnt, int n)
{
 int tam = n;
 double *auxL2 = malloc (n * sizeof(double));
 double normaL2 = 0.0;

 for (int i = 0; i < tam; i++)
 {
  auxL2[i] = x[i] - xAnt[i];
  auxL2[i] = auxL2[i] * auxL2[i];
 }
 
 normaL2 = ABS(somaKahan(auxL2, tam));
 normaL2 = sqrt(normaL2);
 free(auxL2);
 return normaL2; 
}
/*!
    \brief aloca espaço para os vetores a partir do tamanho n determinado na EDO
    
    \param edoeq equação diferencial ordinária
    \param sl sistema linear tri-diagonal 
*/
void aloca_tri_diagonal (Edo *edoeq, SL_Tridiag *sl)
{
    sl->B = malloc ( edoeq->n * sizeof(double));
    sl->D = malloc ( edoeq->n * sizeof(double));
    sl->Di = malloc ( edoeq->n * sizeof(double));
    sl->Ds = malloc ( edoeq->n * sizeof(double));
}

/*!
  \brief algoritmo que gera a matriz Tri-diagonal 

  \param edoeq equação diferencial ordinária
  \param  sl sistema linear tri-diagonal
*/
void gera_tri_diagonal (Edo *edoeq, SL_Tridiag *sl)
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

/*!
    \brief aplica o método de Gauss Seidel para obter a solução de um sistema tri-driagonal
    
    \param D vetor diagonal principal
    \param Di vetor diagonal inferior
    \param Ds vetor diagonal superior
    \param B vetor de termos independentes
    \param x vetor de soluções
    \param n tamanho do sistema
    \param norma retorna última normaL2 do método
    */
void gaussSeidel (double *D, double *Di, double *Ds, 
double *B, double *x, int n, double norma )
{
    double *xAnt = calloc (n , sizeof(double));  //vetor de solucao anterior, começa nulo
    int it = 0;
    norma = 1.0 + FLT_EPSILON;

    //aplica método GaussSeidel
    while (norma > FLT_EPSILON && it < MAXIT) {
    x[0] = (B[0] - Ds[0] * x[1]) / D[0];
    for (int i=1; i < n-1; ++i)
        x[i] = (B[i] - Di[i-1] * x[i-1] - Ds[i] * x[i+1]) / D[i];
    x[n-1] = (B[n-1] - Di[n-2] * x[n-2] ) / D[n-1];
    norma = norma_L2(x, xAnt, n);

    //copia vetor de solução para o vetor de solução anterior
    for (int iterator = 0; iterator < n; iterator++) 
        xAnt[iterator] = x[iterator];
    it++;
    }

    free(xAnt);
}

