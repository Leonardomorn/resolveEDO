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
    edoq-> n = n;
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

void imprime_SL (Edo* edoq, SL_Tridiag* sl)
{
    printf("\n__________________________________________________________");
    printf("\npara n = %d passos, o sistema linear é:", edoq->n);
    printf("\n\nDiagonal inferior: \n");
    for (int i = 1; i < edoq->n; i++)
        printf("%f " ,sl->Di[i]);
    
    printf("\n\nDiagonal principal: \n");
    for (int i = 0; i < edoq->n; i++)
        printf("%f " ,sl->D[i]);
    
    printf("\n\nDiagonal superior : \n");
    for (int i = 0; i < edoq->n-1; i++)
        printf("%f " ,sl->Ds[i]);    
    printf("\n\n o termo independente é: \n");
    for (int i = 0; i < edoq->n; i++)
        printf("%f ", sl->B[i]);
    
    printf("\n ");    

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
  \brief Essa função calcula a norma L2 do residuo do Gauss-Seidel 

  \param x Vetor solução do sistema linear
  \param xAnt Vetor solução do sistema linear da iteração anterior
  \param n Tamanho dos vetores
*/
double norma_L2_residuo(double *x,SL_Tridiag *sl, int n)
{
 int tam = n;
 double normaL2 = 0.0;
 double* residuo = malloc (n * sizeof(double));
 double* auxResiduo = malloc (n * sizeof(double));


 residuo[0] = sl->B[0] - (sl->D[0] * x[0] + sl->Ds[0] * x[1]);
 for (int i = 1; i < n-1; i++)
    residuo[i] = sl->B[i] - ( (sl->Di[i] * x[i-1]) + (sl->D[i] * x[i]) + (sl->Ds[i] * x[i+1]));
 residuo[n-1] = sl->B[n-1] - (sl->Di[n-1] * x[n-2] + sl->D[n-1] * x [n-1]); 
 

 for (int i = 0; i < tam; i++)
 {
  auxResiduo[i] = x[i] - residuo[i];
  auxResiduo[i] = auxResiduo[i] * auxResiduo[i];
 }
 
 normaL2 = ABS(somaKahan(auxResiduo, tam));
 normaL2 = sqrt(normaL2);
 free(auxResiduo);
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

void libera_tri_diagonal (Edo *edoeq, SL_Tridiag *sl)
{
    free (sl->B); free (sl->D); free (sl->Di); free (sl->Ds); 
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
    \brief aplica o método de Gauss Seidel com vetores para obter a solução de um sistema tri-driagonal
    
    \param D vetor diagonal principal
    \param Di vetor diagonal inferior
    \param Ds vetor diagonal superior
    \param B vetor de termos independentes
    \param x vetor de soluções
    \param n tamanho do sistema
    \param norma retorna última normaL2 do método
    */
void gaussSeidel_vetor (double *D, double *Di, double *Ds, 
double *B, double *x, int n, double*norma, SL_Tridiag* sl, double *tempo )
{
    int it = 0;
    *norma = 1.0 + FLT_EPSILON;
    zera_vetor(x, n);
    //aplica método GaussSeidel
    *tempo = timestamp();
    while (*norma > FLT_EPSILON && it < MAXIT) {
    x[0] = (B[0] - Ds[0] * x[1]) / D[0];
    for (int i=1; i < n-1; ++i)
        x[i] = (B[i] - Di[i-1] * x[i-1] - Ds[i] * x[i+1]) / D[i];
    x[n-1] = (B[n-1] - Di[n-2] * x[n-2] ) / D[n-1];
    *norma = norma_L2_residuo(x, sl, n);
    it++;
    }
    *tempo = timestamp() - *tempo;

}
/*!
    \brief Aplica o método de gauss seidel sem utilizar vetores 
*/
void gaussSeidel_direto(Edo *edoeq, SL_Tridiag *sl, double *X, double *norma, double *tempo)
{
    
    int n = edoeq->n, k = 0 , i = 0;

    double xi, h, yi;
    double Ds, D, Di, B; // diagonais temporárias e termo independente temporário 
    h = (edoeq->b - edoeq->a) / (edoeq->n+1.0);  //largura do passo da malha
    *norma = 1.0 + FLT_EPSILON;
    zera_vetor(X, n);
    //aplica método GaussSeidel
    *tempo = timestamp();
    while (k  < MAXIT && (*norma > FLT_EPSILON))
    {
        for ( i = 0; i < n; ++i)
        {
            xi = edoeq->a + (i+1)*h; // ponto da malha
            Di = 1 - h * edoeq->p(xi)/2.0; // diagonal inferior
            D  = -2 + h*h * edoeq->q(xi); // diagonal principal
            Ds = 1 + h * edoeq->p(xi)/2.0; // diagonal superior
            B = h*h * edoeq->r(xi); // termo independente
            if (i == 0)
                B = B - edoeq->ya * (1 - h*edoeq->p(edoeq->a+h)/2.0); //condição de contorno inicial
            else if (i == n-1)
                B = B - edoeq->yb * (1 + h*edoeq->p(edoeq->b-h)/2.0); //condição de contorno final
            else
                B = B - (Ds * X[i+1] + Di * X[i-1]);
            X[i] = B / D ;//calcula incógnita
        }
        *norma = norma_L2_residuo(X, sl, n);
        ++k;
    }
    *tempo = timestamp() - *tempo;
}

/*! 
    \brief imprime o Resumo do GaussSeidel aplicado 
*/
void imprime_resumo_gaussSeidel(int metodo, int n, double tempo, double norma)
{
    printf("\n ****************************************************************\n");
    if(metodo == GS_COM_VETOR)
        printf("Gauss Seidel utilizando vetor para um número de passos N = %d:\n", n);
    else
        printf("Gauss Seidel SEM vetor para um número de passos N = %d:\n", n);
    printf("Tempo de execução: %9g\n", tempo);
    printf("normaL2: %9g\n", norma);
    printf("\n ****************************************************************\n");

}