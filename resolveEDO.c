#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "EDO.h"
#include "auxiliar.h"

int main ()
{
    Edo* edeoq; SL_Tridiag* sl;
    int n;
    double norma; 
    double *X;
    double tempo;
    edeoq = malloc (sizeof(Edo)); sl = malloc (sizeof(SL_Tridiag));
    
    //resolve para número de passos igual a 5
    n = 5;
    X = malloc( n  * sizeof(double));
    inicializa_EDO(edeoq, n);
    aloca_tri_diagonal(edeoq, sl);
    gera_tri_diagonal(edeoq, sl);
    imprime_SL(edeoq, sl);
    

    gaussSeidel_vetor(sl->D, sl->Di, sl->Ds, sl->B, X, n, &norma, sl, &tempo);

    imprime_resumo_gaussSeidel(GS_COM_VETOR, n, tempo, norma);

    imprime_EDO(edeoq);
    gaussSeidel_direto(edeoq, sl, X, &norma, &tempo);
    imprime_resumo_gaussSeidel(GS_SEM_VETOR, n, tempo, norma);

    libera_tri_diagonal(edeoq, sl);
    free(X);
    //resolve para número de passos igual a 10
    n = 10;
    X = malloc( n  * sizeof(double));
    ajusta_n_EDO(edeoq, n);
    aloca_tri_diagonal(edeoq,sl);
    gera_tri_diagonal(edeoq,sl);
    imprime_SL(edeoq, sl);
    
    gaussSeidel_vetor(sl->D, sl->Di, sl->Ds, sl->B, X, n, &norma, sl, &tempo);
    imprime_resumo_gaussSeidel(GS_COM_VETOR, n, tempo, norma);
    
    gaussSeidel_direto(edeoq, sl, X, &norma, &tempo);
    imprime_resumo_gaussSeidel(GS_SEM_VETOR, n, tempo, norma);

    libera_tri_diagonal(edeoq, sl);
    free(X);
    //resolve para número de passos igual a 100
    n = 100;
    ajusta_n_EDO(edeoq, n);
    aloca_tri_diagonal(edeoq,sl);
    gera_tri_diagonal(edeoq,sl);
    imprime_SL(edeoq, sl);
    
    gaussSeidel_vetor(sl->D, sl->Di, sl->Ds, sl->B, X, n, &norma, sl, &tempo);
    imprime_resumo_gaussSeidel(GS_COM_VETOR, n, tempo, norma);
    
    gaussSeidel_direto(edeoq, sl, X, &norma, &tempo);
    imprime_resumo_gaussSeidel(GS_SEM_VETOR, n, tempo, norma);

    libera_tri_diagonal(edeoq, sl);
}
