#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "EDO.h"
#include "auxiliar.h"

int main ()
{
    Edo* edeoq; SL_Tridiag* sl;
    int n = 5;
    double norma; 
    double *x = malloc( n  * sizeof(double));
    edeoq = malloc (sizeof(Edo)); sl = malloc (sizeof(SL_Tridiag));
    inicializa_EDO(edeoq, n);
    //imprime_EDO (edeoq);
    aloca_tri_diagonal(edeoq, sl);
    gera_tri_diagonal(edeoq, sl);
    imprime_SL(edeoq, sl);
    gaussSeidel(sl->D, sl->Di, sl->Ds, sl->B, x, n, &norma, sl);
    
    printf("o vetor solução para n = 5 é:\n");
    imprime_vetor(x,n);
    printf("a norma n = 5 é:%f \n", norma);

}
