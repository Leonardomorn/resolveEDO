#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "EDO.h"

int main ()
{
    Edo* edeoq; SL_Tridiag* sl;
    int n = 5;
    double norma; 
    double *x = malloc( n  * sizeof(double));
    inicializa_EDO(edeoq, n);
    //imprime_EDO (edeoq);
    aloca_tri_diagonal(edeoq, sl);
    gera_tri_diagonal(edeoq, sl);
    gaussSeidel(sl->D, sl->Di, sl->Ds, sl->B, x, n, norma);



}
