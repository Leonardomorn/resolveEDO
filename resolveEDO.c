#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "EDO.h"

int main ()
{
    Edo* edeoq;
    int n = 5;
    inicializa_EDO(edeoq, n);
    imprime_EDO (edeoq);

}
