#include "auxiliar.h"


/*!
    \brief imprime um vetor de double

    \param vetor vetor de double
    \param n tamanho do vetor 
*/
void imprime_vetor(double *vetor, int n)
{
    for (int i = 0; i < n; i++)
    {
        printf("%f ", vetor[i]);
    }
    printf("\n");
}