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
/*!
    \brief zera um vetor de double

    \param vetor vetor de double
    \param n tamanho do vetor 
*/
void zera_vetor(double *vetor, int n)
{
    for (int i = 0; i < n; i++)
        vetor[i] = 0.0;
    
}

/*  Retorna tempo em milisegundos desde EPOCH

    Forma de uso:
 
    double tempo;
    tempo = timestamp();
    <trecho de programa do qual se deseja medir tempo>
    tempo = timestamp() - tempo;
*/

double timestamp (void)
{
  struct timespec tp;
  clock_gettime(CLOCK_MONOTONIC_RAW, &tp);
  return ( (double) tp.tv_sec*1.0e3 + (double) tp.tv_nsec*1.0e-6 );
}