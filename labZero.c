#include <stdio.h>
#include <math.h>

#include "utils.h"
#include "ZeroFuncao.h"

int main (int argc, char **argv){
	double x = 2;
	double px, dpx;
	Polinomio p;

	p.p = malloc(2);

	p.p[0] = 3;
	p.p[1] = 0;
	p.p[2] = 2;
	p.grau = 2;



	calcPolinomio_rapido(p, x, &px, &dpx);
	printf("RAPIDO %lf  %lf\n", px, dpx);


	calcPolinomio_lento(p, x, &px, &dpx);

	/*bisseccao (p, double a, double b, double eps, int *it, double *raiz);
	newtonRaphson (p, double x0, double eps, int *it, double *raiz);
	secante (p, double x0, double x1, double eps, int *it, double *raiz);
*/


	printf("LENTO %lf   %lf\n", px, dpx);

  return 0;
}

