//Annelyse Schatzmann GRR20151731

#include <stdio.h>
#include <math.h>

#include "utils.h"
#include "ZeroFuncao.h"

double bisseccao (Polinomio p, double a, double b, double eps,int *it, double *raiz){
	double xm, xmN;


	xm = a+b/2.0;
	do{
		xm == xmN;
		if((calcPolinomio_rapido(a) * calcPolinomio_rapido(xm)) > 0){
			a= xm;
		}
		else
			b=xm;
		xmN = a+b/2.0;
	}while(fabs((xmN - xm)/xmN) > eps);
}


double newtonRaphson (Polinomio p, double x0, double eps,int *it, double *raiz){
	int k=1;
	double x, erro;



	while(fabs((x-x0)/x) > eps && k < it){
		x0 = x;
		k++;

	}

	*raiz = x;
	if ( k >= it){
		return (-1);
	}
	else
		return (0);

}


double secante (Polinomio p, double x0, double x1, double eps,int *it, double *raiz)
{

}


void calcPolinomio_rapido(Polinomio p, double x, double *px, double *dpx){
 	double b = p.p[p.grau];
 	double c = b;
 	int k;

 	for(k = p.grau - 1; k; k--){
 		b = p.p[k] + b * x;
 		c = b + c * x;
  	}

  	b = p.p[0] + b * x;
  	*px = b;
  	*dpx = c;
}
/* b = polinômio
   c =  derivada */


void calcPolinomio_lento(Polinomio p, double x, double *px, double *dpx){
	double b = p.p[p.grau];
 	double c = b;
 	int k;

 	for(k = p.grau-1; k; k--){
 		b = pow(p.p[k], k) + b * x;
 		c = b + c * x;
  	}

  	b = p.p[0] + b * x;
  	*px = b;
  	*dpx = c;
}










