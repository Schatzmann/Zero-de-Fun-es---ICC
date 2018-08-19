//Annelyse Schatzmann GRR20151731

#include <stdio.h>
#include <math.h>

#include "utils.h"
#include "ZeroFuncao.h"

double bisseccao (Polinomio p, double a, double b, double eps, int *it, double *raiz){
	double xm, pxm, dpxm, xmn, pxmn, dpxmn, pa, dpa;

	xmn = (a + b) / 2;

	do {
		xm = xmn;
		
		if((calcPolinomio_rapido(p, a, &pa, &dpa) * calcPolinomio_rapido(p, xmn, &pxmn, &dpxmn)) > ZERO){
			a = xm;
		} 
		else {
			b = xm;
		}

		xmn = (a + b) / 2;
	} while(fabs((xmn - xm) / xmn) > eps);

	raiz = xmn;

	return 0;
}


double newtonRaphson (Polinomio p, double x0, double eps,int *it, double *raiz){
	double x, px0, dpx0;

	calcPolinomio_rapido(p, x0, &px0, &dpx0);

	if(dpx0 != ZERO){
		
		x = x0 - px0 / dpx0;

		while(fabs((x - x0) / x) > eps){
			x0 = x;
			calcPolinomio_rapido(p, x0, &px0, &dpx0);
			x = x0 - px0 / dpx0;
		}

		raiz = x;

		return 0;
	}

	return -1;
}


double secante (Polinomio p, double x0, double x1, double eps, int *it, double *raiz){
	double x, px0, dpx0, px1, dpx1;

	calcPolinomio_rapido(p, x0, &px0, &dpx0);
	calcPolinomio_rapido(p, x1, &px1, &dpx1);
		
	x = x1 - ((px1 * (x1 - x0)) / (px1 - px0));

	while(fabs((x - x0) / x) > eps){
		x0 = x1;
		x1 = x;
		calcPolinomio_rapido(p, x0, &px0, &dpx0);
		calcPolinomio_rapido(p, x1, &px1, &dpx1);
		
		x = x1 - ((px1 * (x1 - x0)) / (px1 - px0));
	}

	raiz = x;

	return 0;
}



void calcPolinomio_rapido(Polinomio p, double x, double *px, double *dpx){
void calcPolinomio_rapido(Polinomio p, double x, double *px, double *dpx){
	double poli, deri;
	int k;

	poli = p.p[p.grau];
	deri = poli;

	for(k = p.grau - 1; k; k--){
		poli = p.p[k] + poli * x;
		deri = poli + deri * x;
	}

	poli = p.p[0] + poli * x;

	*px = poli;
	*dpx = deri;
}


void calcPolinomio_lento(Polinomio p, double x, double *px, double *dpx){
	double poli = p.p[p.grau];
 	double deri = poli;
 	int k;

 	for(k = p.grau-1; k; k--){
 		poli = pow(p.p[k], k) + poli * x;
 		deri = poli + deri * x;
  	}

  	poli = p.p[0] + poli * x;
  	*px = poli;
  	*dpx = deri;
}