//Annelyse Schatzmann GRR20151731

#include <stdio.h>
#include <math.h>

#include "utils.h"
#include "ZeroFuncao.h"

double bisseccao (Polinomio p, double a, double b, double eps, int *it, double *raiz){
	double xm, pxm, dpxm, xmn, pxmn, dpxmn, pa, dpa, erro;
	int inte = 1; // ja começa com uma iteração

	xmn = (a + b) / 2;

	do {
		xm = xmn;

		calcPolinomio_rapido(p, a, &pa, &dpa);
		calcPolinomio_rapido(p, xmn, &pxmn, &dpxmn);
		
		if((pa * pxmn) > ZERO){
			a = xm;
		} 
		else {
			b = xm;
		}

		xmn = (a + b) / 2;
		erro = fabs((xmn - xm) / xmn);
		inte ++;

	} while(erro > eps && inte < MAXIT);


	*raiz = xmn;
	*it = inte;

	return erro;
}


double newtonRaphson (Polinomio p, double x0, double eps,int *it, double *raiz){
	double x, px0, dpx0, erro;
	int inte = 1;

	calcPolinomio_rapido(p, x0, &px0, &dpx0);

	if(dpx0 != ZERO){  //derivada = 0;
		
		x = x0 - px0 / dpx0;
		erro = fabs((x - x0) / x);

		while(erro > eps && inte < MAXIT){
			x0 = x;
			calcPolinomio_rapido(p, x0, &px0, &dpx0);
			x = x0 - px0 / dpx0;

			erro = fabs((x - x0) / x);
			inte++;
		}

		*raiz = x;
		*it = inte;

		return erro;
	}

	return -1;
}


double secante (Polinomio p, double x0, double x1, double eps, int *it, double *raiz){
	double x, px0, dpx0, px1, dpx1, erro;
	int inte = 1;

	calcPolinomio_rapido(p, x0, &px0, &dpx0);
	calcPolinomio_rapido(p, x1, &px1, &dpx1);
		
	x = x1 - ((px1 * (x1 - x0)) / (px1 - px0));
	erro = fabs((x - x0) / x);

	while(erro > eps && inte < MAXIT){
		x0 = x1;
		x1 = x;
		calcPolinomio_rapido(p, x0, &px0, &dpx0);
		calcPolinomio_rapido(p, x1, &px1, &dpx1);
		
		x = x1 - ((px1 * (x1 - x0)) / (px1 - px0));
		erro = fabs((x - x0) / x);	
		inte++;
	}

	*raiz = x;
	*it = inte;

	return erro;
}


double bisseccao_lento(Polinomio p, double a, double b, double eps, int *it, double *raiz){
	double xm, pxm, dpxm, xmn, pxmn, dpxmn, pa, dpa, erro;
	int inte = 1; // ja começa com uma iteração

	xmn = (a + b) / 2;

	do {
		xm = xmn;

		calcPolinomio_lento(p, a, &pa, &dpa);
		calcPolinomio_lento(p, xmn, &pxmn, &dpxmn);
		
		if((pa * pxmn) > ZERO){
			a = xm;
		} 
		else {
			b = xm;
		}

		xmn = (a + b) / 2;
		erro = fabs((xmn - xm) / xmn);
		inte ++;

	} while(erro > eps && inte < MAXIT);


	*raiz = xmn;
	*it = inte;

	return erro;
}

double newtonRaphson_lento (Polinomio p, double x0, double eps,int *it, double *raiz){
	double x, px0, dpx0, erro;
	int inte = 1;

	calcPolinomio_lento(p, x0, &px0, &dpx0);

	if(dpx0 != ZERO){  //derivada = 0;
		
		x = x0 - px0 / dpx0;
		erro = fabs((x - x0) / x);

		while(erro > eps && inte < MAXIT){
			x0 = x;
			calcPolinomio_lento(p, x0, &px0, &dpx0);
			x = x0 - px0 / dpx0;

			erro = fabs((x - x0) / x);
			inte++;
		}

		*raiz = x;
		*it = inte;

		return erro;
	}
}

double secante_lento (Polinomio p, double x0, double x1, double eps, int *it, double *raiz){
	double x, px0, dpx0, px1, dpx1, erro;
	int inte = 1;

	calcPolinomio_lento(p, x0, &px0, &dpx0);
	calcPolinomio_lento(p, x1, &px1, &dpx1);
		
	x = x1 - ((px1 * (x1 - x0)) / (px1 - px0));
	erro = fabs((x - x0) / x);

	while(erro > eps && inte < MAXIT){
		x0 = x1;
		x1 = x;
		calcPolinomio_lento(p, x0, &px0, &dpx0);
		calcPolinomio_lento(p, x1, &px1, &dpx1);
		
		x = x1 - ((px1 * (x1 - x0)) / (px1 - px0));
		erro = fabs((x - x0) / x);	
		inte++;
	}

	*raiz = x;
	*it = inte;

	return erro;
}


void calcPolinomio_rapido(Polinomio p, double x, double *px, double *dpx){
	double poli, deri;
	int k;

	poli = p.p[p.grau];
	deri = poli;

	for(k = p.grau - 1; k; --k){
		poli = p.p[k] + poli * x;
		deri = poli + deri * x;
	}

	poli = p.p[0] + poli * x;

	*px = poli;
	*dpx = deri;
}


void calcPolinomio_lento(Polinomio p, double x, double *px, double *dpx){
	double poli = 0;
 	double deri = poli;
 	int k;

 	for(k = p.grau; k; --k){
 		poli = p.p[k] * pow(x, k) + poli;
 		deri = (p.p[k] * k) * pow(x, k-1) + deri;
  }

  poli = p.p[0] * pow(x,0) + poli;

  *px = poli;
  *dpx = deri;
}