//Annelyse Schatzmann GRR20151731

#include <stdio.h>
#include <math.h>

#include "utils.h"
#include "ZeroFuncao.h"

int main (int argc, char **argv){
	double x= 2;
	double px, dpx, raizB, raizN, raizS, raizBL, raizNL, raizSL;
	int itB, itN, itS, itBL, itNL, itSL;
	double x0,x1,a,b;
	Polinomio p;

	printf("Digite grau do polinômio: ");
    scanf("%d", &p.grau);

  p.p = (double*)malloc(sizeof(double*) * (p.grau + 1)); // porque tem o grau 0
  
  printf("Digite os coeficientes do polinômio: ");
  for (int i = p.grau; i > -1; --i){ 
    scanf("%lf", &p.p[i]);
  }
	
	printf("Digite o valor de x0: ");
    scanf("%lf", &x0);
  printf("Digite o valor de x1: ");
    scanf("%lf", &x1);
  printf("Digite o valor de a: ");
    scanf("%lf", &a);
  printf("Digite o valor de b: ");
    scanf("%lf", &b);
  printf("\n");


/*	calcPolinomio_rapido(p, x, &px, &dpx);
	printf("RAPIDO %lf  %lf\n", px, dpx);

	calcPolinomio_lento(p, x, &px, &dpx);
	printf("LENTO %lf %lf\n", px, dpx);*/


	double temp_ib = timestamp();
	double erroB = bisseccao (p, a, b, EPS, &itB, &raizB);
	double temp_fb = timestamp();
	double tempoB = temp_fb - temp_ib;
	

	double temp_in = timestamp();
	double erroN = newtonRaphson (p, x0, EPS, &itN, &raizN);
	double temp_fn = timestamp();
	double tempoN = temp_fn - temp_in;


	double temp_is = timestamp();
	double erroS = secante (p, x0, x1, EPS, &itS, &raizS);
	double temp_fs = timestamp();
	double tempoS = temp_fs - temp_is;

	printf("+---------------+--------------------------+--------------------------+-------------+-----------+\n");
  printf("|                                  ::Cálculo Poli Rápido::                                      |\n");
	printf("+---------------+--------------------------+--------------------------+-------------+-----------+\n");
  printf("|    Método     |          Raíz            |          Erro            |  Iterações  |   Tempo   |\n");
  printf("+---------------+--------------------------+--------------------------+-------------+-----------+\n");
  printf("|  Bissecção    | %.22lf | %.22lf |     %4d    |  %lf |\n",raizB, erroB, itB, tempoB);
	printf("| NewtonRaphson | %.22lf | %.22lf |     %4d    |  %lf |\n",raizN, erroN, itN, tempoN);
  printf("|    Secante    | %.22lf | %.22lf |     %4d    |  %lf |\n",raizS, erroS, itS, tempoS);
	printf("+---------------+--------------------------+--------------------------+-------------+-----------+\n\n");


	double tempL_ib = timestamp();
	double erroBL = bisseccao_lento(p, a, b, EPS, &itBL, &raizBL);
	double tempL_fb = timestamp();
	double tempoBL = tempL_fb - tempL_ib;
	

	double tempL_in = timestamp();
	double erroNL = newtonRaphson_lento(p, x0, EPS, &itNL, &raizNL);
	double tempL_fn = timestamp();
	double tempoNL = tempL_fn - tempL_in;


	double tempL_is = timestamp();
	double erroSL = secante_lento(p, x0, x1, EPS, &itSL, &raizSL);
	double tempL_fs = timestamp();
	double tempoSL = tempL_fs - tempL_is;


	printf("+---------------+--------------------------+--------------------------+-------------+-----------+\n");
  printf("|                                  ::Cálculo Poli Lento::                                       |\n");
	printf("+---------------+--------------------------+--------------------------+-------------+-----------+\n");
  printf("|    Método     |          Raíz            |          Erro            |  Iterações  |   Tempo   |\n");
  printf("+---------------+--------------------------+--------------------------+-------------+-----------+\n");
  printf("|  Bissecção    | %.22lf | %.22lf |     %4d    |  %lf |\n",raizBL, erroBL, itBL, tempoBL);
	printf("| NewtonRaphson | %.22lf | %.22lf |     %4d    |  %lf |\n",raizNL, erroNL, itNL, tempoNL);
  printf("|    Secante    | %.22lf | %.22lf |     %4d    |  %lf |\n",raizSL, erroSL, itSL, tempoSL);
	printf("+---------------+--------------------------+--------------------------+-------------+-----------+\n\n");
 
  return 0;
}

