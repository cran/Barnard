#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void BarnardW(int *a, int *b, int *c, int *d, int *one_sided, double *dp, double *nuisance_vector_x, double *nuisance_vector_y, double *wald_statistic) {
  int i,j;
  int c1  = (*a)+(*c);
  int c2  = (*b)+(*d);
  int n   = c1+c2;
  double irat = (1.0 / (double)(c1) + 1.0 / (double)(c2));
  double pxo = (double)((*a)+(*b)) / (double)(n);
  (*wald_statistic) = (pxo<=0 || pxo>=1) ? 0 : (((double)((*b))/(double)(c2))-((double)((*a))/(double)(c1))) / sqrt(pxo*(1.0-pxo)*(irat));

  double ps  = 1.0+1.0/(*dp);
  double *IJ = (double*)calloc(2.0*(c1+1)*(c2+1),sizeof(double));

  double txo = (*one_sided) ? (*wald_statistic) : fabs(*wald_statistic);
  double tx;
  double px;
  int ccc=0;
  if (*one_sided) {
    if (txo<0) {
      for (i=0; i<=c1; i++) for (j=0; j<=c2; j++) {
	  px = (double)(i+j)/(double)(n);
	  tx = (px<=0 || px>=1) ? 0 : (((double)j/(double)c2)-((double)i/(double)c1)) / sqrt(px*(1.0-px)*(irat));
	  if (tx<=txo) {IJ[ccc++] = i; IJ[ccc++] = j;}
	}
    } else {
      for (i=0; i<=c1; i++) for (j=0; j<=c2; j++) {
	  px = (double)(i+j)/(double)(n);
	  tx = (px<=0 || px>=1) ? 0 : (((double)j/(double)c2)-((double)i/(double)c1)) / sqrt(px*(1.0-px)*(irat));
	  if (tx>=txo) {IJ[ccc++] = i; IJ[ccc++] = j;}
	}
    }
  } else {
    for (i=0; i<=c1; i++) for (j=0; j<=c2; j++) {
	px = (double)(i+j)/(double)(n);
	tx = (px<=0 || px>=1) ? 0 : (((double)j/(double)c2)-((double)i/(double)c1)) / sqrt(px*(1.0-px)*(irat));
	if (fabs(tx)>=txo) {IJ[ccc++] = i; IJ[ccc++] = j;}
      }
  }

  double n1  = lgamma(c1+1);
  double n2  = lgamma(c2+1);
  double p;
  int k, ii;
  for (k=0; k<ps; k++) {
    p = (double)(k)*(*dp);
    nuisance_vector_x[k] = p;
    nuisance_vector_y[k] = 0;

    for (ii=0; ii<ccc; ii+=2) {
      i = IJ[ii];
      j = IJ[ii+1];
      nuisance_vector_y[k] += exp(n1+n2+(double)(i+j)*log(p)+(double)(n-i-j)*log(1.0-p)-(lgamma(i+1)+lgamma(j+1)+lgamma(c1-i+1)+lgamma(c2-j+1)));
    }
  }

  free(IJ);
}

