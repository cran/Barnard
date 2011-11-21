#include "R.h"
#include "Rmath.h"

void BarnardW(int *a, int *b, int *c, int *d, double *dp, double *nuisance_vector_x, double *nuisance_vector_y0, double *nuisance_vector_y1, double *wald_statistic, double *nuisance_parameter, double *p_value) {
  int i,j;
  int c1  = (*a)+(*c);
  int c2  = (*b)+(*d);
  int n   = c1+c2;
  double irat = (1.0 / (double)(c1) + 1.0 / (double)(c2));
  double pxo = (double)((*a)+(*b)) / (double)(n);
  (*wald_statistic) = (pxo<=0 || pxo>=1) ? 0 : (((double)((*b))/(double)(c2))-((double)((*a))/(double)(c1))) / sqrt(pxo*(1.0-pxo)*(irat));

  double ps  = 1.0+1.0/(*dp);
  double *IJ = (double*)Calloc(3.0*(c1+1)*(c2+1),double);

  double txo = *wald_statistic;
  double tx;
  double px;
  int ccc=0;

  for (i=0; i<=c1; i++) for (j=0; j<=c2; j++) {
      px = (double)(i+j)/(double)(n);
      tx = (px<=0 || px>=1) ? 0 : (((double)j/(double)c2)-((double)i/(double)c1)) / sqrt(px*(1.0-px)*(irat));
      if (i==(*a) && j==(*b)) {
	IJ[ccc++] = i; IJ[ccc++] = j; IJ[ccc++] = 1;
      } else if (c1==c2 && i==(*b) && j==(*a)) {
	IJ[ccc++] = i; IJ[ccc++] = j; IJ[ccc++] = 0;
      } else if (txo<0) {
	if (tx<txo) {IJ[ccc++] = i; IJ[ccc++] = j; IJ[ccc++] = 1;}
	if (tx>-txo) {IJ[ccc++] = i; IJ[ccc++] = j; IJ[ccc++] = 0;}
      } else {
	if (tx>txo) {IJ[ccc++] = i; IJ[ccc++] = j; IJ[ccc++] = 1;}
	if (tx<-txo) {IJ[ccc++] = i; IJ[ccc++] = j; IJ[ccc++] = 0;}
      }
    }

  double n1  = lgamma(c1+1);
  double n2  = lgamma(c2+1);
  double p, ad=0;
  int k, ii;

  nuisance_parameter[0] = 0.0; nuisance_parameter[1] = 0.0;
  p_value[0] = 0.0; p_value[1] = 0.0;
  for (k=0; k<ps; k++) {
    p = (double)(k)*(*dp);
    nuisance_vector_x[k] = p;
    nuisance_vector_y0[k] = 0;
    nuisance_vector_y1[k] = 0;

    for (ii=0; ii<ccc; ii+=3) {
      i = IJ[ii];
      j = IJ[ii+1];
      ad = exp(n1+n2+(double)(i+j)*log(p)+(double)(n-i-j)*log(1.0-p)-(lgamma(i+1)+lgamma(j+1)+lgamma(c1-i+1)+lgamma(c2-j+1)));
      if (IJ[ii+2]) nuisance_vector_y0[k] += ad;
      nuisance_vector_y1[k] += ad;
    }

    if (nuisance_vector_y0[k] >= p_value[0]) {
      nuisance_parameter[0] = nuisance_vector_x[k];
      p_value[0] = nuisance_vector_y0[k];
    }
    if (nuisance_vector_y1[k] >= p_value[1]) {
      nuisance_parameter[1] = nuisance_vector_x[k];
      p_value[1] = nuisance_vector_y1[k];
    }
  }

  Free(IJ);
}

