#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include "utils.h"

void medrms(double *arr, long int np)
{
    long int i,iter,maxiter,nn,k;
    double *barr,med,rms,rms0,diff,atemp,an;

    nn = np;

    barr=(double *) malloc(np*sizeof(double));
    for (i=0;i<np;i++) barr[i] = arr[i];
    med = dmedian(barr,nn);
    free(barr);

    atemp = 0.0;
    for (i=0;i<np;i++) atemp = atemp + (arr[i]-med)*(arr[i]-med);
    rms = sqrt(atemp/np);
    iter=0;
    maxiter=10;
    rms0=rms*2.0;
    while (iter<maxiter && fabs((rms0/rms)-1.0) > 0.05){
      rms0 = rms;
      atemp = 0.0;
      an = 0.0;
      for (i=0;i<np;i++){
        diff = arr[i]-med;
        if(fabs(diff) <= 3.5*rms0){
          an = an + 1.0;
          atemp = atemp + diff*diff;
        }
      }
      rms = sqrt(atemp/an);
      iter = iter+1;
    }
    arr[np] = med;
    arr[np+1] = rms;

}



