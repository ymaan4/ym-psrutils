/*
 * Copyright (c) 2018  Yogesh Maan <ymaan4@gmail.com>
 * 
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of version 2 of the GNU General Public License as
 * published by the Free Software Foundation.
 * 
 * This program is distributed in the hope that it would be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * 
 * Further, this software is distributed without any warranty that it is
 * free of the rightful claim of any third person regarding infringement
 * or the like.  Any license provided herein, whether implied or
 * otherwise, applies only to this software file.  Patent licenses, if
 * any, provided herein do not apply to combinations of this program with
 * other software, or any other product whatsoever.
 * 
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write the Free Software Foundation, Inc., 59
 * Temple Place - Suite 330, Boston MA 02111-1307, USA.
 * 
 */

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "utils.h"


/*------------------ Preamble functions -------------------------*/
int help_required(char *string)
{
  if (strings_equal(string,"--help")) return(1);
  if (strings_equal(string,"-h")) return(1);
  return(0);
}
int file_exists(char *filename)
{
  if ((fopen(filename,"rb"))==NULL) { return(0);}
  else { return(1); }
}
/*----------*/
FILE *open_file(char *filename, char *descriptor)
{
  FILE *fopen(), *fptr;
  if ((fptr=fopen(filename,descriptor)) == NULL) {
    fprintf(stderr,"Error in opening file: %s\n",filename);
    exit(1);
  }
  return fptr;
}
/*----------*/
void error_message(char *message)
{
  fprintf(stderr,"ERROR: %s\n",message);
  exit(1);
}
/*-------------------------------------- ------------------------*/

//============== Conventional mean rms ===============================
/* TO COMPUTE MEAN and RMS OF 'arr' */

void simple_meanrms(double *arr, long int np)
{
        long int i;
        double amean,rms,diff,avar;

        amean = 0.0;
        for (i=0; i<np; i++){
          amean=arr[i]+amean;
        }
        amean=amean/np;

        rms=0.0;
        for (i=0;i<np;i++){
          diff=arr[i]-amean;
          rms=rms+diff*diff;
        }
        avar=rms/(np-1.0);
        rms=sqrt(rms/(np-1.0));

        arr[np] = amean;
        arr[np+1] = rms;
        arr[np+2] = avar;
        return;
}
//============================================================================



//============== robust mean rms ===============================
/* TO COMPUTE MEAN and RMS OF 'arr' BY EXCLUDING WHAT MAY BE
SOME CONTRIBUTION FROM INTERFERENCE */

void robust_meanrms(double *arr, long int np)
{
        long int i,iter,maxiter;
        double an,amean,amean0,rms,rms0,diff,thresh;

        simple_meanrms(arr,np);
        amean=arr[np];
        rms=arr[np+1];

        iter=0;
        maxiter=10;
        amean0 = amean;
        amean = 0.0;
        an = 0.0;
        for (iter=1; iter<maxiter; iter++) {
          thresh = 4.0*rms;
          for (i=0; i<np; i++){
             diff=fabs(arr[i]-amean0);
             if(diff<=thresh){
               amean=arr[i]+amean;
               an=an+1.0;
             }
          }
          if(an>0)amean=amean/an;
          rms0=rms;

          rms=0.0;
          an=0.0;
          for (i=0;i<np;i++){
            diff=fabs(arr[i]-amean);
            if(diff <= thresh){
              rms=rms+diff*diff;
              an=an+1.0;
            }
          }
          if(an > 0.0)rms=sqrt(rms/an);
          amean0=amean;

          arr[np] = amean;
          arr[np+1] = rms;
          if(rms == 0.0){
             arr[np+1] = 100000.0;
             return;
          }
          if(fabs((rms0/rms)-1.0) < 0.05) return;

        }
        return;
}
//==============================================================


