/*
 * Copyright (c) 2018  Yogesh Maan <maan@astron.nl>
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
#include <stdint.h>
#include <math.h>
#include "utils.h"


void normByRMS_help()
{
  puts("");
  puts("normByRMS - Normalize the individual channels data by rms.\n");
  puts("usage: normByRMS -{options} {input-filename(def: stdin)} \n");
  puts("options:\n");
  puts("-o filename  - specify output filename (def=stdout)");
  puts("-t blocksize - block-size to process (def: 4096) ");
//  puts("-n numbits  - specify output number of bits (def=input-size)");
  puts("");
}


//--------------------------------------------------------------
int all_same(float *arr, int np)
{
   int i;
   float atemp;
   atemp = arr[0];
   for (i=1; i<np; i++){
     if(atemp != arr[i]) return 0;
   }
   return 1;
}

//============================================================================
//--------------------------------------------------------------

void main (int argc, char *argv[])
{
  int i, j, nc, headersize, opened_input=0,opened_output=0;
  int64_t ichan,isamp,isampmin,isampmax,k;
  long istart=0,ii,jj;
  float tstart=0.0,tend,dur=1.0,lg;
  char string[80];
  float *fblock,min,max;
  double *chandata;
  unsigned short *sblock;
  unsigned char  *cblock;
  int nsaved=0,ns=0,nsblk,opened=0,nout,iter;
  long int itemp, isum;

  /* set up default global variables */
  obits=nsamp=0;
  naddt=4096;
  input=stdin;
  strcpy(inpfile,"stdin");
  output=stdout;
  strcpy(outfile,"stdout");

  if (argc > 1) {
    i=1;
    while (i<argc) {
      if (strings_equal(argv[i],"-o")) {
	output=fopen(argv[++i],"wb");
        opened_output=1;
      } else if (strings_equal(argv[i],"-n")) {
	i++;
	obits=atoi(argv[i]);
      } else if (strings_equal(argv[i],"-t")) {
        i++;
        naddt=atoi(argv[i]);
      } else if (strings_equal(argv[i],"-headerless")) {
        headerless=1;
      } else if (help_required(argv[1])) {
	normByRMS_help();
	exit(0);
      } else if (file_exists(argv[i])) {
	strcpy(inpfile,argv[i]);
	input=open_file(inpfile,"r+b");
        opened_input=1;
      } else {
	normByRMS_help();
	sprintf(string,"unknown argument (%s) passed to normByRMS",argv[i]);
	error_message(string);
      }
      i++;
    }
  }
  else {
   normByRMS_help();
   exit(0);
  }
  if (!opened_input) {
    /* no input file selected, use standard input */
    input=stdin;
    strcpy(inpfile,"stdin");
  }
  if (!opened_output) {
    /* no output file selected, use standard output */
    output=stdout;
    strcpy(outfile,"stdout");
  }


  if ((headersize=read_header(input))) {
    //totsamp = nsamples(inpfile,headersize,nbits,nifs,nchans);
  } else {
    error_message("input data file is of unknown origin!!!");
  }
  if (obits == 0) obits=nbits;
  if (!headerless) bcast_header();

  nsblk=nchans*nifs*naddt;
  fblock=(float *) malloc(nsblk*sizeof(float));
  sblock=(unsigned short *) malloc(nsblk*sizeof(unsigned short));
  cblock=(unsigned char *) malloc(nsblk*sizeof(unsigned short));
  chandata=(double *) malloc((naddt+4)*sizeof(double));
  min=0.0;
  max=(float) pow(2.0,(double)obits) -1.0;

  while((ns=read_block_orig(input,nbits,fblock,nsblk))>0) {
    n0 = ns/nchans;
    //----------------------------------------------
    for (ichan=0;ichan<nchans;ichan++){
      for (ii=0, jj=ichan;ii<n0; ii++,jj+=nchans) chandata[ii]=fblock[jj];
      medrms(chandata,n0);
      for (ii=0, jj=ichan;ii<n0; ii++,jj+=nchans){
        fblock[jj] = (float)(chandata[ii]-chandata[n0])/chandata[n0+1];
      }
    }
    nout=ns;
    switch (obits) {
    case 32:
      fwrite(fblock,sizeof(float),nout,output);
      break;
    case 16:
      float2short(fblock,nout,min,max,sblock);
      fwrite(sblock,sizeof(unsigned short),nout,output);
      break;
    case 8:
      float2char(fblock,nout,min,max,cblock);
      fwrite(cblock,sizeof(unsigned char),nout,output);
      break;
    } 
  }

  free (fblock);
  free (sblock);
  free (cblock);
  free (chandata);
  fclose(input);
  fclose(output);
  exit(0);
}


