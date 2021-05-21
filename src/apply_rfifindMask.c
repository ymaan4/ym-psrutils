/*
 * 
 */

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <unistd.h>
#include "utils.h"
#include "utils.h"

struct mask {
  double timesigma,freqsigma,mjd,dtint,lofreq,dfreq;
  int numchan,numint,ptsperint,num_zap_chans,num_zap_ints;
  int *num_chans_per_int;
  int *zap_chans,*zap_ints,**chans;
};
struct mask read_mask(char *filename);
float magic_val=-9999.0;


void apply_rfifindMask_help()
{
  puts("*** This program in-place modifies the data corresponding ***");
  puts("*** to the rfifind-mask. The original data can not be     ***");
  puts("*** recovered after this, so proceed with caution.        ***");
  puts("***                                                       ***");
  puts("");
  puts("");
  puts("apply_rfifindMask - apply rfifind-mask to filterbank data\n");
  puts("usage: apply_rfifindMask -{options} {input-filename} \n");
  puts("options:\n");
  puts("-mask <file> - rfifind generated mask-file name ");
  puts("-rmean       - replace the rejected sections by running mean (instead of 0) ");
  puts("-rmed        - replace the rejected sections by running median (instead of 0) ");
//  puts("-n numbits  - specify output number of bits (def=input-size)");
//  puts("-o filename - specify output filename (def=stdout)");
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

//============== Conventional mean rms ===============================
/* TO COMPUTE MEAN and RMS OF 'arr' */

void simple_meanrms_fl(float *arr, long int np)
{
        long int i;
        float amean,rms,diff,avar;

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

float simple_mean(float *arr, long int np)
{
        long int i;
        float amean;

        amean = 0.0;
        for (i=0; i<np; i++){
          amean=arr[i]+amean;
        }
        amean=amean/np;
        return amean;
}
//============================================================================


//============== robust mean rms ===============================
/* TO COMPUTE MEAN and RMS OF 'arr' BY EXCLUDING WHAT MAY BE
SOME CONTRIBUTION FROM INTERFERENCE -- improved coding */

void robust_meanrms_new(float *arr, long int np)
{
        long int i,iter,maxiter;
        float an,amean,amean0,rms,rms0,diff,thresh;

        simple_meanrms_fl(arr,np);
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
//============== spectrum cleaning =========================================
void spfind(float *cdata, long int npts, float thresh, float *wt)
{
  long int i,j,k,n;
  float anpt,ath,m1,r1,m2;

  //--------------------------
  n = npts-1;
  for (i=0;i<n;i++) ai[i]=fabs(cdata[i+1]-cdata[i]) ;
  ai[n] = 0.0;
  ai[n+1] = 0.0;
  robust_meanrms_new(ai,n);
  m1 = ai[n];
  r1 = ai[n+1];
  if(r1<=0.0) r1 = 10000.0;

  ath = thresh*r1;
  for (i=0;i<n;i++) {
    if( fabs(ai[i]-m1)>ath){
      wt[i] = -1.0;
      wt[i+1] = -1.0;
    }
  }
  //--------------------------
  n = npts-2;
  for (i=0;i<n;i++) ai[i]=fabs(cdata[i+2]-cdata[i]) ;
  ai[n] = 0.0;
  ai[n+1] = 0.0;
  robust_meanrms_new(ai,n);
  m1 = ai[n];
  r1 = ai[n+1];
  if(r1<=0.0) r1 = 10000.0;
  ath = thresh*r1;
  for (i=0;i<n;i++) {
    if( fabs(ai[i]-m1)>ath){
      wt[i] = -1.0;
      wt[i+1] = -1.0;
      wt[i+2] = -1.0;
    }
  }
  //--------------------------
  n = npts-3;
  for (i=0;i<n;i++) ai[i]=fabs(cdata[i+3]-cdata[i]) ;
  ai[n] = 0.0;
  ai[n+1] = 0.0;
  robust_meanrms_new(ai,n);
  m1 = ai[n];
  r1 = ai[n+1];
  if(r1<=0.0) r1 = 10000.0;
  ath = thresh*r1;
  for (i=0;i<n;i++) {
    if( fabs(ai[i]-m1)>ath){
      wt[i] = -1.0;
      wt[i+1] = -1.0;
      wt[i+2] = -1.0;
      wt[i+3] = -1.0;
    }
  }
  //--------------------------
  n = npts-4;
  for (i=0;i<n;i++) ai[i]=fabs(cdata[i+4]-cdata[i]) ;
  ai[n] = 0.0;
  ai[n+1] = 0.0;
  robust_meanrms_new(ai,n);
  m1 = ai[n];
  r1 = ai[n+1];
  if(r1<=0.0) r1 = 10000.0;
  ath = thresh*r1;
  for (i=0;i<n;i++) {
    if( fabs(ai[i]-m1)>ath){
      wt[i] = -1.0;
      wt[i+1] = -1.0;
      wt[i+2] = -1.0;
      wt[i+3] = -1.0;
      wt[i+4] = -1.0;
    }
  }
  //--------------------------
}
void tsfind(float *ttdata, long int npts, float thresh, float *wt)
{
  long int i,j,n;
  float anpt,ath;
  float m1,r1,m2;


  for (i=0;i<npts;i++) ai[i]=ttdata[i] ;
  ai[npts] = 0.0;
  ai[npts+1] = 0.0;
  robust_meanrms_new(ai,npts);
  m1 = ai[npts];
  r1 = ai[npts+1];
  if(r1<=0.0) r1 = 10000.0;

  ath = thresh*r1;
  for (i=0;i<npts;i++) {
    if( fabs((ai[i]-m1))>ath){
      wt[i] = -1.0;
    }
  }
  wt[npts] = m1;

}
//--------------------------------------------------------------
// To show a minimal status bar 
// (adapted from: https://cboard.cprogramming.com/c-programming/99580-progress-bar.html)
int show_status (double percent)
{   
    int x;
    printf("["); 
    for(x = 0; x < percent/2; x++) printf("=");
    for(x = percent/2; x < 100/2; x++) printf(" ");
    
    printf("]");

    printf(" [%d%%]\r", (int)percent);
    fflush(stdout);
    return(EXIT_SUCCESS);
}
//--------------------------------------------------------------


void main (int argc, char *argv[])
{
  int i, j, nc, headersize, channum ;
  int64_t ichan,isamp,isampmin,isampmax,k;
  long istart=0,ii,jj;
  float tstart=0.0,tend,dur=1.0,lg, atemp;
  char string[80];
  float *fblock,min,max, *mmspec,*mspec,*wspec,*replace,*last_replace;
  float *chandata;
  unsigned short *sblock;
  unsigned char  *cblock;
  int nsaved=0,ns=0,nsblk,opened=0,nout,iter,rmean=0,rmed=0,nbl=10;
  long int itemp, isum;
  struct mask msk;

  /* set up default global variables */
  obits=naddt=nsamp=0;

  if (argc > 2) {
    i=1;
    while (i<argc) {
      if (strings_equal(argv[i],"-n")) {
	i++;
	obits=atoi(argv[i]);
      } else if (strings_equal(argv[i],"-mask")) {
	i++;
	msk=read_mask(argv[i]);
      } else if (strings_equal(argv[i],"-rmean")) {
	rmean=1;
      } else if (strings_equal(argv[i],"-rmed")) {
	rmed=1;
      } else if (help_required(argv[1])) {
	apply_rfifindMask_help();
	exit(0);
      } else if (file_exists(argv[i])) {
	strcpy(inpfile,argv[i]);
	input=open_file(inpfile,"r+b");
      } else {
	apply_rfifindMask_help();
	sprintf(string,"unknown argument (%s) passed to apply_rfifindMask",argv[i]);
	error_message(string);
      }
      i++;
    }
  }
  else {
   apply_rfifindMask_help();
   exit(0);
  }

  if (rmean==1 && rmed==1){
    printf("Both 'rmean' and 'rmed' options specified!\n");
    printf("Now proceeding with running median.\n");
  }

  if ((headersize=read_header(input))) {
    totsamp = nsamples(inpfile,headersize,nbits,nifs,nchans);
  } else {
    error_message("input data file is of unknown origin!!!");
  }

  naddt = msk.ptsperint ; 
  if (obits == 0) obits=nbits;

  nsblk=nchans*nifs*naddt;
  fblock=(float *) malloc(nsblk*sizeof(float));
  sblock=(unsigned short *) malloc(nsblk*sizeof(unsigned short));
  cblock=(unsigned char *) malloc(nsblk*sizeof(unsigned char)); 
  wspec=(float *) malloc(2*nchans*sizeof(float));
  mspec=(float *) malloc(2*nchans*sizeof(float));
  mmspec=(float *) malloc(2*nchans*nbl*sizeof(float));
  replace=(float *) malloc(2*nchans*sizeof(float));
  last_replace=(float *) malloc(2*nchans*sizeof(float));
  if(naddt>nchans){
    chandata=(float *) malloc(naddt*2*sizeof(float));
    ai = (float *) malloc(naddt*2*sizeof(float));
  } else {
    chandata=(float *) malloc(nchans*2*sizeof(float));
    ai = (float *) malloc(nchans*2*sizeof(float));
  }
  min=0.0;
  max=(float) pow(2.0,(double)obits) -1.0;
  for (i=0; i<nchans; i++) replace[i] = 0.0;
  for (i=0; i<nchans; i++) last_replace[i] = 0.0;

  if (rmean>0 || rmed>0){
    printf("Determining a representative spectrum using first %d blocks...",nbl);
    for (i=0; i<nbl; i++){
      istart = (long) (headersize + (long)(naddt*nifs*i*nchans*(nbits/8.0)));
      if ((ns=read_block(input,nbits,fblock,nsblk,istart))>0) {
        j=i;
        for (channum=0; channum<nchans; channum++) {
          /* Select the correct channel */
          for (ii=0, jj=channum; ii<naddt; ii++, jj+=nchans) chandata[ii]=(float) fblock[jj];
          if(rmed>0) {
            mmspec[i*nchans+channum] = median(chandata,naddt);
          } else {
            mmspec[i*nchans+channum] = simple_mean(chandata,naddt);
          }
        }
      }
    }
    // get the median/mean of the median/mean spectra
    for (channum=0; channum<nchans; channum++) {
      /* Select the correct channel */
      for (ii=0, jj=channum; ii<j; ii++, jj+=nchans) chandata[ii]=(float) mmspec[jj];
      if(rmed>0) {
        mspec[channum] = median(chandata,j);
        last_replace[channum] = mspec[channum];
      } else {
        mspec[channum] = simple_mean(chandata,j);
        last_replace[channum] = mspec[channum];
      }
    }
    // clean the spiky RFI from the median/mean spectrum 
    atemp=3.0;
    for (i=0;i<nchans; i++) chandata[i] = (float) mspec[i];
    for (i=0; i<nchans; i++) wspec[i]=+1.0;
    spfind(chandata,nchans,atemp,wspec);
    tsfind(chandata,nchans,atemp,wspec);
    for (i=0;i<nchans; i++){
      if(wspec[i] > 0.0) {
        lg = mspec[i];
        break;
      }
    }
    for (i=0;i<nchans;i++){
      if(wspec[i]<0.0) mspec[i] = lg;
      if(wspec[i]>0.0) lg = mspec[i];
    }
    // write out, only for test purposes
    {
    FILE *testOut;
    testOut = fopen("out.ascii","w");
    for (i=0;i<nchans;i++) fprintf(testOut, "%d    %f   %f\n",i,mspec[i],last_replace[i]);
    fclose(testOut); }

    printf(" ...done!\n");
    printf("\n");
  }

  for (i=0; i<msk.numint; i++){
    isampmin = msk.ptsperint*(int64_t) i;
    isampmax = msk.ptsperint*(int64_t) (i+1);
    if (isampmax > totsamp){
       isampmax = totsamp;
       naddt = isampmax - isampmin + 1;
    }
    istart = (long) (headersize + (long)(isampmin*nchans*(nbits/8.0)));
 
    if ((ns=read_block(input,nbits,fblock,nsblk,istart))>0) {
      //----------------------------------------------
      //--- compute optimum replacement values from the current and previous blocks ---
      if(rmean>0 || rmed>0){
        for (ichan=0;ichan<nchans;ichan++){
          for (ii=0,jj=ichan; ii<naddt; ii++,jj+=nchans) chandata[ichan]=fblock[jj];
            if(rmean>0){ replace[ichan] = simple_mean(chandata,naddt);}
            else { replace[ichan] = median(chandata,naddt);}
        }
        for (j=0; j<msk.num_chans_per_int[i]; j++){
          ichan = nchans - msk.chans[i][j] - 1; 
          replace[ichan] = last_replace[ichan];
        }
        for (ichan=0;ichan<nchans;ichan++)last_replace[ichan] = replace[ichan];
      }
      //----------------------------------------------
      //--- now the actual replacements
      for (j=0; j<msk.num_chans_per_int[i]; j++){
        ichan = nchans - msk.chans[i][j] - 1; 
        for (isamp=0; isamp<naddt;isamp++){
          k = ichan + nchans*isamp;
          fblock[k] = replace[ichan];
        }
      }

      nout=ns;
      switch (obits) {
      case 32:
        fseek(input, istart, SEEK_SET);
        fwrite(fblock,sizeof(float),nout,input);
        break;
      case 16:
        fseek(input, istart, SEEK_SET);
        float2short(fblock,nout,min,max,sblock);
        fwrite(sblock,sizeof(unsigned short),nout,input);
        break;
      case 8:
        fseek(input, istart, SEEK_SET);
        float2char(fblock,nout,min,max,cblock);
        fwrite(cblock,sizeof(unsigned char),nout,input);
        break;
      } 
    } 
    show_status(i*100.0/msk.numint);
  }

  printf("\n");
  printf("\n");

  free (fblock);
  free (sblock);
  free (cblock);
  free (ai);
  free (mspec);
  free (mmspec);
  free (wspec);
  free (replace);
  free (last_replace);
  free (chandata);
  fclose(input);
  exit(0);
}

// Read PRESTO rfifind mask (*.mask)
struct mask read_mask(char *filename)
{
  struct mask msk;
  FILE *file;
  int i,j;

  // Open file
  file=fopen(filename,"rb");

  // Read information
  fread(&msk.timesigma,sizeof(double),1,file);
  fread(&msk.freqsigma,sizeof(double),1,file);
  fread(&msk.mjd,sizeof(double),1,file);
  fread(&msk.dtint,sizeof(double),1,file);
  fread(&msk.lofreq,sizeof(double),1,file);
  fread(&msk.dfreq,sizeof(double),1,file);
  fread(&msk.numchan,sizeof(int),1,file);
  fread(&msk.numint,sizeof(int),1,file);
  fread(&msk.ptsperint,sizeof(int),1,file);

  // Channels
  fread(&msk.num_zap_chans,sizeof(int),1,file);
  msk.zap_chans=(int *) malloc(sizeof(int)*msk.num_zap_chans);
  fread(msk.zap_chans,sizeof(int),msk.num_zap_chans,file);

  // Subints
  fread(&msk.num_zap_ints,sizeof(int),1,file);
  msk.zap_ints=(int *) malloc(sizeof(int)*msk.num_zap_ints);
  fread(msk.zap_ints,sizeof(int),msk.num_zap_ints,file);

  // Full mask
  msk.num_chans_per_int=(int *) malloc(sizeof(int)*msk.numint);
  fread(msk.num_chans_per_int,sizeof(int),msk.numint,file);

  msk.chans=(int **) malloc(sizeof(int *)*msk.numint);
  for (i=0;i<msk.numint;i++) {
    if (msk.num_chans_per_int[i]>0 && msk.num_chans_per_int[i]<msk.numchan) {
      msk.chans[i]=(int *) malloc(sizeof(int)*msk.num_chans_per_int[i]);
      fread(msk.chans[i],sizeof(int),msk.num_chans_per_int[i],file);
    } else if (msk.num_chans_per_int[i]==msk.numchan) {
      msk.chans[i]=(int *) malloc(sizeof(int)*msk.num_chans_per_int[i]);
      for (j=0;j<msk.numchan;j++)
        msk.chans[i][j]=j;
    }
  }

  // Close file
  fclose(file);

  return msk;
}

