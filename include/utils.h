#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "header.h"
//#include <fftw3.h>
long int nsamp,naddc,naddt,totsamp,tot_now,nout,n0;
int headerless,obits,iflip,nred, nints,nsub,nsubchans,tnorm,fnorm,nsect;
char inpfile[128], outfile[128], gminfofile[128],gmhdrfile[128],psfile[100], listfile[128];
FILE *input, *output, *output2, *gminfo, *gmhdr, *list;
double tempra,tempdec;
float *ai;

int strings_equal (char *string1, char *string2);
void simple_meanrms(double *arr, long int np);
void robust_meanrms(double *arr, long int np);
void bcast_header();
int help_required(char *string);
double nrselect (long int k, long int n, double arr[]);
double dmedian(double *arr, long int np);
float median(float *arr, long int np);
void medrms(double *arr, long int np);
int file_exists(char *filename);
FILE *open_file(char *filename, char *descriptor);

/* include these from sigproc-4.3 */
void slaCldj ( int iy, int im, int id, double *djm, int *j );
char *backend_name (int machine_id) ;
char *data_category (int data_type) ;
char *headername (char *filename) ;
char *telescope_name (int telescope_id) ;
char tempo_site(int telescope_id) ;
double mjd(int year, int month, int day) ;
int read_block(FILE *input, int nbits, float *block, int nread, long istart);
int read_block_orig(FILE *input, int nbits, float *block, int nread);
int read_header(FILE *inputfile) ;
int read_gmheader(char gminfofile[], char gmhdrfile[]) ;
int typeof_inputdata(FILE *fptr, char *filename) ;
long long nsamples(char *filename,int headersize, int nbits, int nifs, int nchans) ;
long long sizeof_file(char name[]) ;
void error_message(char *message) ;
unsigned char charof2ints (int i, int j) ;
void char2ints (unsigned char c, int *i, int *j) ;
void char2fourints (unsigned char c, int *i, int *j, int *k, int *l);
void float2char(float *f, int n, float min, float max, unsigned char *c) ;
void float2four(float *f, int n, float min, float max, unsigned char *c) ;
void float2int(float *f, int n, int b, float min, float max, int *i) ;
void float2short(float *f, int n, float min, float max, unsigned short *s) ;
void get_string(FILE *inputfile, int *nbytes, char string[]) ;
void int2float(int *i, int n, int b, float min, float max, float *f) ;
void send_coords(double raj, double dej, double az, double za) ;
void send_double (char *name, double double_precision) ;
void send_float(char *name,float floating_point) ;
void send_int(char *name, int integer) ;
void send_long(char *name, long integer) ;
void send_string(char *string) ;
void swap_double( double *pd ) ;
void swap_float( float *pf ) ;
void swap_int( int *pi ) ;
void swap_longlong( long long *pl ) ;
void swap_long( long *pi ) ;
void swap_short( unsigned short *ps ) ;
void swap_ulong( unsigned long *pi ) ;
