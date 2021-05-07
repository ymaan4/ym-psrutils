#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "header.h"
long int nsamp,naddc,naddt,totsamp;
int headerless,obits,iflip;
char inpfile[128], outfile[128] ;
FILE *input, *output ;
double tempra,tempdec ;
float *ai;

/* include these from sigproc-4.3 */
int strings_equal (char *string1, char *string2);
void slaCldj ( int iy, int im, int id, double *djm, int *j );
int read_block(FILE *input, int nbits, float *block, int nread, long istart) ;
int read_header(FILE *inputfile) ;
long long nsamples(char *filename,int headersize, int nbits, int nifs, int nchans) ;
long long sizeof_file(char name[]) ;
double nrselect (long int k, long int n, double arr[]);
float median(float *arr, long int np);
unsigned char charof2ints (int i, int j) ;
void char2ints (unsigned char c, int *i, int *j) ;
void char2fourints (unsigned char c, int *i, int *j, int *k, int *l);
void error_message(char *message) ;
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
