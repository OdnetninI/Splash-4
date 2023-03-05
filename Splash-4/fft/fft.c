/****************************************************************************/
/*                                                                          */
/*  Copyright (c) 2023 Eduardo Jose Gomez-Hernandez (University of Murcia)  */       
/*  Copyright (c) 1994 Stanford University                                  */
/*                                                                          */
/*  All rights reserved.                                                    */
/*                                                                          */
/*  Permission is given to use, copy, and modify this software for any      */
/*  non-commercial purpose as long as this copyright notice is not          */
/*  removed.  All other uses, including redistribution in whole or in       */
/*  part, are forbidden without prior written permission.                   */
/*                                                                          */
/*  This software is provided with absolutely no warranty and no            */
/*  support.                                                                */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/*  Perform 1D fast Fourier transform using six-step FFT method             */
/*                                                                          */
/*  1) Performs staggered, blocked transposes for cache-line reuse          */
/*  2) Roots of unity rearranged and distributed for only local             */
/*     accesses during application of roots of unity                        */
/*  3) Small set of roots of unity elements replicated locally for          */
/*     1D FFTs (less than root N elements replicated at each node)          */
/*  4) Matrix data structures are padded to reduce cache mapping            */
/*     conflicts                                                            */
/*                                                                          */
/*  Command line options:                                                   */
/*                                                                          */
/*  -mM : M = even integer; 2**M total complex data points transformed.     */
/*  -pNumThreads : NumThreads = number of processors; Must be a power of 2. */
/*  -nN : N = number of cache lines.                                        */
/*  -lL : L = Log base 2 of cache line length in bytes.                     */
/*  -y  : Enable touching the data before operations.                       */
/*        It is not useful when measuring the full application or the ROI.  */
/*  -t  : Perform FFT and inverse FFT.  Test output by comparing the        */
/*        integral of the original data to the integral of the data         */
/*        that results from performing the FFT and inverse FFT.             */
/*  -o  : Print out complex data points.                                    */
/*  -h  : Print out command line options.                                   */
/*                                                                          */
/****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <unistd.h>
#include <stdbool.h>

#include "../common/common.h"

#define PAGE_SIZE               4096
#define NUM_CACHE_LINES        65536 
#define LOG2_LINE_SIZE             4
#define PI                         3.1416
#define DEFAULT_M                 10
#define DEFAULT_THREADS            1

MAIN_ENV();

#define SWAP_VALS(a,b) {			\
    double tmp;					\
    tmp=a;					\
    a=b;					\
    b=tmp;					\
  }						\

#define IF_EXIT(c, msg) {			\
    if (c) {					\
      fprintf(stderr, msg);			\
      exit(-1);					\
    }						\
  }						\

// Global data
struct GlobalMemory {
  long id;
  BARDEC(barrier)
} *Global;

// Running Options
bool test_result = false;
bool doprint = false;
bool is_output = false;
bool doprefetch = false;
long NumThreads = DEFAULT_THREADS;

// Array sizes
long M = DEFAULT_M;
long N;                  /* N = 2^M                                */
long rootN;              /* rootN = N^1/2                          */

// Data arrays
double* x;              /* x is the original time-domain data     */
double* trans;          /* trans is used as scratch space         */
double* umain;          /* umain is roots of unity for 1D FFTs    */
double* umain2;         /* umain2 is entire roots of unity matrix */

long orig_num_lines = NUM_CACHE_LINES;     /* number of cache lines */
long num_cache_lines = NUM_CACHE_LINES;    /* number of cache lines */
long log2_line_size = LOG2_LINE_SIZE;
long pad_length;

// Const
const size_t size_of_an_element = sizeof(double) << 1;

/*
  Forward declarations
*/
void ArgumentParser(int argc, char* argv[]);
void ParallelMain(void);
double TouchArray(double* x, double* scratch, double* u, double* upriv, long MyFirst, long MyLast);
double CheckSum(double* x);
void InitX(double* x);
void InitU(long N, double* u);
void InitU2(long N, double* u, long n1);
long BitReverse(long M, long k);
void FFT1D(long direction, long M, long N, double* x, double* scratch, double* upriv, double* umain2, long thread_id, long MyFirst, long MyLast, long pad_length, bool test_result);
void TwiddleOneCol(long direction, long n1, long j, double* u, double* x, long pad_length);
void Scale(long n1, long N, double* x);
void Transpose(long n1, double* src, double* dest, long thread_id, long MyFirst, long MyLast, long pad_length);
void CopyColumn(long n1, double* src, double* dest);
void Reverse(long N, long M, double* x);
void FFT1DOnce(long direction, long M, long N, double* u, double* x);
void PrintArray(long N, double* x);
void printerr(char *str);

/*
  Main
*/
int main(int argc, char *argv[]) {    
  ArgumentParser(argc, argv);
  long m1 = M >> 1;
  
  MAIN_INITENV(,80000000);
  
  N = 1 << M;
  rootN = 1 << (M >> 1);
  long rowsperproc = rootN / NumThreads;
  IF_EXIT(rowsperproc == 0, "Matrix not large enough. 2**(M/2) must be >= NumThreads\n");
  
  long line_size = 1 << log2_line_size;
  if (line_size < size_of_an_element) {
    printf("WARNING: Each element is a complex double (%ld bytes)\n",size_of_an_element);
    printf("  => Less than one element per cache line\n");
    printf("     Computing transpose blocking factor\n");
    long factor = size_of_an_element / line_size;
    num_cache_lines = orig_num_lines / factor;
  }
  
  pad_length = line_size / size_of_an_element;
  if (line_size <= size_of_an_element) {
    pad_length = 1;
  }
  
  if (rowsperproc * rootN * size_of_an_element >= PAGE_SIZE) {
    long pages = (pad_length * size_of_an_element * rowsperproc) / PAGE_SIZE;
    if (pages * PAGE_SIZE != pad_length * size_of_an_element * rowsperproc) {
      pages ++;
    }
    pad_length = (pages * PAGE_SIZE) / (size_of_an_element * rowsperproc);
  }
  else {
    pad_length = (PAGE_SIZE - (rowsperproc * rootN * size_of_an_element)) / (size_of_an_element * rowsperproc);
    IF_EXIT(pad_length * (size_of_an_element * rowsperproc) !=
	    (PAGE_SIZE - (rowsperproc * rootN * size_of_an_element)),
	    "Padding algorithm unsuccessful\n");
  }

  Global = (struct GlobalMemory *) G_MALLOC(sizeof(struct GlobalMemory));

  x = (double* ) G_MALLOC((N+rootN*pad_length)*size_of_an_element+PAGE_SIZE);
  trans = (double* ) G_MALLOC((N+rootN*pad_length)*size_of_an_element+PAGE_SIZE);
  umain = (double* ) G_MALLOC(rootN*size_of_an_element);  
  umain2 = (double* ) G_MALLOC((N+rootN*pad_length)*size_of_an_element+PAGE_SIZE);
  
  x = (double* ) (((unsigned long) x) + PAGE_SIZE - ((unsigned long) x) % PAGE_SIZE);
  trans = (double* ) (((unsigned long) trans) + PAGE_SIZE - ((unsigned long) trans) % PAGE_SIZE);
  umain2 = (double* ) (((unsigned long) umain2) + PAGE_SIZE - ((unsigned long) umain2) % PAGE_SIZE);

   /*
     In order to optimize data distribution, the data structures x, trans, 
     and umain2 have been aligned so that each begins on a page boundary. 
     This ensures that the amount of padding calculated by the program is 
     such that each processor's partition ends on a page boundary, thus 
     ensuring that all data from these structures that are needed by a 
     processor can be allocated to its local memory
   */

   /*
     POSSIBLE ENHANCEMENT:  Here is where one might distribute the x,
     trans, and umain2 data structures across physically distributed 
     memories as desired.
     
     One way to place data is as follows:
     
     long i = ((N/NumThreads)+(rootN/NumThreads)*pad_length)*2;
     double* base = &(x[0]);
     for (long j=0;j<NumThreads;j++) {
      Place all addresses x such that (base <= x < base+i) on node j
      base += i;
     }
     
     The trans and umain2 data structures can be placed in a similar manner.
   */

  printf("\n");
  printf("FFT with Blocking Transpose\n");
  printf("   %ld Complex Doubles\n",N);
  printf("   %ld Processors\n", NumThreads);
  if (num_cache_lines != orig_num_lines) {
    printf("   %ld Cache lines\n",orig_num_lines);
    printf("   %ld Cache lines for blocking transpose\n",num_cache_lines);
  } else {
    printf("   %ld Cache lines\n",num_cache_lines);
  }
  printf("   %d Byte line size\n",(1 << log2_line_size));
  printf("   %d Bytes per page\n",PAGE_SIZE);
  printf("\n");

  BARINIT(Global->barrier, NumThreads);
  
  Global->id = 0;
  InitX(x);                  /* place random values in x */

  double ck1 = 0;
  if (test_result) {
    ck1 = CheckSum(x);
  }
  if (doprint) {
    printf("Original data values:\n");
    PrintArray(N, x);
  }

  InitU(N,umain);               /* initialize u arrays*/
  InitU2(N,umain2,rootN);

  // Parallel Region
  SPLASH4_ROI_BEGIN();
  CREATE(ParallelMain, NumThreads);
  WAIT_FOR_END(NumThreads);
  SPLASH4_ROI_END();

  // Clossing and printing
  if (doprint) {
    printf("Data values after%s FFT:\n", test_result ? " inverse" : "");
    PrintArray(N, x);
  }

  if (test_result) {
    double ck3 = CheckSum(x);
    printf("              INVERSE FFT TEST RESULTS\n");
    printf("Checksum difference is %.3f (%.3f, %.3f)\n",
	   ck1-ck3, ck1, ck3);
    if (fabs(ck1-ck3) < 0.001) {
      printf("TEST PASSED\n");
    } else {
      printf("TEST FAILED\n");
    }
  }

  MAIN_END();
}

void ArgumentParser(int argc, char* argv[]) {
  long c;
  extern char *optarg;
  
  while ((c = getopt(argc, argv, "p:m:n:l:stoh")) != -1) {
    switch(c) {
      case 'p':
	NumThreads = atoi(optarg);
	IF_EXIT(NumThreads < 1, "NumThreads must be >= 1\n");
	IF_EXIT(NumThreads != (1 << LOG2(NumThreads)), "NumThreads must be a power of 2\n");
	break;  
      case 'm':
	M = atoi(optarg);
	long m1 = M >> 1;
	IF_EXIT((m1 << 1) != M, "M must be even\n");
	break;  
      case 'n':
	num_cache_lines = atoi(optarg); 
	orig_num_lines = num_cache_lines;
	IF_EXIT(num_cache_lines < 1, "Number of cache lines must be >= 1\n")
	break;  
      case 'l':
	log2_line_size = atoi(optarg); 
	IF_EXIT(log2_line_size < 0, "Log base 2 of cache line length in bytes must be >= 0\n");
	break;  
      case 'y':
	doprefetch = true; 
	break;
      case 't':
	test_result = true; 
	break;
      case 'o':
	doprint = true; 
	break;
      case 'h':
	printf("Usage: FFT <options>\n\n");
	printf("options:\n");
	printf("  -mM : M = even integer; 2**M total complex data points transformed.\n");
	printf("  -pNumThreads : NumThreads = number of processors; Must be a power of 2.\n");
	printf("  -nN : N = number of cache lines.\n");
	printf("  -lL : L = Log base 2 of cache line length in bytes.\n");
	printf("  -y  : Enable touching the data before operations.\n");
	printf("        It is not useful when measuring the full application or the ROI.\n");
	printf("  -t  : Perform FFT and inverse FFT.  Test output by comparing the\n");
	printf("        integral of the original data to the integral of the data that\n");
	printf("        results from performing the FFT and inverse FFT.\n");
	printf("  -o  : Print out complex data points.\n");
	printf("  -h  : Print out command line options.\n\n");
	printf("Default: FFT -m%1d -p%1d -n%1d -l%1d\n",
	       DEFAULT_M,DEFAULT_THREADS,NUM_CACHE_LINES,LOG2_LINE_SIZE);
	exit(0);
	break;
    }
  }
}


void ParallelMain() {  
  long thread_id = FETCH_ADD(Global->id, 1);
  
  double* upriv = (double* ) G_MALLOC((rootN-1) * size_of_an_element);
  
  for (long i = 0; i < ((rootN-1) << 1); i++) {
    upriv[i] = umain[i];
  }   

  long MyFirst = rootN * thread_id / NumThreads;
  long MyLast = rootN * (thread_id+1) / NumThreads;

  if (doprefetch) {
    TouchArray(x, trans, umain2, upriv, MyFirst, MyLast);
  }
  
  /* perform forward FFT */
  FFT1D(1, M, N, x, trans, upriv, umain2, thread_id, MyFirst, MyLast, pad_length, test_result);

  /* perform backward FFT */
  if (test_result) {
    FFT1D(-1, M, N, x, trans, upriv, umain2, thread_id, MyFirst, MyLast, pad_length, test_result);
  }  
}


double TouchArray(double* x, double* scratch, double* u, double* upriv, long MyFirst, long MyLast) {
  double tot = 0.0;
  /* touch my data */
  for (long j = 0; j < ((rootN-1) << 1); j++) {
    tot += upriv[j];
  }
  /* touch global arrays */
  for (long j = MyFirst; j < MyLast; j++) {
    long k = j * (rootN + pad_length);
    for (long i = 0; i < rootN; i++) {
      tot += x[(k+i) << 1] + x[((k+i) << 1)+1] + 
	     scratch[(k+i) << 1] + scratch[((k+i) << 1)+1] +
	     u[(k+i) << 1] + u[((k+i) << 1)+1];
    }
  }  
  return tot;
}


double CheckSum(double* x) {
  double cks = 0.0;
  for (long j = 0; j < rootN; j++) {
    long k = j * (rootN + pad_length);
    for (long i = 0; i < rootN; i++) {
      cks += x[(k+i) << 1] + x[((k+i) << 1)+1];
    }
  }
  return cks;
}


void InitX(double* x) {
  srand48(0);
  for (long j = 0; j < rootN; j++) {
    long k = j * (rootN + pad_length);
    for (long i = 0; i < rootN; i++) {
      x[(k+i) << 1] = drand48();
      x[((k+i) << 1)+1] = drand48();
    }
  }
}


void InitU(long N, double* u) {
  for (long q = 0; 1<<q<N; q++) {  
    long n1 = 1<<q;
    long base = n1-1;
    for (long j = 0; j < n1; j++) {
      if (base+j > rootN-1) return;
      u[(base+j) << 1] = cos(2.0*PI*j/(n1 << 1));
      u[((base+j) << 1)+1] = -sin(2.0*PI*j/(n1 << 1));
    }
  }
}


void InitU2(long N, double* u, long n1) {
  for (long j = 0; j < n1; j++) {  
    long k = j*(rootN+pad_length);
    for (long i = 0; i < n1; i++) {  
      u[(k+i) << 1] = cos(2.0*PI*i*j/(N));
      u[((k+i) << 1)+1] = -sin(2.0*PI*i*j/(N));
    }
  }
}


long BitReverse(long M, long k) {
  long j = 0;
  long tmp = k;
  for (long i = 0; i < M; i++) {
    j = (j << 1) + (tmp&0x1);
    tmp = tmp>>1;
  }
  return j;
}


void FFT1D(long direction, long M, long N, double* x, double* scratch, double* upriv, double* umain2, long thread_id, long MyFirst, long MyLast, long pad_length, bool test_result) {
  long m1 = M >> 1;
  long n1 = 1 << m1;

  if(is_output){
    printf("iter_num = %lu\n", MyLast - MyFirst);
  }

  /* transpose from x into scratch */
  Transpose(n1, x, scratch, thread_id, MyFirst, MyLast, pad_length);

  /* do n1 1D FFTs on columns */
  for (long j = MyFirst; j < MyLast; j++) {
    FFT1DOnce(direction, m1, n1, upriv, &scratch[(j << 1)*(n1+pad_length)]);
    TwiddleOneCol(direction, n1, j, umain2, &scratch[(j << 1)*(n1+pad_length)], pad_length);
  }  
  
  BARRIER(Global->barrier, NumThreads);

  /* transpose */
  Transpose(n1, scratch, x, thread_id, MyFirst, MyLast, pad_length);

  /* do n1 1D FFTs on columns again */
  for (long j = MyFirst; j < MyLast; j++) {
    FFT1DOnce(direction, m1, n1, upriv, &x[(j << 1)*(n1+pad_length)]);
    if (direction == -1)
      Scale(n1, N, &x[(j << 1)*(n1+pad_length)]);
  }
  
  BARRIER(Global->barrier, NumThreads);

  /* transpose back */
  Transpose(n1, x, scratch, thread_id, MyFirst, MyLast, pad_length);

  /* copy columns from scratch to x */
  if ((test_result) || (doprint)) {
    // All threads has finished, and the data can be copied to the output
    BARRIER(Global->barrier, NumThreads);
    for (long j = MyFirst; j < MyLast; j++) {
      CopyColumn(n1, &scratch[(j << 1)*(n1+pad_length)], &x[(j << 1)*(n1+pad_length)]); 
    }  
  }

  // Data is copied and all threads are done
  BARRIER(Global->barrier, NumThreads);
}


void TwiddleOneCol(long direction, long n1, long j, double* u, double* x, long pad_length) {
  for (long i = 0; i < n1; i++) {
    double omega_r = u[(j*(n1+pad_length)+i) << 1];
    double omega_c = direction*u[((j*(n1+pad_length)+i) << 1)+1];  
    double x_r = x[i << 1]; 
    double x_c = x[(i << 1)+1];
    x[i << 1] = omega_r*x_r - omega_c*x_c;
    x[(i << 1)+1] = omega_r*x_c + omega_c*x_r;
  }
}


void Scale(long n1, long N, double* x) {
  for (long i = 0; i < n1; i++) {
    x[i << 1] /= N;
    x[(i << 1)+1] /= N;
  }
}


void Transpose(long n1, double* src, double* dest, long thread_id, long MyFirst, long MyLast, long pad_length) {
  long iter_num = 0;
  long firstfirst = MyFirst;  
  long row_count = n1 / NumThreads;
  long n1p = n1 + pad_length;
  long blksize = MyLast - MyFirst;
  long numblks = (blksize << 1)/num_cache_lines;
  
  if (numblks * num_cache_lines != (blksize << 1)) {
    numblks ++;
  }
  blksize = blksize / numblks;

  for (long l = thread_id + 1; l < NumThreads; l++) {
    for (long k = 0; k < numblks; k++) {
      long v_off = l * row_count + k * blksize;
      for (long m = 0; m < numblks; m++) {
	long h_off = firstfirst + m * blksize;
        for (long i = 0; i < blksize; i++) {
	  long v = v_off + i; 
          for (long j = 0; j < blksize; j++) {
            iter_num ++;
	    long h = h_off + j;
            dest[(h*n1p+v) << 1] = src[(v*n1p+h) << 1];
            dest[((h*n1p+v) << 1)+1] = src[((v*n1p+h) << 1)+1];
          }
        }	
      }
    }
  }

  if(is_output){
    printf("Transpose: iter_num = %lu\n", iter_num);
  }

  for (long l = 0; l < thread_id; l++) {
    for (long k = 0; k < numblks; k++) {
      long v_off = l * row_count + k * blksize;
      for (long m = 0; m < numblks; m++) {
	long h_off = firstfirst + m * blksize;
        for (long i = 0; i < blksize; i++) {
	  long v = v_off + i;
          for (long j = 0; j < blksize; j++) {
            long h = h_off + j;
            dest[(h*n1p+v) << 1] = src[(v*n1p+h) << 1];
            dest[((h*n1p+v) << 1)+1] = src[((v*n1p+h) << 1)+1];
          }
        }
      }
    }
  }

  for (long k = 0; k < numblks; k++) {
    long v_off = thread_id * row_count + k * blksize;
    for (long m = 0; m < numblks; m++) {
      long h_off = firstfirst + m * blksize;
      for (long i = 0; i < blksize; i++) {
        long v = v_off + i;
        for (long j = 0; j < blksize; j++) {
          long h = h_off + j;
          dest[(h*n1p+v) << 1] = src[(v*n1p+h) << 1];
          dest[((h*n1p+v) << 1)+1] = src[((v*n1p+h) << 1)+1];
	}
      }
    }
  }
}


void CopyColumn(long n1, double* src, double* dest) {
  for (long i = 0; i < n1; i++) {
    dest[i << 1] = src[i << 1];
    dest[(i << 1)+1] = src[(i << 1)+1];
  }
}


void Reverse(long N, long M, double* x) {
  for (long k = 0; k < N; k++) {
    long j = BitReverse(M, k);
    if (j > k) {
      SWAP_VALS(x[j << 1], x[k << 1]);
      SWAP_VALS(x[(j << 1)+1], x[(k << 1)+1]);
    }
  }
}


void FFT1DOnce(long direction, long M, long N, double* u, double* x) {
  long iter_num = 0;

  Reverse(N, M, x);

  for (long q = 1; q <= M; q++) {
    long L = 1<<q;
    long r = N/L;
    long Lstar = L/2;
    double* u1 = &u[(Lstar-1) << 1];
    for (long k = 0; k < r; k++) {
      double* x1 = &x[(k*L) << 1];
      double* x2 = &x[(k*L+Lstar) << 1];
      for (long j = 0; j < Lstar; j++) {
        iter_num ++;
	double omega_r = u1[j << 1]; 
        double omega_c = direction*u1[(j << 1)+1];
	double x_r = x2[j << 1]; 
        double x_c = x2[(j << 1)+1];
	double tau_r = omega_r*x_r - omega_c*x_c;
	double tau_c = omega_r*x_c + omega_c*x_r;
	x_r = x1[j << 1]; 
        x_c = x1[(j << 1)+1];
	x2[j << 1] = x_r - tau_r;
	x2[(j << 1)+1] = x_c - tau_c;
	x1[j << 1] = x_r + tau_r;
	x1[(j << 1)+1] = x_c + tau_c;
      }
    }
  }
  if(is_output){
      printf("FFt1DOnce: iter_num = %lu\n", iter_num);
      is_output = false;
  }
}


void PrintArray(long N, double* x) {
  for (long i = 0; i < rootN; i++) {
    long k = i * (rootN+pad_length);
    for (long j = 0; j < rootN; j++) {
      printf(" %4.2f %4.2f", x[(k+j) << 1], x[((k+j) << 1)+1]);
      if (i*rootN+j != N-1) {
        printf(",");
      }
      if ((i*rootN+j+1) % 8 == 0) {
        printf("\n");
      }
    }
  }
  printf("\n");
  printf("\n");
}


void printerr(char *str) {
  fprintf(stderr,"ERROR: %s\n",str);
}

