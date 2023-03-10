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
/*  Parallel dense blocked LU factorization (no pivoting)                   */
/*                                                                          */
/*  This version contains one dimensional arrays in which the matrix        */
/*  to be factored is stored.                                               */
/*                                                                          */
/*  Command line options:                                                   */
/*                                                                          */
/*  -nN : Decompose NxN matrix.                                             */
/*  -pNumThreads : NumThreads = number of processors.                       */
/*  -bB : Use a block size of B. BxB elements should fit in cache for       */
/*        good performance. Small block sizes (B=8, B=16) work well.        */
/*  -y  : Enable touching the data before operations.                       */
/*        It is not useful when measuring the full application or the ROI.  */
/*  -t  : Test output.                                                      */
/*  -o  : Print out matrix values.                                          */
/*  -h  : Print out command line options.                                   */
/*                                                                          */
/****************************************************************************/

#include "../common/common.h"

MAIN_ENV();

#define MAXRAND                         32767.0
#define DEFAULT_N                         128
#define DEFAULT_THREADS                     1
#define DEFAULT_B                          16
#define min(a,b) ((a) < (b) ? (a) : (b))

#define ERROR_ALLOWED                0.000001

struct GlobalMemory {
  long id;
  BARDEC(start)
} *Global;

long n = DEFAULT_N;          /* The size of the matrix */
long NumThreads = DEFAULT_THREADS;          /* Number of processors */
long block_size = DEFAULT_B; /* Block dimension */
long nblocks;                /* Number of blocks in each dimension */
long num_rows;               /* Number of processors per row of processor grid */
long num_cols;               /* Number of processors per col of processor grid */
double *a;                   /* a = lu; l and u both placed back in a */
double *rhs;
long *proc_bytes;            /* Bytes to malloc per processor to hold blocks of A*/

bool test_result = false;        /* Test result of factorization? */
bool doprint = false;            /* Print out matrix values? */
bool doprefetch = false;     

/*
  Forward declarations
*/
void ArgumentParser(int argc, char* argv[]);
void SlaveStart(void);
void OneSolve(long n, long block_size, long MyNum);
void lu0(double *a, long n, long stride);
void bdiv(double *a, double *diag, long stride_a, long stride_diag, long dimi, long dimk);
void bmodd(double *a, double *c, long dimi, long dimj, long stride_a, long stride_c);
void bmod(double *a, double *b, double *c, long dimi, long dimj, long dimk, long stride);
void daxpy(double *a, double *b, long n, double alpha);
long BlockOwner(long I, long J);
long BlockOwnerColumn(long I, long J);
long BlockOwnerRow(long I, long J);
void lu(long n, long bs, long MyNum);
void InitA(double *rhs);
double TouchA(long bs, long MyNum);
void PrintA(void);
void CheckResult(long n, double *a, double *rhs);
void printerr(char *s);

/*
  Main
 */
int main(int argc, char *argv[]) {
  ArgumentParser(argc, argv);

  MAIN_INITENV(,150000000);

  printf("\n");
  printf("Blocked Dense LU Factorization\n");
  printf("     %ld by %ld Matrix\n",n,n);
  printf("     %ld Processors\n",NumThreads);
  printf("     %ld by %ld Element Blocks\n",block_size,block_size);
  printf("\n");
  printf("\n");

  long num_rows = (long) sqrt((double) NumThreads);
  for (;;) {
    long num_cols = NumThreads/num_rows;
    if (num_rows*num_cols == NumThreads)
      break;
    num_rows--;
  }
  long nblocks = n/block_size;
  if (block_size * nblocks != n) {
    nblocks++;
  }

  a = (double *) G_MALLOC(n*n*sizeof(double));
  rhs = (double *) G_MALLOC(n*sizeof(double));

  Global = (struct GlobalMemory *) G_MALLOC(sizeof(struct GlobalMemory));

  /* POSSIBLE ENHANCEMENT:  Here is where one might distribute the a
  matrix data across physically distributed memories in a 
  round-robin fashion as desired. */

  BARINIT(Global->start, NumThreads);
  
  Global->id = 0;

  InitA(rhs);
  if (doprint) {
    printf("Matrix before decomposition:\n");
    PrintA();
  }

  SPLASH4_ROI_BEGIN();
  CREATE(SlaveStart, NumThreads);
  WAIT_FOR_END(NumThreads);

  SPLASH4_ROI_END();

  if (doprint) {
    printf("\nMatrix after decomposition:\n");
    PrintA();
  }

  if (test_result) {
    printf("                             TESTING RESULTS\n");
    CheckResult(n, a, rhs);
  }

  MAIN_END();
}

void ArgumentParser(int argc, char* argv[]) {
  long ch;
  extern char *optarg;
  while ((ch = getopt(argc, argv, "n:p:b:cstoh")) != -1) {
    switch(ch) {
      case 'n': n = atoi(optarg); break;
      case 'p': NumThreads = atoi(optarg);
          IF_EXIT(NumThreads < 1, "NumThreads must be >= 1\n");
          break;
      case 'b': block_size = atoi(optarg); break;
      case 'y': doprefetch = true; break;
      case 't': test_result = true; break;
      case 'o': doprint = true; break;
      case 'h': printf("Usage: LU <options>\n\n");
          printf("options:\n");
          printf("  -nN : Decompose NxN matrix.\n");
          printf("  -pNumThreads : NumThreads = number of processors.\n");
          printf("  -bB : Use a block size of B. BxB elements should fit in cache for \n");
          printf("        good performance. Small block sizes (B=8, B=16) work well.\n");
          printf("  -y  : Enable touching the data before operations.\n");
          printf("        It is not useful when measuring the full application or the ROI.\n");
          printf("  -c  : Copy non-locally allocated blocks to local memory before use.\n");
          printf("  -t  : Test output.\n");
          printf("  -o  : Print out matrix values.\n");
          printf("  -h  : Print out command line options.\n\n");
          printf("Default: LU -n%1d -p%1d -b%1d\n", DEFAULT_N,DEFAULT_THREADS,DEFAULT_B);
          exit(0);
          break;
    }
  }
}

void SlaveStart() {
  long MyNum = FETCH_ADD(Global->id, 1);
  OneSolve(n, block_size, MyNum);
}

void OneSolve(long n, long block_size, long MyNum) {
  if (doprefetch) {
    TouchA(block_size, MyNum);
  }

  /* POSSIBLE ENHANCEMENT:  Here is where one might reset the
     statistics that one is measuring about the parallel execution */

  lu(n, block_size, MyNum);

  BARRIER(Global->start, NumThreads);
}


void lu0(double *a, long n, long stride) {
  for (long k = 0; k < n; k++) {
    /* modify subsequent columns */
    for (long j = k + 1; j < n; j++) {
      a[k+j*stride] /= a[k+k*stride];
      double alpha = -a[k+j*stride];
      long length = n-k-1;
      daxpy(&a[k+1+j*stride], &a[k+1+k*stride], n-k-1, alpha);
    }
  }
}


void bdiv(double *a, double *diag, long stride_a, long stride_diag, long dimi, long dimk) {
  for (long k = 0; k < dimk; k++) {
    for (long j = k + 1; j < dimk; j++) {
      double alpha = -diag[k+j*stride_diag];
      daxpy(&a[j*stride_a], &a[k*stride_a], dimi, alpha);
    }
  }
}


void bmodd(double *a, double *c, long dimi, long dimj, long stride_a, long stride_c) {
  for (long k = 0; k < dimi; k++)
    for (long j = 0; j < dimj; j++) {
      c[k+j*stride_c] /= a[k+k*stride_a];
      double alpha = -c[k+j*stride_c];
      long length = dimi - k - 1;
      daxpy(&c[k+1+j*stride_c], &a[k+1+k*stride_a], dimi-k-1, alpha);
    }
}


void bmod(double *a, double *b, double *c, long dimi, long dimj, long dimk, long stride) {
  for (long k = 0; k < dimk; k++) {
    for (long j = 0; j < dimj; j++) {
      double alpha = -b[k+j*stride];
      daxpy(&c[j*stride], &a[k*stride], dimi, alpha);
    }
  }
}


void daxpy(double *a, double *b, long n, double alpha) {
  for (long i = 0; i < n; i++) {
    a[i] += alpha*b[i];
  }
}


long BlockOwner(long I, long J) {
	return((I + J*nblocks) % NumThreads);
}

long BlockOwnerColumn(long I, long J) {
	return(I % NumThreads);
}

long BlockOwnerRow(long I, long J) {
	return(((J % NumThreads) + (NumThreads / 2)) % NumThreads);
}

void lu(long n, long bs, long MyNum) {
  long strI = n;
  for (long k = 0, K = 0; k < n; k += bs, K++) {
    long kl = k+bs; 
    if (kl>n) {
      kl = n;
    }

    double* A;

    /* factor diagonal block */
    if (BlockOwner(K, K) == MyNum) {
      A = &(a[k+k*n]); 
      lu0(A, kl-k, strI);
    }

    BARRIER(Global->start, NumThreads);

    /* divide column k by diagonal block */
    double* D = &(a[k+k*n]);
    for (long i = kl, I = K + 1; i < n; i += bs, I++) {
      if (BlockOwner/*Column*/(I, K) == MyNum) {  /* parcel out blocks */
        long il = i + bs;
        if (il > n) {
          il = n;
        }
        A = &(a[i+k*n]);
        bdiv(A, D, strI, n, il-i, kl-k);
      }
    }
    /* modify row k by diagonal block */
    for (long j = kl, J = K + 1; j < n; j += bs, J++) {
      if (BlockOwner/*Row*/(K, J) == MyNum) {  /* parcel out blocks */
        long jl = j+bs;
        if (jl > n) {
          jl = n;
        }
        A = &(a[k+j*n]);
        bmodd(D, A, kl-k, jl-j, n, strI);
      }
    }
    
    BARRIER(Global->start, NumThreads);

    /* modify subsequent block columns */
    for (long i = kl, I = K + 1; i < n; i += bs, I++) {
      long il = i + bs;
      if (il > n) {
        il = n;
      }
      A = &(a[i+k*n]);
      for (long j = kl, J = K + 1; j < n; j += bs, J++) {
        long jl = j + bs;
        if (jl > n) {
          jl = n;
        }
        if (BlockOwner(I, J) == MyNum) {  /* parcel out blocks */
          double* B = &(a[k+j*n]);
          double* C = &(a[i+j*n]);
          bmod(A, B, C, il-i, jl-j, kl-k, n);
        }
      }
    }
  }
}


void InitA(double *rhs) {
  srand48((long) 1);
  for (long j = 0; j < n; j++) {
    for (long i = 0; i < n; i++) {
      a[i+j*n] = (double) lrand48()/MAXRAND;
      if (i == j) {
	a[i+j*n] *= 10;
      }
    }
  }

  for (long j = 0; j < n; j++) {
    rhs[j] = 0.0;
  }
  for (long j = 0; j < n; j++) {
    for (long i = 0; i < n; i++) {
      rhs[i] += a[i+j*n];
    }
  }
}


double TouchA(long bs, long MyNum) {
  double tot = 0.0;
  for (long J = 0; J*bs < n; J++) {
    for (long I = 0; I*bs < n; I++) {
      if (BlockOwner(I, J) == MyNum) {
        for (long j = J*bs; j < (J+1)*bs && j < n; j++) {
          for (long i = I*bs; i < (I+1)*bs && i < n; i++) {
            tot += a[i+j*n];
          }
        }
      }
    }
  }
  return tot;
}


void PrintA() {
  for (long i = 0; i < n; i++) {
    for (long j = 0; j < n; j++) {
      printf("%8.1f ", a[i+j*n]);
    }
    printf("\n");
  }
}


void CheckResult(long n, double *a, double *rhs) {
  double* y = (double *) G_MALLOC(n*sizeof(double));
  
  for (long j = 0; j < n; j++) {
    y[j] = rhs[j];
  }
  for (long j = 0; j < n; j++) {
    y[j] = y[j]/a[j+j*n];
    for (long i = j+1; i < n; i++) {
      y[i] -= a[i+j*n]*y[j];
    }
  }

  for (long j = n-1; j >= 0; j--) {
    for (long i = 0; i < j; i++) {
      y[i] -= a[i+j*n]*y[j];
    }
  }

  bool bogus = false;
  double max_diff = 0.0;
  for (long j=0; j < n; j++) {
    double diff = y[j] - 1.0;
    if (fabs(diff) > 0.00001) {
      bogus = true;
      max_diff = diff;
    }
  }
  
  if (bogus) printf("TEST FAILED: (%.5f diff)\n", max_diff);
  else printf("TEST PASSED (%.5f diff)\n", max_diff);
  free(y);
}


void printerr(char *s) {
  fprintf(stderr,"ERROR: %s\n",s);
}

