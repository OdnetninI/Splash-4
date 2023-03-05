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
/*  This version contains two dimensional arrays in which the first         */
/*  dimension is the block to be operated on, and the second contains       */
/*  all data points in that block.  In this manner, all data points in      */
/*  a block (which are operated on by the same processor) are allocated     */
/*  contiguously and locally, and false sharing is eliminated.              */
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

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdbool.h>

#include "../common/common.h"

MAIN_ENV();

#define MAXRAND                         32767.0
#define DEFAULT_N                         512
#define DEFAULT_THREADS                     1
#define DEFAULT_B                          16
#define min(a,b) ((a) < (b) ? (a) : (b))

#define ERROR_ALLOWED                0.000001

struct GlobalMemory {
  long id;
  BARDEC(barrier)
} *Global;

long n = DEFAULT_N;          /* The size of the matrix */
long NumThreads = DEFAULT_THREADS;          /* Number of processors */
long block_size = DEFAULT_B; /* Block dimension */
long nblocks;                /* Number of blocks in each dimension */
double **a;                  /* a = lu; l and u both placed back in a */
double *rhs;
long *proc_bytes;            /* Bytes to malloc per processor to hold blocks of A*/
double **last_malloc;        /* Starting point of last block of A */

bool test_result = false;        /* Test result of factorization? */
bool doprint = false;            /* Print out matrix values? */
bool doprefetch = false;

/*
  Forward declarations
*/
void ArgumentParser(int argc, char* argv[]);
void ParallelMain(void);
void OneSolve(long n, long block_size, long thread_id);
void lu0(double *a, long n, long stride);
void bdiv(double *a, double *diag, long stride_a, long stride_diag, long dimi, long dimk);
void bmodd(double *a, double *c, long dimi, long dimj, long stride_a, long stride_c);
void bmod(double *a, double *b, double *c, long dimi, long dimj, long dimk, long stridea, long strideb, long stridec);
void daxpy(double *a, double *b, long n, double alpha);
long BlockOwner(long I, long J);
long BlockOwnerColumn(long I, long J);
long BlockOwnerRow(long I, long J);
void lu(long n, long bs, long thread_id);
void InitA(double *rhs);
double TouchA(long bs, long thread_id);
void PrintA(void);
void CheckResult(long n, double **a, double *rhs);
void printerr(char *str);

/*
  Main
 */
int main(int argc, char *argv[]) {
  ArgumentParser(argc, argv);

  MAIN_INITENV(,150000000)

  printf("\n");
  printf("Blocked Dense LU Factorization\n");
  printf("     %ld by %ld Matrix\n", n, n);
  printf("     %ld Processors\n", NumThreads);
  printf("     %ld by %ld Element Blocks\n", block_size, block_size);
  printf("\n");
  printf("\n");

  long num_rows = (long) sqrt((double) NumThreads);
  for (;;) {
    long num_cols = NumThreads/num_rows;
    if (num_rows * num_cols == NumThreads)
      break;
    num_rows--;
  }
  
  nblocks = n/block_size;
  if (block_size * nblocks != n) {
    nblocks++;
  }

  long edge = n%block_size;
  if (edge == 0) {
    edge = block_size;
  }

  proc_bytes = (long *) malloc(NumThreads*sizeof(long));
  if (proc_bytes == NULL) {
	  fprintf(stderr,"Could not malloc memory for proc_bytes.\n");
	  exit(-1);
  }
  last_malloc = (double **) G_MALLOC(NumThreads*sizeof(double *));
  if (last_malloc == NULL) {
	  fprintf(stderr,"Could not malloc memory for last_malloc.\n");
	  exit(-1);
  }
  
  for (long i = 0; i < NumThreads; i++) {
    proc_bytes[i] = 0;
    last_malloc[i] = NULL;
  }
  
  for (long i = 0; i < nblocks; i++) {
    for (long j = 0; j < nblocks; j++) {
      long proc_num = BlockOwner(i,j);
      long size = block_size*block_size;
      if ((i == nblocks-1) && (j == nblocks-1)) {
        size = edge*edge;
      } else if ((i == nblocks-1) || (j == nblocks-1)) {
        size = edge*block_size;
      }
      proc_bytes[proc_num] += size * sizeof(double);
    }
  }
  
  for (long i = 0; i < NumThreads; i++) {
    last_malloc[i] = (double *) G_MALLOC(proc_bytes[i] + PAGE_SIZE);
    if (last_malloc[i] == NULL) {
      fprintf(stderr,"Could not malloc memory blocks for proc %ld\n",i);
      exit(-1);
    } 
    last_malloc[i] = (double *) (((unsigned long) last_malloc[i]) + PAGE_SIZE -
                     ((unsigned long) last_malloc[i]) % PAGE_SIZE);

/* Note that this causes all blocks to start out page-aligned, and that
   for block sizes that exceed cache line size, blocks start at cache-line
   aligned addresses as well.  This reduces false sharing */

  }
  a = (double **) G_MALLOC(nblocks*nblocks*sizeof(double *));
  if (a == NULL) {
    printerr("Could not malloc memory for a\n");
    exit(-1);
  }
  
  for (long i = 0; i < nblocks; i++) {
    for (long j = 0; j < nblocks; j++) {
      long proc_num = BlockOwner(i,j);
      long size = block_size*block_size;
      a[i+j*nblocks] = last_malloc[proc_num];
      if ((i == nblocks-1) && (j == nblocks-1)) {
        size = edge*edge;
      } else if ((i == nblocks-1) || (j == nblocks-1)) {
        size = edge*block_size;
      }
      last_malloc[proc_num] += size;
    }
  }

  rhs = (double *) G_MALLOC(n*sizeof(double));
  Global = (struct GlobalMemory *) G_MALLOC(sizeof(struct GlobalMemory));

/* POSSIBLE ENHANCEMENT:  Here is where one might distribute the a[i]
   blocks across physically distributed memories as desired.

   One way to do this is as follows:

   for (i=0;i<nblocks;i++) {
     for (j=0;j<nblocks;j++) {
       proc_num = BlockOwner(i,j);
       if ((i == nblocks-1) && (j == nblocks-1)) {
         size = edge*edge;
       } else if ((i == nblocks-1) || (j == nblocks-1)) {
         size = edge*block_size;
       } else {
         size = block_size*block_size;
       }

       Place all addresses x such that 
       (&(a[i+j*nblocks][0]) <= x < &(a[i+j*nblocks][size-1])) 
       on node proc_num
     }
   }
*/

  BARINIT(Global->barrier, NumThreads);
  
  Global->id = 0;

  InitA(rhs);
  if (doprint) {
    printf("Matrix before decomposition:\n");
    PrintA();
  }

  SPLASH4_ROI_BEGIN();
  CREATE(ParallelMain, NumThreads);
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
    case 'p': NumThreads = atoi(optarg); break;
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
              printf("Default: LU -n%1d -p%1d -b%1d\n",
		     DEFAULT_N,DEFAULT_THREADS,DEFAULT_B);
              exit(0);
              break;
    }
  }
}


void ParallelMain() {
  long thread_id = FETCH_ADD(Global->id, 1);
  OneSolve(n, block_size, thread_id);
}


void OneSolve(long n, long block_size, long thread_id) {
  if (doprefetch) {
    TouchA(block_size, thread_id);
  }

  /* POSSIBLE ENHANCEMENT:  Here is where one might reset the
     statistics that one is measuring about the parallel execution */

  lu(n, block_size, thread_id);

  BARRIER(Global->barrier, NumThreads);
}


void lu0(double *a, long n, long stride) {
  for (long k = 0; k < n; k++) {
    /* modify subsequent columns */
    for (long j = k + 1; j < n; j++) {
      a[k+j*stride] /= a[k+k*stride];
      double alpha = -a[k+j*stride];
      long length = n - k - 1;
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
  for (long k = 0; k < dimi; k++) {
    for (long j = 0; j < dimj; j++) {
      c[k+j*stride_c] /= a[k+k*stride_a];
      double alpha = -c[k+j*stride_c];
      long length = dimi - k - 1;
      daxpy(&c[k+1+j*stride_c], &a[k+1+k*stride_a], dimi-k-1, alpha);
    }
  }
}


void bmod(double *a, double *b, double *c, long dimi, long dimj, long dimk, long stridea, long strideb, long stridec) {
  for (long k = 0; k < dimk; k++) {
    for (long j = 0; j < dimj; j++) {
      double alpha = -b[k+j*strideb]; 
      daxpy(&c[j*stridec], &a[k*stridea], dimi, alpha);
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

void lu(long n, long bs, long thread_id) {
  for (long k = 0, K = 0; k < n; k += bs, K++) {
    long kl = k + bs;
    long strK = bs;
    if (kl > n) {
      kl = n;
      strK = kl - k;
    }
    /* factor diagonal block */
    if (BlockOwner(K, K) == thread_id) {
      double* A = a[K+K*nblocks]; 
      lu0(A, strK, strK);
    }

    BARRIER(Global->barrier, NumThreads);

    /* divide column k by diagonal block */
    double* D = a[K+K*nblocks];
    for (long i = kl, I = K + 1; i < n; i += bs, I++) {
      if (BlockOwnerColumn(I, K) == thread_id) {  /* parcel out blocks */
	long il = i + bs;
	long strI = bs;
	if (il > n) {
	  il = n;
          strI = il - i;
        }
	double* A = a[I+K*nblocks]; 
	bdiv(A, D, strI, strK, strI, strK);  
      }
    }
    /* modify row k by diagonal block */
    for (long j = kl, J = K + 1; j < n; j += bs, J++) {
      if (BlockOwnerRow(K, J) == thread_id) {  /* parcel out blocks */
	long jl = j+bs;
	long  strJ = bs;
	if (jl > n) {
	  jl = n;
          strJ = jl - j;
        }
        double* A = a[K+J*nblocks];
	bmodd(D, A, strK, strJ, strK, strK);
      }
    }

    BARRIER(Global->barrier, NumThreads);

    /* modify subsequent block columns */
    for (long i = kl, I = K + 1; i < n; i += bs, I++) {
      long il = i+bs;
      long strI = bs;
      if (il > n) {
	il = n;
        strI = il - i;
      }
      double* A = a[I+K*nblocks]; 
      for (long j = kl, J = K + 1; j < n; j += bs, J++) {
	long jl = j + bs;
	long strJ = bs;
	if (jl > n) {
	  jl = n;
          strJ= jl - j;
        }
	if (BlockOwner(I, J) == thread_id) {  /* parcel out blocks */
	  double* B = a[K+J*nblocks]; 
	  double* C = a[I+J*nblocks];
	  bmod(A, B, C, strI, strJ, strK, strI, strK, strI);
	}
      }
    }
  }
}


void InitA(double *rhs) {
  srand48((long) 1);
  long edge = n%block_size;
  for (long j = 0; j < n; j++) {
    for (long i = 0; i < n; i++) {
      long ibs = block_size;
      long jbs = block_size;
      long skip = block_size;
      if ((n - i) <= edge) {
	ibs = n-edge;
	skip = edge;
      }
      if ((n - j) <= edge) {
	jbs = n-edge;
      }
      long ii = (i/block_size) + (j/block_size)*nblocks;
      long jj = (i%ibs)+(j%jbs)*skip;
      a[ii][jj] = ((double) lrand48())/MAXRAND;
      if (i == j) {
	a[ii][jj] *= 10;
      }
    }
  }

  for (long j = 0; j < n; j++) {
    rhs[j] = 0.0;
  }
  
  for (long j = 0; j < n; j++) {
    for (long i = 0; i < n; i++) {
      long ibs = block_size;
      long jbs = block_size;
      long skip = block_size;
      if ((n - i) <= edge) {
	ibs = edge;
	ibs = n-edge;
	skip = edge;
      }
      if ((n - j) <= edge) {
	jbs = edge;
	jbs = n-edge;
      }
      long ii = (i/block_size) + (j/block_size)*nblocks;
      long jj = (i%ibs)+(j%jbs)*skip;
      rhs[i] += a[ii][jj];
    }
  }
}


double TouchA(long bs, long thread_id) {
  double tot = 0.0;
  
  /* touch my portion of A[] */
  for (long J = 0; J < nblocks; J++) {
    for (long I = 0; I < nblocks; I++) {
      if (BlockOwner(I, J) == thread_id) {
	long jbs = bs;
	long ibs = bs;
	if (J == nblocks-1) {
	  jbs = n%bs;
	  if (jbs == 0) {
	    jbs = bs;
          }
	}
	if (I == nblocks-1) {
	  ibs = n%bs;
	  if (ibs == 0) {
	    ibs = bs;
          }
	}
	for (long j = 0; j < jbs; j++) {
	  for (long i = 0; i < ibs; i++) {
	    tot += a[I+J*nblocks][i+j*ibs];
          }
	}
      }
    }
  } 
  return tot;
}


void PrintA() {
  long edge = n%block_size;
  for (long i = 0; i < n; i++) {
    for (long j = 0; j < n; j++) {
      long ibs = block_size;
      long skip = block_size;
      long jbs = block_size;
      if ((n - i) <= edge) {
	ibs = n-edge;
        skip = edge;
      }
      if ((n - j) <= edge) {
	jbs = n-edge;
      }
      long ii = (i/block_size) + (j/block_size)*nblocks;
      long jj = (i%ibs)+(j%jbs)*skip;
      printf("%8.1f ", a[ii][jj]);   
    }
    printf("\n");
  }
  fflush(stdout);
}


void CheckResult(long n, double **a, double *rhs) {
  long edge = n % block_size;
  double* y = (double *) G_MALLOC(n*sizeof(double));  
  
  for (long j = 0; j < n; j++) {
    y[j] = rhs[j];
  }
  
  for (long j = 0; j < n; j++) {
    long jbs = block_size;
    long skip = block_size;
    if ((n - j) <= edge) {
      jbs = n-edge;
      skip = edge;
    }
    long ii = (j/block_size) + (j/block_size)*nblocks;
    long jj = (j%jbs)+(j%jbs)*skip;

    y[j] = y[j] / a[ii][jj];
    
    for (long i = j + 1; i < n; i++) {
      long ibs = block_size;
      long skip = block_size;
      if ((n - i) <= edge) {
        ibs = n-edge;
        skip = edge;
      }
      long iii = (i/block_size) + (j/block_size)*nblocks;
      long jjj = (i%ibs)+(j%jbs)*skip;

      y[i] -= a[iii][jjj] * y[j];
    }
  }

  for (long j = n - 1; j >= 0; j--) {
    for (long i = 0; i < j; i++) {
      long ibs = block_size;
      long jbs = block_size;
      long skip = block_size;
      if ((n - i) <= edge) {
        ibs = n-edge;
        skip = edge;
      }
      if ((n - j) <= edge) {
        jbs = n-edge;
      }
      long ii = (i/block_size) + (j/block_size)*nblocks;
      long jj = (i%ibs)+(j%jbs)*skip;
      y[i] -= a[ii][jj] * y[j];
    }
  }

  bool bogus = false;
  double max_diff = 0.0;
  for (long j = 0; j < n; j++) {
    double diff = y[j] - 1.0;
    if (fabs(diff) > ERROR_ALLOWED) {
      bogus = true;
      max_diff = diff;
    }
  }
  if (bogus) printf("TEST FAILED: (%.5f diff)\n", max_diff);
  else printf("TEST PASSED (%.5f diff)\n", max_diff);
  free(y);
}


void printerr(char *str) {
  fprintf(stderr,"ERROR: %s\n",str);
}

