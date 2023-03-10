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
/*  Integer radix sort of non-negative integers.                            */
/*                                                                          */
/*  Command line options:                                                   */
/*                                                                          */
/*  -pNumThreads : pNumThreads = number of processors.                      */
/*  -rR : R = radix for sorting.  Must be power of 2.                       */
/*  -nN : N = number of keys to sort.                                       */
/*  -mM : M = maximum key value.  Integer keys k will be generated such     */
/*        that 0 <= k <= M.                                                 */
/*  -t  : Check to make sure all keys are sorted correctly.                 */
/*  -o  : Print out sorted keys.                                            */
/*  -h  : Print out command line options.                                   */
/*                                                                          */
/*  Default: RADIX -p1 -n262144 -r1024 -m524288                             */
/*                                                                          */
/*  Note: This version works under both the FORK and SPROC models           */
/*                                                                          */
/****************************************************************************/

#include "../common/common.h"

#define DEFAULT_THREADS              1
#define DEFAULT_N               262144
#define DEFAULT_R                 1024 
#define DEFAULT_M             67108864
#define MAX_PROCESSORS            1024
#define RADIX_S                8388608.0e0
#define RADIX           70368744177664.0e0
#define SEED                 314159265.0e0
#define RATIO               1220703125.0e0
#define PAGE_SIZE                 4096
#define PAGE_MASK     (~(PAGE_SIZE-1))
#define MAX_RADIX                 4096

MAIN_ENV();

struct prefix_node {
   long densities[MAX_RADIX];
   long ranks[MAX_RADIX];
   PAUSEDEC(done)
   char pad[PAGE_SIZE];
};

struct global_memory {
   long Index;                             /* process ID */
   long final;
   BARDEC(barrier_rank)                   /* for ranking process */
   BARDEC(barrier_key)                    /* for key sorting process */
   struct prefix_node prefix_tree[2 * MAX_PROCESSORS];
} *global;

struct global_private {
  char pad[PAGE_SIZE];
  long *rank_ff;         /* overall processor ranks */
} gp[MAX_PROCESSORS];

long *key[2];            /* sort from one index into the other */
long **rank_me;          /* individual processor ranks */
long *key_partition;     /* keys a processor works on */
long *rank_partition;    /* ranks a processor works on */

long number_of_processors = DEFAULT_THREADS;
long max_num_digits;
long radix = DEFAULT_R;
long num_keys = DEFAULT_N;
long max_key = DEFAULT_M;
long log2_radix;
long log2_keys;

bool test_result = false;
bool doprint = false;

/*
  Forward declarations
*/
void ArgumentParser(int argc, char* argv[]);
void slave_sort(void);
double product_mod_46(double t1, double t2);
double ran_num_init(unsigned long k, double b, double t);
long get_max_digits(long max_key);
long get_log2_radix(long rad);
long get_log2_keys(long num_keys);
void printerr(char *s);
void init(long key_start, long key_stop, long from);
void test_sort(long final);
void printout(void);

int main(int argc, char *argv[]) {
  ArgumentParser(argc, argv);
  MAIN_INITENV(,80000000);
  
  log2_radix = LOG2(radix); 
  log2_keys = LOG2(num_keys);
  
  global = (struct global_memory *) G_MALLOC(sizeof(struct global_memory));
   
  key[0] = (long *) G_MALLOC(num_keys*sizeof(long));
  key[1] = (long *) G_MALLOC(num_keys*sizeof(long));
  key_partition = (long *) G_MALLOC((number_of_processors+1)*sizeof(long));
  rank_partition = (long *) G_MALLOC((number_of_processors+1)*sizeof(long));;
  long size = number_of_processors*(radix*sizeof(long)+sizeof(long *));
  rank_me = (long **) G_MALLOC(size);

  long** temp = rank_me;
  long** temp2 = temp + number_of_processors;
  long* a = (long *) temp2;

  for (long i = 0; i < number_of_processors; i++) {
    *temp = (long *) a;
    temp++;
    a += radix;
  }
  for (long i = 0; i < number_of_processors; i++) {
    gp[i].rank_ff = (long *) G_MALLOC(radix*sizeof(long)+PAGE_SIZE);
  }

  BARINIT(global->barrier_rank, number_of_processors)
  BARINIT(global->barrier_key, number_of_processors)

  for (long i = 0; i < number_of_processors << 1; i++) {
    PAUSEINIT(global->prefix_tree[i].done);
  }

  global->Index = 0;
  max_num_digits = get_max_digits(max_key);
  printf("\n");
  printf("Integer Radix Sort\n");
  printf("     %ld Keys\n",num_keys);
  printf("     %ld Processors\n",number_of_processors);
  printf("     Radix = %ld\n",radix);
  printf("     Max key = %ld\n",max_key);
  printf("\n");

  long quotient = num_keys / number_of_processors;
  long remainder = num_keys % number_of_processors;
  long sum_i = 0;
  long sum_f = 0;
  long p = 0;
  while (sum_i < num_keys) {
    key_partition[p] = sum_i;
    p++;
    sum_i = sum_i + quotient;
    sum_f = sum_f + remainder;
    sum_i = sum_i + sum_f / number_of_processors;
    sum_f = sum_f % number_of_processors;
  }
  key_partition[p] = num_keys;

  quotient = radix / number_of_processors;
  remainder = radix % number_of_processors;
  sum_i = 0;
  sum_f = 0;
  p = 0;
  while (sum_i < radix) {
    rank_partition[p] = sum_i;
    p++;
    sum_i = sum_i + quotient;
    sum_f = sum_f + remainder;
    sum_i = sum_i + sum_f / number_of_processors;
    sum_f = sum_f % number_of_processors;
  }
  rank_partition[p] = radix;

  /* POSSIBLE ENHANCEMENT:  Here is where one might distribute the key,
   rank_me, rank, and gp data structures across physically 
   distributed memories as desired. 
   
   One way to place data is as follows:  
      
   for (i=0;i<number_of_processors;i++) {
     Place all addresses x such that:
       &(key[0][key_partition[i]]) <= x < &(key[0][key_partition[i+1]])
            on node i
       &(key[1][key_partition[i]]) <= x < &(key[1][key_partition[i+1]])
            on node i
       &(rank_me[i][0]) <= x < &(rank_me[i][radix-1]) on node i
       &(gp[i]) <= x < &(gp[i+1]) on node i
       &(gp[i].rank_ff[0]) <= x < &(gp[i].rank_ff[radix]) on node i
   }
   start_p = 0;
   i = 0;

   for (toffset = 0; toffset < number_of_processors; toffset ++) {
     offset = toffset;
     level = number_of_processors >> 1;
     base = number_of_processors;
     while ((offset & 0x1) != 0) {
       offset >>= 1;
       index = base + offset;
       Place all addresses x such that:
         &(global->prefix_tree[index]) <= x < 
              &(global->prefix_tree[index + 1]) on node toffset
       base += level;
       level >>= 1;
     }  
   }  */

  /* Fill the random-number array. */
  SPLASH4_ROI_BEGIN();
  CREATE(slave_sort, number_of_processors);
  WAIT_FOR_END(number_of_processors);
  SPLASH4_ROI_END();

  if (doprint) {
    printout();
  }
  if (test_result) {
    test_sort(global->final);  
  }
  
  MAIN_END();
}

void ArgumentParser(int argc, char* argv[]) {
  long c;
  extern char *optarg;

  while ((c = getopt(argc, argv, "p:r:n:m:stoh")) != -1) {
  switch(c) {
    case 'p': number_of_processors = atoi(optarg);
              IF_EXIT(number_of_processors < 1, "NumThreads must be >= 1\n");
              IF_EXIT(number_of_processors > MAX_PROCESSORS, "Maximum processors (MAX_PROCESSORS) exceeded\n");
              break;
    case 'r': radix = atoi(optarg);
              IF_EXIT(radix < 1, "Radix must be a power of 2 greater than 0\n");
              IF_EXIT(radix != (1 << LOG2(radix)), "Radix must be a power of 2\n");
              break;
    case 'n': num_keys = atoi(optarg);
              IF_EXIT(num_keys < 1, "Number of keys must be >= 1\n");
              break;
    case 'm': max_key = atoi(optarg);
              IF_EXIT(max_key < 1, "Maximum key must be >= 1\n");
              break;
    case 't': test_result = true; break;
    case 'o': doprint = true; break;
    case 'h': printf("Usage: RADIX <options>\n\n");
              printf("   -pPpNumThreads : pNumThreads = number of processors.\n");
              printf("   -rR : R = radix for sorting.  Must be power of 2.\n");
              printf("   -nN : N = number of keys to sort.\n");
              printf("   -mM : M = maximum key value.  Integer keys k will be generated such\n");
              printf("         that 0 <= k <= M.\n");
              printf("   -t  : Check to make sure all keys are sorted correctly.\n");
              printf("   -o  : Print out sorted keys.\n");
              printf("   -h  : Print out command line options.\n\n");
              printf("Default: RADIX -p%1d -n%1d -r%1d -m%1d\n",
                      DEFAULT_THREADS,DEFAULT_N,DEFAULT_R,DEFAULT_M);
  exit(0);
    }
  }
}

void slave_sort() {
  long MyNum = FETCH_ADD(global->Index, 1);

  /* POSSIBLE ENHANCEMENT:  Here is where one might pin processes to
  processors to avoid migration */

  long* key_density = (long *) G_MALLOC(radix*sizeof(long));

   /* Fill the random-number array. */

  long key_start = key_partition[MyNum];
  long key_stop = key_partition[MyNum + 1];
  long rank_start = rank_partition[MyNum];
  long rank_stop = rank_partition[MyNum + 1];
  if (rank_stop == radix) {
    rank_stop--;
  }

  long from = 0;
  long to = 1;
  init(key_start,key_stop,from);

  BARRIER(global->barrier_key, number_of_processors) 

  /* POSSIBLE ENHANCEMENT:  Here is where one might reset the
  statistics that one is measuring about the parallel execution */

  BARRIER(global->barrier_key, number_of_processors) 

  /* Do 1 iteration per digit.  */

  long* rank_me_mynum = rank_me[MyNum];
  long* rank_ff_mynum = gp[MyNum].rank_ff;
  for (long loopnum = 0; loopnum < max_num_digits; loopnum++) {
    long shiftnum = (loopnum * log2_radix);
    long bb = (radix-1) << shiftnum;

    /* generate histograms based on one digit */

    for (long i = 0; i < radix; i++) {
      rank_me_mynum[i] = 0;
    }

    long* key_from = (long *) key[from];
    long* key_to = (long *) key[to];
    for (long i = key_start; i < key_stop; i++) {
      long my_key = key_from[i] & bb;
      my_key = my_key >> shiftnum;  
      rank_me_mynum[my_key]++;
    }

    key_density[0] = rank_me_mynum[0]; 
    for (long i = 1; i < radix; i++) {
      key_density[i] = key_density[i-1] + rank_me_mynum[i];  
    }

    BARRIER(global->barrier_key, number_of_processors) 

    struct prefix_node* n = &(global->prefix_tree[MyNum]);
    for (long i = 0; i < radix; i++) {
      n->densities[i] = key_density[i];
      n->ranks[i] = rank_me_mynum[i];
    }

    long offset = MyNum;
    long level = number_of_processors >> 1;
    long base = number_of_processors;
    if ((MyNum & 0x1) == 0) {
      SETPAUSE(global->prefix_tree[base + (offset >> 1)].done);
    }

    while ((offset & 0x1) != 0) {
      offset >>= 1;
      struct prefix_node* r = n;
      struct prefix_node* l = n - 1;
      long index = base + offset;
      n = &(global->prefix_tree[index]);
      WAITPAUSE(n->done);
      CLEARPAUSE(n->done);
      if (offset != (level - 1)) {
        for (long i = 0; i < radix; i++) {
          n->densities[i] = r->densities[i] + l->densities[i];
          n->ranks[i] = r->ranks[i] + l->ranks[i];
        }
      } else {
        for (long i = 0; i < radix; i++) {
          n->densities[i] = r->densities[i] + l->densities[i];
        }
      }
      base += level;
      level >>= 1;
      if ((offset & 0x1) == 0) {
        SETPAUSE(global->prefix_tree[base + (offset >> 1)].done);
      }
    }
    
    BARRIER(global->barrier_rank, number_of_processors);

    struct prefix_node* my_node;
    struct prefix_node* their_node;

    if (MyNum != (number_of_processors - 1)) {
      offset = MyNum;
      level = number_of_processors;
      base = 0;
      while ((offset & 0x1) != 0) {
        offset >>= 1;
        base += level;
        level >>= 1;
      }
      my_node = &(global->prefix_tree[base + offset]);
      offset >>= 1;
      base += level;
      level >>= 1;
      while ((offset & 0x1) != 0) {
        offset >>= 1;
        base += level;
        level >>= 1;
      }
      their_node = &(global->prefix_tree[base + offset]);
      WAITPAUSE(my_node->done);
      CLEARPAUSE(my_node->done);
      for (long i = 0; i < radix; i++) {
        my_node->densities[i] = their_node->densities[i];
      }
    } else {
      my_node = &(global->prefix_tree[(2 * number_of_processors) - 2]);
    }

    offset = MyNum;
    level = number_of_processors;
    base = 0;
    while ((offset & 0x1) != 0) {
      SETPAUSE(global->prefix_tree[base + offset - 1].done);
      offset >>= 1;
      base += level;
      level >>= 1;
    }

    offset = MyNum;
    level = number_of_processors;
    base = 0;
    for(long i = 0; i < radix; i++) {
      rank_ff_mynum[i] = 0;
    }

    while (offset != 0) {
      if ((offset & 0x1) != 0) {
        /* Add ranks of node to your left at this level */
        struct prefix_node* l = &(global->prefix_tree[base + offset - 1]);
        for (long i = 0; i < radix; i++) {
          rank_ff_mynum[i] += l->ranks[i];
        }
      }
      base += level;
      level >>= 1;
      offset >>= 1;
    }

    for (long i = 1; i < radix; i++) {
      rank_ff_mynum[i] += my_node->densities[i - 1];
    }
    
    BARRIER(global->barrier_rank, number_of_processors) 

    /* put it in order according to this digit */

    for (long i = key_start; i < key_stop; i++) {  
      long this_key = key_from[i] & bb;
      this_key = this_key >> shiftnum;  
      long tmp = rank_ff_mynum[this_key];
      key_to[tmp] = key_from[i];
      rank_ff_mynum[this_key]++;
    }   /*  i */  

    if (loopnum != max_num_digits-1) {
      from = from ^ 0x1;
      to = to ^ 0x1;
    }
    
    BARRIER(global->barrier_rank, number_of_processors) 
  } /* for */

  if (MyNum == 0) global->final = to;
}

/*
 * product_mod_46() returns the product (mod 2^46) of t1 and t2.
 */
double product_mod_46(double t1, double t2) {			
   double a1 = (double)((long)(t1 / RADIX_S));    /* Decompose the arguments.  */
   double a2 = t1 - a1 * RADIX_S;
   double b1 = (double)((long)(t2 / RADIX_S));
   double b2 = t2 - b1 * RADIX_S;
   t1 = a1 * b2 + a2 * b1;      /* Multiply the arguments.  */
   t2 = (double)((long)(t1 / RADIX_S));
   t2 = t1 - t2 * RADIX_S;
   t1 = t2 * RADIX_S + a2 * b2;
   t2 = (double)((long)(t1 / RADIX));

   return (t1 - t2 * RADIX);    /* Return the product.  */
}

/*
 * finds the (k)th random number, given the seed, b, and the ratio, t.
 */
double ran_num_init(unsigned long k, double b, double t) {
  while (k != 0) {  /* while() is executed m times such that 2^m > k.  */
    unsigned long j = k >> 1;
    if ((j << 1) != k) {
      b = product_mod_46(b, t);
    }
    t = product_mod_46(t, t);
    k = j;
  }
  return b;
}

long get_max_digits(long max_key) {
  long temp = 1;
  long key_val = max_key;
  while (true) {
    key_val = key_val / radix;
    if (key_val == 0) break;
    temp++;
  }
  return temp;
}

long get_log2_radix(long rad) {
  if (rad == (1 << LOG2(rad))) return LOG2(rad);
  fprintf(stderr,"ERROR: Radix %ld not a power of 2\n", rad);
  exit(-1);
}

long get_log2_keys(long num_keys) {
  if (num_keys == (1 << LOG2(num_keys))) return LOG2(num_keys);
  fprintf(stderr,"ERROR: Number of keys %ld not a power of 2\n", num_keys);
  exit(-1);
}

void printerr(char *s) {
  fprintf(stderr,"ERROR: %s\n",s);
}

void init(long key_start, long key_stop, long from) {
   double ran_num = ran_num_init((key_start << 2) + 1, SEED, RATIO);
   double sum = ran_num / RADIX;
   long* key_from = (long *) key[from];
   for (long i = key_start; i < key_stop; i++) {
      ran_num = product_mod_46(ran_num, RATIO);
      sum = sum + ran_num / RADIX;
      ran_num = product_mod_46(ran_num, RATIO);
      sum = sum + ran_num / RADIX;
      ran_num = product_mod_46(ran_num, RATIO);
      sum = sum + ran_num / RADIX;
      key_from[i] = (long) ((sum / 4.0) *  max_key);
      ran_num = product_mod_46(ran_num, RATIO);
      sum = ran_num / RADIX;
   }
}

void test_sort(long final) {
  long mistake = 0;

  printf("\n");
  printf("                  TESTING RESULTS\n");
  long* key_final = key[final];
  for (long i = 0; i < num_keys-1; i++) {
    if (key_final[i] > key_final[i + 1]) {
      fprintf(stderr,"error with key %ld, value %ld %ld \n", i,key_final[i],key_final[i + 1]);
      mistake++;
    }
  }

  if (mistake > 0) printf("FAILED: %ld keys out of place.\n", mistake);
  else printf("PASSED: All keys in place.\n");
  printf("\n");
}

void printout() {
  long* key_final = (long *) key[global->final];
  printf("\n");
  printf("                 SORTED KEY VALUES\n");
  printf("%8ld ",key_final[0]);
  for (long i = 0; i < num_keys-1; i++) {
    printf("%8ld ",key_final[i+1]);
    if ((i+2)%5 == 0) printf("\n");
  }
  printf("\n");
}

