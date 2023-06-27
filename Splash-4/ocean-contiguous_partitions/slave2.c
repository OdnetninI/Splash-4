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

/*    ****************
      subroutine slave2
      ****************  */
#include "../common/common.h"

EXTERN_ENV();

#include "decs.h"

void slave2(long procid, long firstrow, long lastrow, long numrows, long firstcol, long lastcol, long numcols) {
  double ressqr = lev_res[numlev-1] * lev_res[numlev-1];
  long i_off = gp[procid].rownum*numrows;
  long j_off = gp[procid].colnum*numcols;

  long up_neighbor = gp[procid].neighbors[UP];
  long down_neighbor = gp[procid].neighbors[DOWN];
  long left_neighbor = gp[procid].neighbors[LEFT];
  long right_neighbor = gp[procid].neighbors[RIGHT];

  /*   ***************************************************************

          f i r s t     p h a s e   (of timestep calculation)

     ***************************************************************/
  
  initializeWithBorders(ga[procid], gp[procid].neighbors, firstrow, lastrow, firstcol, lastcol, 0.0D);
  initializeWithBorders(gb[procid], gp[procid].neighbors, firstrow, lastrow, firstcol, lastcol, 0.0D);

  /* put the laplacian of psi{1,3} in work1{1,2}
   note that psi(i,j,2) represents the psi3 array in
   the original equations  */

  for(long psiindex = 0; psiindex <= 1; psiindex++) {
    double** work1_local = (double **) work1[procid][psiindex];
    if ((up_neighbor == -1) && (left_neighbor == -1)) {
      work1_local[0][0] = 0;
    }
    if ((down_neighbor == -1) && (left_neighbor == -1)) {
      work1_local[im-1][0] = 0;
    }
    if ((up_neighbor == -1) && (right_neighbor == -1)) {
      work1_local[0][jm-1] = 0;
    }
    if ((down_neighbor == -1) && (right_neighbor == -1)) {
      work1_local[im-1][jm-1] = 0;
    }
    laplacalc(procid,psi,work1,psiindex,firstrow,lastrow,firstcol,lastcol);
  }

  /* set values of work2 array to psi1 - psi3   */
  subValues(work2[procid], psi[procid][0], psi[procid][1], gp[procid].neighbors, firstrow, lastrow, firstcol, lastcol);

  /* set values of work3 array to h3/h * psi1 + h1/h * psi3  */
  double hh1 = h1/h;
  double hh3 = h3/h;
  addValuesWith2Weight(work3[procid], hh3, work3[procid], hh1, psi[procid][1], gp[procid].neighbors, firstrow, lastrow, firstcol, lastcol);

  /* set values of temparray{1,3} to psim{1,3}  */

  for(long psiindex = 0; psiindex <= 1; psiindex++) {
    copyValues(temparray[procid][psiindex], psi[procid][psiindex], gp[procid].neighbors, firstrow, lastrow, firstcol, lastcol);
  }

  /*     *******************************************************

              s e c o n d   p h a s e

       *******************************************************

   set values of psi{1,3} to psim{1,3}   */

  for(long psiindex = 0; psiindex <= 1; psiindex++) {
    copyValues(psi[procid][psiindex], psim[procid][psiindex], gp[procid].neighbors, firstrow, lastrow, firstcol, lastcol);
  }

/* put the laplacian of the psim array
   into the work7 array; first part of a three-laplacian
   calculation to compute the friction terms  */

  for(long psiindex = 0; psiindex <= 1; psiindex++) {
    double** work7_local = (double **) work7[procid][psiindex];
    if ((up_neighbor == -1) && (left_neighbor == -1)) {
      work7_local[0][0] = 0;
    }
    if ((down_neighbor == -1) && (left_neighbor == -1)) {
      work7_local[im-1][0] = 0;
    }
    if ((up_neighbor == -1) && (right_neighbor == -1)) {
      work7_local[0][jm-1] = 0;
    }
    if ((down_neighbor == -1) && (right_neighbor == -1)) {
      work7_local[im-1][jm-1] = 0;
    }
    laplacalc(procid,psim,work7,psiindex,firstrow,lastrow,firstcol,lastcol);
  }

  /* to the values of the work1{1,2} arrays obtained from the
   laplacians of psi{1,2} in the previous phase, add to the
   elements of every column the corresponding value in the
   one-dimenional f array  */

  for(long psiindex = 0; psiindex <= 1; psiindex++) {
    double** work1_local = (double **) work1[procid][psiindex];
    if ((up_neighbor == -1) && (left_neighbor == -1)) {
      work1_local[0][0] = work1_local[0][0] + f[0];
    }
    if ((down_neighbor == -1) && (left_neighbor == -1)) {
      work1_local[im-1][0] = work1_local[im-1][0] + f[0];
    }
    if ((up_neighbor == -1) && (right_neighbor == -1)) {
      work1_local[0][jm-1] = work1_local[0][jm-1] + f[jmx[numlev-1]-1];
    }
    if ((down_neighbor == -1) && (right_neighbor == -1)) {
      work1_local[im-1][jm-1]=work1_local[im-1][jm-1] + f[jmx[numlev-1]-1];
    }
    if (up_neighbor == -1) {
      for(long j = firstcol; j <= lastcol; j++) {
        work1_local[0][j] = work1_local[0][j] + f[j+j_off];
      }
    }
    if (down_neighbor == -1) {
      for(long j = firstcol; j <= lastcol; j++) {
        work1_local[im-1][j] = work1_local[im-1][j] + f[j+j_off];
      }
    }
    if (left_neighbor == -1) {
      for(long j = firstrow; j <= lastrow; j++) {
        work1_local[j][0] = work1_local[j][0] + f[j+i_off];
      }
    }
    if (right_neighbor == -1) {
      for(long j = firstrow; j <= lastrow; j++) {
        work1_local[j][jm-1] = work1_local[j][jm-1] + f[j+i_off];
      }
    }
    for(long i = firstrow; i <= lastrow; i++) {
      double* work1_local_i = (double *) work1_local[i];
      for(long col = firstcol; col <= lastcol; col++) {
        work1_local_i[col] = work1_local_i[col] + f[col+j_off];
      }
    }
  }

  BARRIER(bars->sl_phase_2,nprocs);

  /* 	*******************************************************

                 t h i r d   p h a s e

 	*******************************************************

   put the jacobian of the work1{1,2} and psi{1,3} arrays
   (the latter currently in temparray) in the work5{1,2} arrays  */

  for(long psiindex = 0; psiindex <= 1; psiindex++) {
    jacobcalc2(work1,temparray,work5,psiindex,procid,firstrow,lastrow,firstcol,lastcol);
  }

  /* set values of psim{1,3} to temparray{1,3}  */

  for(long psiindex = 0; psiindex <= 1; psiindex++) {
    copyValues(psim[procid][psiindex], temparray[procid][psiindex], gp[procid].neighbors, firstrow, lastrow, firstcol, lastcol);
  }

  /* put the laplacian of the work7{1,2} arrays in the work4{1,2}
   arrays; second step in the three-laplacian friction calculation  */

  for(long psiindex = 0; psiindex <= 1; psiindex++) {
    laplacalc(procid,work7,work4,psiindex,firstrow,lastrow,firstcol,lastcol);
  }

  BARRIER(bars->sl_phase_3,nprocs);

  /*     *******************************************************

                f o u r t h   p h a s e

       *******************************************************

   put the jacobian of the work2 and work3 arrays in the work6
   array  */

  jacobcalc(work2,work3,work6,procid,firstrow,lastrow,firstcol,lastcol);

  /* put the laplacian of the work4{1,2} arrays in the work7{1,2}
   arrays; third step in the three-laplacian friction calculation  */

  for(long psiindex = 0; psiindex <= 1; psiindex++) {
    laplacalc(procid,work4,work7,psiindex,firstrow,lastrow,firstcol,lastcol);
  }
  
  BARRIER(bars->sl_phase_4,nprocs);

  /*     *******************************************************

                f i f t h   p h a s e

       *******************************************************

   use the values of the work5, work6 and work7 arrays
   computed in the previous time-steps to compute the
   ga and gb arrays   */

  double hinv = 1.0/h;
  double h1inv = 1.0/h1;

  double** ga_local = (double **) ga[procid];
  double** gb_local = (double **) gb[procid];
  double** work5_local_0 = (double **) work5[procid][0];
  double** work5_local_1 = (double **) work5[procid][1];
  double** work7_local_0 = (double **) work7[procid][0];
  double** work7_local_1 = (double **) work7[procid][1];
  double** work6_local = (double **) work6[procid];
  double** tauz_local = (double **) tauz[procid];
  if ((up_neighbor == -1) && (left_neighbor == -1)) {
    ga_local[0][0] = work5_local_0[0][0]-work5_local_1[0][0] + eig2*work6_local[0][0]+h1inv*tauz_local[0][0] + lf*work7_local_0[0][0]-lf*work7_local_1[0][0];
    gb_local[0][0] = hh1*work5_local_0[0][0]+hh3*work5_local_1[0][0] + hinv*tauz_local[0][0]+lf*hh1*work7_local_0[0][0] + lf*hh3*work7_local_1[0][0];
  }
  if ((down_neighbor == -1) && (left_neighbor == -1)) {
    ga_local[im-1][0] = work5_local_0[im-1][0]-work5_local_1[im-1][0] + eig2*work6_local[im-1][0] + h1inv*tauz_local[im-1][0] + lf*work7_local_0[im-1][0] - lf*work7_local_1[im-1][0];
    gb_local[im-1][0] = hh1*work5_local_0[im-1][0] + hh3*work5_local_1[im-1][0] + hinv*tauz_local[im-1][0] + lf*hh1*work7_local_0[im-1][0] + lf*hh3*work7_local_1[im-1][0];
  }
  if ((up_neighbor == -1) && (right_neighbor == -1)) {
    ga_local[0][jm-1] = work5_local_0[0][jm-1]-work5_local_1[0][jm-1]+ eig2*work6_local[0][jm-1]+h1inv*tauz_local[0][jm-1] + lf*work7_local_0[0][jm-1]-lf*work7_local_1[0][jm-1];
    gb_local[0][jm-1] = hh1*work5_local_0[0][jm-1] + hh3*work5_local_1[0][jm-1]+hinv*tauz_local[0][jm-1] + lf*hh1*work7_local_0[0][jm-1]+lf*hh3*work7_local_1[0][jm-1];
  }
  if ((down_neighbor == -1) && (right_neighbor == -1)) {
    ga_local[im-1][jm-1] = work5_local_0[im-1][jm-1] - work5_local_1[im-1][jm-1]+eig2*work6_local[im-1][jm-1] + h1inv*tauz_local[im-1][jm-1]+lf*work7_local_0[im-1][jm-1] - lf*work7_local_1[im-1][jm-1];
    gb_local[im-1][jm-1] = hh1*work5_local_0[im-1][jm-1] + hh3*work5_local_1[im-1][jm-1]+hinv*tauz_local[im-1][jm-1] + lf*hh1*work7_local_0[im-1][jm-1] + lf*hh3*work7_local_1[im-1][jm-1];
  }
  if (up_neighbor == -1) {
    double* ga_local_0 = (double *) ga_local[0];
    double* gb_local_0 = (double *) gb_local[0];
    double* t1c = (double *) work5_local_0[0];
    double* t1d = (double *) work5_local_1[0];
    double* t1e = (double *) work7_local_0[0];
    double* t1f = (double *) work7_local_1[0];
    double* t1g = (double *) work6_local[0];
    double* t1h = (double *) tauz_local[0];
    for(long col = firstcol; col <= lastcol; col++) {
      ga_local_0[col] = t1c[col]-t1d[col] + eig2*t1g[col]+h1inv*t1h[col] + lf*t1e[col]-lf*t1f[col];
      gb_local_0[col] = hh1*t1c[col] + hh3*t1d[col]+hinv*t1h[col] + lf*hh1*t1e[col]+lf*hh3*t1f[col];
    }
  }
  if (down_neighbor == -1) {
    double* ga_local_last = (double *) ga_local[im-1];
    double* gb_local_last = (double *) gb_local[im-1];
    double* t1c = (double *) work5_local_0[im-1];
    double* t1d = (double *) work5_local_1[im-1];
    double* t1e = (double *) work7_local_0[im-1];
    double* t1f = (double *) work7_local_1[im-1];
    double* t1g = (double *) work6_local[im-1];
    double* t1h = (double *) tauz_local[im-1];
    for(long col = firstcol; col <= lastcol; col++) {
      ga_local_last[col] = t1c[col] - t1d[col]+eig2*t1g[col] + h1inv*t1h[col]+lf*t1e[col] - lf*t1f[col];
      gb_local_last[col] = hh1*t1c[col] + hh3*t1d[col]+hinv*t1h[col] + lf*hh1*t1e[col]+lf*hh3*t1f[col];
    }
  }
  if (left_neighbor == -1) {
    for(long row = firstrow; row <= lastrow; row++) {
      ga_local[row][0] = work5_local_0[row][0]-work5_local_1[row][0] + eig2*work6_local[row][0]+h1inv*tauz_local[row][0] + lf*work7_local_0[row][0]-lf*work7_local_1[row][0];
      gb_local[row][0] = hh1*work5_local_0[row][0] + hh3*work5_local_1[row][0]+hinv*tauz_local[row][0] + lf*hh1*work7_local_0[row][0]+lf*hh3*work7_local_1[row][0];
    }
  }
  if (right_neighbor == -1) {
    for(long row = firstrow; row <= lastrow; row++) {
      ga_local[row][jm-1] = work5_local_0[row][jm-1] - work5_local_1[row][jm-1]+eig2*work6_local[row][jm-1] + h1inv*tauz_local[row][jm-1]+lf*work7_local_0[row][jm-1] - lf*work7_local_1[row][jm-1];
      gb_local[row][jm-1] = hh1*work5_local_0[row][jm-1] + hh3*work5_local_1[row][jm-1]+hinv*tauz_local[row][jm-1] + lf*hh1*work7_local_0[row][jm-1]+lf*hh3*work7_local_1[row][jm-1];
    }
  }

  for(long i = firstrow; i <= lastrow; i++) {
    double* ga_local_i = (double *) ga_local[i];
    double* gb_local_i = (double *) gb_local[i];
    double* t1c = (double *) work5_local_0[i];
    double* t1d = (double *) work5_local_1[i];
    double* t1e = (double *) work7_local_0[i];
    double* t1f = (double *) work7_local_1[i];
    double* t1g = (double *) work6_local[i];
    double* t1h = (double *) tauz_local[i];
    for(long col = firstcol; col <= lastcol; col++) { 
      ga_local_i[col] = t1c[col] - t1d[col]+eig2*t1g[col] + h1inv*t1h[col]+lf*t1e[col] - lf*t1f[col];
      gb_local_i[col] = hh1*t1c[col] + hh3*t1d[col]+hinv*t1h[col] + lf*hh1*t1e[col] + lf*hh3*t1f[col];
    }
  }
   
  BARRIER(bars->sl_phase_5,nprocs);

  /*     *******************************************************

               s i x t h   p h a s e

       *******************************************************  */

  long istart = 1;
  long iend = istart + gp[procid].rel_num_y[numlev-1] - 1;
  long jstart = 1;
  long jend = jstart + gp[procid].rel_num_x[numlev-1] - 1;
  long ist = istart;
  long ien = iend;
  long jst = jstart;
  long jen = jend;

  if (up_neighbor == -1) istart = 0;
  if (left_neighbor == -1) jstart = 0;
  if (down_neighbor == -1) iend = im-1;
  if (right_neighbor == -1) jend = jm-1;
  
  double** rhs_multi_local = (double **) rhs_multi[procid][numlev-1];
  double** t2c = (double **) oldga[procid];
  double** t2d = (double **) q_multi[procid][numlev-1];
  Matrix_copy_coef(double, rhs_multi[procid][numlev-1], ga_local, ressqr, istart, iend + 1, jstart, jend + 1);

  if (up_neighbor == -1) { Vector_copy(double, q_multi[procid][numlev-1][0], ga_local[0], jstart, jend + 1); }
  if (down_neighbor == -1) { Vector_copy(double, q_multi[procid][numlev-1][im-1], ga_local[im-1], jstart, jend + 1); }
  if (left_neighbor == -1) { Matrix_copy_col(double, q_multi[procid][numlev-1], ga_local, istart, iend + 1, 0, 0); }
  if (right_neighbor == -1) { Matrix_copy_col(double, q_multi[procid][numlev-1], ga_local, istart, iend + 1, jm-1, jm-1); }
  Matrix_copy(double, q_multi[procid][numlev-1], oldga[procid], ist, ien + 1, jst, jen+1);

  multig(procid);

  /* the shared sum variable psiai is initialized to 0 at
   every time-step  */

  if (procid == MASTER) {
    global->psiai=0.0;
  }

  /*  copy the solution for use as initial guess in next time-step  */
  Matrix_copy_to_two_dest(double, ga_local, oldga[procid], q_multi[procid][numlev-1], istart, iend + 1, jstart, jend + 1);

  BARRIER(bars->sl_phase_6,nprocs);

  /*     *******************************************************

                s e v e n t h   p h a s e

       *******************************************************

   every process computes the running sum for its assigned portion
   in a private variable psiaipriv   */

  double psiaipriv = calculatePsiipriv(ga[procid], gp[procid].neighbors, firstrow, lastrow, firstcol, lastcol);

  /* after computing its private sum, every process adds that to the
   shared running sum psiai  */

  FETCH_ADD_DOUBLE(&global->psiai, psiaipriv);

  BARRIER(bars->sl_phase_7,nprocs);

  /*      *******************************************************

                e i g h t h   p h a s e

        *******************************************************

   augment ga(i,j) with [-psiai/psibi]*psib(i,j) */
  addValuesWithWeight(ga[procid], ga[procid], (-global->psiai)/(global->psibi), psib[procid], gp[procid].neighbors, firstrow, lastrow, firstcol, lastcol);

  t2c = (double **) oldgb[procid];
  t2d = (double **) q_multi[procid][numlev-1];
  Matrix_copy_coef(double, rhs_multi_local, gb_local, ressqr, istart, iend + 1, jstart, jend + 1);

  if (up_neighbor == -1) { Vector_copy(double, q_multi[procid][numlev-1][0], gb_local[0], jstart, jend + 1); }
  if (down_neighbor == -1) { Vector_copy(double, q_multi[procid][numlev-1][im-1], gb_local[im-1], jstart, jend + 1); }
  if (left_neighbor == -1) { Matrix_copy_col(double, q_multi[procid][numlev-1], gb_local, istart, iend + 1, 0, 0); }
  if (right_neighbor == -1) { Matrix_copy_col(double, q_multi[procid][numlev-1], gb_local, istart, iend + 1, jm-1, jm-1); }
  Matrix_copy(double, q_multi[procid][numlev-1], oldgb[procid], ist, ien + 1, jst, jen+1);

  multig(procid);

  Matrix_copy_to_two_dest(double, gb_local, oldgb[procid], q_multi[procid][numlev-1], istart, iend + 1, jstart, jend + 1);

/*      *******************************************************

                n i n t h   p h a s e

        *******************************************************

   put appropriate linear combinations of ga and gb in work2 and work3;
   note that here (as in most cases) the constant multipliers are made
   private variables; the specific order in which things are done is
   chosen in order to hopefully reuse things brought into the cache

   note that here again we choose to have all processes share the work
   on both matrices despite the fact that the work done per element
   is the same, because the operand matrices are the same in both cases */

  addValuesWithWeight(work3[procid], gb[procid], hh3, ga[procid], gp[procid].neighbors, firstrow, lastrow, firstcol, lastcol);
  subValuesWithWeight(work2[procid], gb[procid], hh1, ga[procid], gp[procid].neighbors, firstrow, lastrow, firstcol, lastcol);

  /*      *******************************************************

                t e n t h    p h a s e

        *******************************************************/
  /* update the psi{1,3} matrices by adding 2*dtau*work3 to each */
  addValuesWithWeight(psi[procid][0], psi[procid][0], 2 * dtau, work3[procid], gp[procid].neighbors, firstrow, lastrow, firstcol, lastcol);
  addValuesWithWeight(psi[procid][1], psi[procid][1], 2 * dtau, work2[procid], gp[procid].neighbors, firstrow, lastrow, firstcol, lastcol);

   BARRIER(bars->sl_phase_10,nprocs)   
}
