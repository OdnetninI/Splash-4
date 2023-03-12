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

  /*   ***************************************************************

          f i r s t     p h a s e   (of timestep calculation)

     ***************************************************************/
  
  initializeWithBorders(ga[procid], gp[procid].neighbors, firstrow, lastrow, firstcol, lastcol, 0.0D);
  initializeWithBorders(gb[procid], gp[procid].neighbors, firstrow, lastrow, firstcol, lastcol, 0.0D);

  /* put the laplacian of psi{1,3} in work1{1,2}
   note that psi(i,j,2) represents the psi3 array in
   the original equations  */

  for(long psiindex = 0; psiindex <= 1; psiindex++) {
    double** t2a = (double **) work1[procid][psiindex];
    if ((gp[procid].neighbors[UP] == -1) && (gp[procid].neighbors[LEFT] == -1)) {
      t2a[0][0] = 0;
    }
    if ((gp[procid].neighbors[DOWN] == -1) && (gp[procid].neighbors[LEFT] == -1)) {
      t2a[im-1][0] = 0;
    }
    if ((gp[procid].neighbors[UP] == -1) && (gp[procid].neighbors[RIGHT] == -1)) {
      t2a[0][jm-1] = 0;
    }
    if ((gp[procid].neighbors[DOWN] == -1) && (gp[procid].neighbors[RIGHT] == -1)) {
      t2a[im-1][jm-1] = 0;
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
    double** t2a = (double **) work7[procid][psiindex];
    if ((gp[procid].neighbors[UP] == -1) && (gp[procid].neighbors[LEFT] == -1)) {
      t2a[0][0] = 0;
    }
    if ((gp[procid].neighbors[DOWN] == -1) && (gp[procid].neighbors[LEFT] == -1)) {
      t2a[im-1][0] = 0;
    }
    if ((gp[procid].neighbors[UP] == -1) && (gp[procid].neighbors[RIGHT] == -1)) {
      t2a[0][jm-1] = 0;
    }
    if ((gp[procid].neighbors[DOWN] == -1) && (gp[procid].neighbors[RIGHT] == -1)) {
      t2a[im-1][jm-1] = 0;
    }
    laplacalc(procid,psim,work7,psiindex,firstrow,lastrow,firstcol,lastcol);
  }

  /* to the values of the work1{1,2} arrays obtained from the
   laplacians of psi{1,2} in the previous phase, add to the
   elements of every column the corresponding value in the
   one-dimenional f array  */

  for(long psiindex = 0; psiindex <= 1; psiindex++) {
    double** t2a = (double **) work1[procid][psiindex];
    if ((gp[procid].neighbors[UP] == -1) && (gp[procid].neighbors[LEFT] == -1)) {
      t2a[0][0] = t2a[0][0] + f[0];
    }
    if ((gp[procid].neighbors[DOWN] == -1) && (gp[procid].neighbors[LEFT] == -1)) {
      t2a[im-1][0] = t2a[im-1][0] + f[0];
    }
    if ((gp[procid].neighbors[UP] == -1) && (gp[procid].neighbors[RIGHT] == -1)) {
      t2a[0][jm-1] = t2a[0][jm-1] + f[jmx[numlev-1]-1];
    }
    if ((gp[procid].neighbors[DOWN] == -1) && (gp[procid].neighbors[RIGHT] == -1)) {
      t2a[im-1][jm-1]=t2a[im-1][jm-1] + f[jmx[numlev-1]-1];
    }
    if (gp[procid].neighbors[UP] == -1) {
      for(long j = firstcol; j <= lastcol; j++) {
        t2a[0][j] = t2a[0][j] + f[j+j_off];
      }
    }
    if (gp[procid].neighbors[DOWN] == -1) {
      for(long j = firstcol; j <= lastcol; j++) {
        t2a[im-1][j] = t2a[im-1][j] + f[j+j_off];
      }
    }
    if (gp[procid].neighbors[LEFT] == -1) {
      for(long j = firstrow; j <= lastrow; j++) {
        t2a[j][0] = t2a[j][0] + f[j+i_off];
      }
    }
    if (gp[procid].neighbors[RIGHT] == -1) {
      for(long j = firstrow; j <= lastrow; j++) {
        t2a[j][jm-1] = t2a[j][jm-1] + f[j+i_off];
      }
    }
    for(long i = firstrow; i <= lastrow; i++) {
      double* t1a = (double *) t2a[i];
      for(long iindex = firstcol; iindex <= lastcol; iindex++) {
        t1a[iindex]=t1a[iindex] + f[iindex+j_off];
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

  double** t2a = (double **) ga[procid];
  double** t2b = (double **) gb[procid];
  double** t2c = (double **) work5[procid][0];
  double** t2d = (double **) work5[procid][1];
  double** t2e = (double **) work7[procid][0];
  double** t2f = (double **) work7[procid][1];
  double** t2g = (double **) work6[procid];
  double** t2h = (double **) tauz[procid];
  if ((gp[procid].neighbors[UP] == -1) && (gp[procid].neighbors[LEFT] == -1)) {
    t2a[0][0] = t2c[0][0]-t2d[0][0] + eig2*t2g[0][0]+h1inv*t2h[0][0] + lf*t2e[0][0]-lf*t2f[0][0];
    t2b[0][0] = hh1*t2c[0][0]+hh3*t2d[0][0] + hinv*t2h[0][0]+lf*hh1*t2e[0][0] + lf*hh3*t2f[0][0];
  }
  if ((gp[procid].neighbors[DOWN] == -1) && (gp[procid].neighbors[LEFT] == -1)) {
    t2a[im-1][0] = t2c[im-1][0]-t2d[im-1][0] + eig2*t2g[im-1][0] + h1inv*t2h[im-1][0] + lf*t2e[im-1][0] - lf*t2f[im-1][0];
    t2b[im-1][0] = hh1*t2c[im-1][0] + hh3*t2d[im-1][0] + hinv*t2h[im-1][0] + lf*hh1*t2e[im-1][0] + lf*hh3*t2f[im-1][0];
  }
  if ((gp[procid].neighbors[UP] == -1) && (gp[procid].neighbors[RIGHT] == -1)) {
    t2a[0][jm-1] = t2c[0][jm-1]-t2d[0][jm-1]+ eig2*t2g[0][jm-1]+h1inv*t2h[0][jm-1] + lf*t2e[0][jm-1]-lf*t2f[0][jm-1];
    t2b[0][jm-1] = hh1*t2c[0][jm-1] + hh3*t2d[0][jm-1]+hinv*t2h[0][jm-1] + lf*hh1*t2e[0][jm-1]+lf*hh3*t2f[0][jm-1];
  }
  if ((gp[procid].neighbors[DOWN] == -1) && (gp[procid].neighbors[RIGHT] == -1)) {
    t2a[im-1][jm-1] = t2c[im-1][jm-1] - t2d[im-1][jm-1]+eig2*t2g[im-1][jm-1] + h1inv*t2h[im-1][jm-1]+lf*t2e[im-1][jm-1] - lf*t2f[im-1][jm-1];
    t2b[im-1][jm-1] = hh1*t2c[im-1][jm-1] + hh3*t2d[im-1][jm-1]+hinv*t2h[im-1][jm-1] + lf*hh1*t2e[im-1][jm-1] + lf*hh3*t2f[im-1][jm-1];
  }
  if (gp[procid].neighbors[UP] == -1) {
    double* t1a = (double *) t2a[0];
    double* t1b = (double *) t2b[0];
    double* t1c = (double *) t2c[0];
    double* t1d = (double *) t2d[0];
    double* t1e = (double *) t2e[0];
    double* t1f = (double *) t2f[0];
    double* t1g = (double *) t2g[0];
    double* t1h = (double *) t2h[0];
    for(long j = firstcol; j <= lastcol; j++) {
      t1a[j] = t1c[j]-t1d[j] + eig2*t1g[j]+h1inv*t1h[j] + lf*t1e[j]-lf*t1f[j];
      t1b[j] = hh1*t1c[j] + hh3*t1d[j]+hinv*t1h[j] + lf*hh1*t1e[j]+lf*hh3*t1f[j];
    }
  }
  if (gp[procid].neighbors[DOWN] == -1) {
    double* t1a = (double *) t2a[im-1];
    double* t1b = (double *) t2b[im-1];
    double* t1c = (double *) t2c[im-1];
    double* t1d = (double *) t2d[im-1];
    double* t1e = (double *) t2e[im-1];
    double* t1f = (double *) t2f[im-1];
    double* t1g = (double *) t2g[im-1];
    double* t1h = (double *) t2h[im-1];
    for(long j = firstcol; j <= lastcol; j++) {
      t1a[j] = t1c[j] - t1d[j]+eig2*t1g[j] + h1inv*t1h[j]+lf*t1e[j] - lf*t1f[j];
      t1b[j] = hh1*t1c[j] + hh3*t1d[j]+hinv*t1h[j] + lf*hh1*t1e[j]+lf*hh3*t1f[j];
    }
  }
  if (gp[procid].neighbors[LEFT] == -1) {
    for(long j = firstrow; j <= lastrow; j++) {
      t2a[j][0] = t2c[j][0]-t2d[j][0] + eig2*t2g[j][0]+h1inv*t2h[j][0] + lf*t2e[j][0]-lf*t2f[j][0];
      t2b[j][0] = hh1*t2c[j][0] + hh3*t2d[j][0]+hinv*t2h[j][0] + lf*hh1*t2e[j][0]+lf*hh3*t2f[j][0];
    }
  }
  if (gp[procid].neighbors[RIGHT] == -1) {
    for(long j = firstrow; j <= lastrow; j++) {
      t2a[j][jm-1] = t2c[j][jm-1] - t2d[j][jm-1]+eig2*t2g[j][jm-1] + h1inv*t2h[j][jm-1]+lf*t2e[j][jm-1] - lf*t2f[j][jm-1];
      t2b[j][jm-1] = hh1*t2c[j][jm-1] + hh3*t2d[j][jm-1]+hinv*t2h[j][jm-1] + lf*hh1*t2e[j][jm-1]+lf*hh3*t2f[j][jm-1];
    }
  }

  for(long i = firstrow; i <= lastrow; i++) {
    double* t1a = (double *) t2a[i];
    double* t1b = (double *) t2b[i];
    double* t1c = (double *) t2c[i];
    double* t1d = (double *) t2d[i];
    double* t1e = (double *) t2e[i];
    double* t1f = (double *) t2f[i];
    double* t1g = (double *) t2g[i];
    double* t1h = (double *) t2h[i];
    for(long iindex = firstcol; iindex <= lastcol; iindex++) { 
      t1a[iindex] = t1c[iindex] - t1d[iindex]+eig2*t1g[iindex] + h1inv*t1h[iindex]+lf*t1e[iindex] - lf*t1f[iindex];
      t1b[iindex] = hh1*t1c[iindex] + hh3*t1d[iindex]+hinv*t1h[iindex] + lf*hh1*t1e[iindex] + lf*hh3*t1f[iindex];
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

  if (gp[procid].neighbors[UP] == -1) {
    istart = 0;
  }
  if (gp[procid].neighbors[LEFT] == -1) {
    jstart = 0;
  }
  if (gp[procid].neighbors[DOWN] == -1) {
    iend = im-1;
  }
  if (gp[procid].neighbors[RIGHT] == -1) {
    jend = jm-1;
  }
  t2a = (double **) rhs_multi[procid][numlev-1];
  t2b = (double **) ga[procid];
  t2c = (double **) oldga[procid];
  t2d = (double **) q_multi[procid][numlev-1];
  for(long i = istart; i <= iend; i++) {
    double* t1a = (double *) t2a[i];
    double* t1b = (double *) t2b[i];
    for(long j = jstart; j <= jend; j++) {
      t1a[j] = t1b[j] * ressqr;
    }
  }

  if (gp[procid].neighbors[UP] == -1) {
    double* t1d = (double *) t2d[0];
    double* t1b = (double *) t2b[0];
    for(long j = jstart; j <= jend; j++) {
      t1d[j] = t1b[j];
    }
  }
  if (gp[procid].neighbors[DOWN] == -1) {
    double* t1d = (double *) t2d[im-1];
    double* t1b = (double *) t2b[im-1];
    for(long j = jstart; j <= jend; j++) {
      t1d[j] = t1b[j];
    }
  }
  if (gp[procid].neighbors[LEFT] == -1) {
    for(long i = istart; i <= iend; i++) {
      t2d[i][0] = t2b[i][0];
    }
  }
  if (gp[procid].neighbors[RIGHT] == -1) {
    for(long i = istart; i <= iend; i++) {
      t2d[i][jm-1] = t2b[i][jm-1];
    }
  }

  for(long i = ist; i <= ien; i++) {
    double* t1d = (double *) t2d[i];
    double* t1c = (double *) t2c[i];
    for(long j = jst; j <= jen; j++) {
      t1d[j] = t1c[j];
    }
  }

  multig(procid);

  /* the shared sum variable psiai is initialized to 0 at
   every time-step  */

  if (procid == MASTER) {
    global->psiai=0.0;
  }

  /*  copy the solution for use as initial guess in next time-step  */

  for(long i = istart; i <= iend; i++) {
    double* t1b = (double *) t2b[i];
    double* t1c = (double *) t2c[i];
    double* t1d = (double *) t2d[i];
    for(long j = jstart; j <= jend; j++) {
      t1b[j] = t1d[j];
      t1c[j] = t1d[j];
    }
  }

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

  t2a = (double **) rhs_multi[procid][numlev-1];
  t2b = (double **) gb[procid];
  t2c = (double **) oldgb[procid];
  t2d = (double **) q_multi[procid][numlev-1];
  for(long i = istart; i <= iend; i++) {
    double* t1a = (double *) t2a[i];
    double* t1b = (double *) t2b[i];
    for(long j = jstart; j <= jend; j++) {
      t1a[j] = t1b[j] * ressqr;
    }
  }
  if (gp[procid].neighbors[UP] == -1) {
    double* t1d = (double *) t2d[0];
    double* t1b = (double *) t2b[0];
    for(long j = jstart; j <= jend; j++) {
      t1d[j] = t1b[j];
    }
  }
  if (gp[procid].neighbors[DOWN] == -1) {
    double* t1d = (double *) t2d[im-1];
    double* t1b = (double *) t2b[im-1];
    for(long j = jstart; j <= jend; j++) {
      t1d[j] = t1b[j];
    }
  }
  if (gp[procid].neighbors[LEFT] == -1) {
    for(long i = istart; i <= iend; i++) {
      t2d[i][0] = t2b[i][0];
    }
  }
  if (gp[procid].neighbors[RIGHT] == -1) {
    for(long i = istart; i <= iend; i++) {
      t2d[i][jm-1] = t2b[i][jm-1];
    }
  }

  for(long i = ist; i <= ien; i++) {
    double* t1d = (double *) t2d[i];
    double* t1c = (double *) t2c[i];
    for(long j = jst; j <= jen; j++) {
      t1d[j] = t1c[j];
    }
  }

  multig(procid);

  for(long i = istart; i <= iend; i++) {
    double* t1b = (double *) t2b[i];
    double* t1c = (double *) t2c[i];
    double* t1d = (double *) t2d[i];
    for(long j = jstart; j <= jend; j++) {
      t1b[j] = t1d[j];
      t1c[j] = t1d[j];
    }
  }

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
