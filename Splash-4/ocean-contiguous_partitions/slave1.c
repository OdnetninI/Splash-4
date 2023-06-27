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
      subroutine slave
      ****************  */
#include "../common/common.h"

EXTERN_ENV();

#include "decs.h"

void slave() {
  double ressqr = lev_res[numlev-1] * lev_res[numlev-1];
  long procid = FETCH_ADD(global->id, 1);

  long up_neighbor = gp[procid].neighbors[UP];
  long down_neighbor = gp[procid].neighbors[DOWN];
  long left_neighbor = gp[procid].neighbors[LEFT];
  long right_neighbor = gp[procid].neighbors[RIGHT];

  /* POSSIBLE ENHANCEMENT:  Here is where one might pin processes to
   processors to avoid migration. */

  /* POSSIBLE ENHANCEMENT:  Here is where one might distribute
   data structures across physically distributed memories as
   desired. */

  Matrix_multi2_set(double, 0.0, oldga[procid], oldgb[procid], 0, im, 0, jm);

  long firstcol = 1;
  long lastcol = firstcol + gp[procid].rel_num_x[numlev-1] - 1;
  long firstrow = 1;
  long lastrow = firstrow + gp[procid].rel_num_y[numlev-1] - 1;
  long numcols = gp[procid].rel_num_x[numlev-1];
  long numrows = gp[procid].rel_num_y[numlev-1];
  long j_off = gp[procid].colnum*numcols;

  /* every process gets its own copy of the timing variables to avoid
   contention at shared memory locations.  here, these variables
   are initialized.  */

  double ttime = 0.0;
  double dhour = 0.0;
  long nstep = 0 ;
  double day = 0.0;
  double ysca1 = 0.5*ysca;

  if (procid == MASTER) {
    Vector_set(double, (f0 + beta * ((((double) _i) * res) - ysca1)), f, 0, jmx[numlev-1]-1);
  }

  initializeWithBorders(psium[procid], gp[procid].neighbors, firstrow, lastrow, firstcol, lastcol, 0.0D);
  initializeWithBorders(psilm[procid], gp[procid].neighbors, firstrow, lastrow, firstcol, lastcol, 0.0D);
  initializeWithBorders(psib[procid], gp[procid].neighbors, firstrow, lastrow, firstcol, lastcol, 1.0D);

  /* compute psib array (one-time computation) and integrate into psibi */

  long istart = 1;
  long iend = istart + gp[procid].rel_num_y[numlev-1] - 1;
  long jstart = 1;
  long jend = jstart + gp[procid].rel_num_x[numlev-1] - 1;
  long ist = istart;
  long ien = iend;
  long jst = jstart;
  long jen = jend;

  if (up_neighbor == -1)  istart = 0;
  if (left_neighbor == -1) jstart = 0;
  if (down_neighbor == -1)  iend = im-1;
  if (right_neighbor == -1)  jend = jm-1;

  double** rhs_multi_local = (double **) rhs_multi[procid][numlev-1];
  double** psib_local = (double **) psib[procid];
  double** q_multi_local = (double **) q_multi[procid][numlev-1];

  Matrix_copy_coef(double, rhs_multi_local, psib_local, ressqr, istart, iend, jstart, jend);

  if (up_neighbor == -1) { Vector_copy(double, q_multi_local[0], psib_local[0], jstart, jend); }
  else { Vector_copy(double, psib_local[0], psib[up_neighbor][im-2], 1, jm - 1); }

  if (down_neighbor == -1) { Vector_copy(double, q_multi_local[0], psib_local[im-1], jstart, jend); }
  else { Vector_copy(double, psib_local[im-1], psib[down_neighbor][1], 1, jm - 1); }

  if (left_neighbor == -1) { Matrix_copy_col(double, q_multi_local, psib_local, istart, iend, 0, 0); }
  else { Matrix_copy_col(double, psib_local, psib[left_neighbor], 1, im-1, 0, jm-2); }

  if (right_neighbor == -1) { Matrix_copy_col(double, q_multi_local, psib_local, istart, iend, jm-1, jm-1); }
  else { Matrix_copy_col(double, psib_local, psib[right_neighbor], 1, im-1, jm-1, 1); }

  double fac = 1.0 / (4.0 - ressqr*eig2);
  for(long i = ist; i <= ien; i++) {
    double* q_multi_local_i = (double *) q_multi_local[i];
    double* psib_local_i = (double *) psib_local[i];
    double* psib_local_prev = (double *) psib_local[i-1];
    double* psib_local_next = (double *) psib_local[i+1];
    for(long j = jst; j <= jen; j++) {
      q_multi_local_i[j] = fac * (psib_local_next[j]+psib_local_prev[j]+psib_local_i[j+1]+psib_local_i[j-1] - ressqr*psib_local_i[j]);
    }
  }

  multig(procid);

  Matrix_copy(double, q_multi_local, psib_local, istart, iend, jstart, jend);

  /* update the local running sum psibipriv by summing all the resulting
   values in that process's share of the psib matrix   */

  double psibipriv = calculatePsiipriv(psib[procid], gp[procid].neighbors, firstrow, lastrow, firstcol, lastcol);

  /* update the shared variable psibi by summing all the psibiprivs
   of the individual processes into it.  note that this combined
   private and shared sum method avoids accessing the shared
   variable psibi once for every element of the matrix.  */

  FETCH_ADD_DOUBLE(&global->psibi, psibipriv);
   
  /* initialize psim matrices
   if there is more than one process, then split the processes
   between the two psim matrices; otherwise, let the single process
   work on one first and then the other   */

  for(long psiindex = 0; psiindex <= 1; psiindex++) {
    initializeWithBorders(psim[procid][psiindex], gp[procid].neighbors, firstrow, lastrow, firstcol, lastcol, 0.0D);
  }

  /* initialize psi matrices the same way  */

  for(long psiindex = 0; psiindex <= 1; psiindex++) {
    initializeWithBorders(psi[procid][psiindex], gp[procid].neighbors, firstrow, lastrow, firstcol, lastcol, 0.0D);
  }

  /* compute input curl of wind stress */

  double** tauz_local = (double **) tauz[procid];
  ysca1 = .5*ysca;
  double factor= -t0*pi/ysca1;

  if (left_neighbor == -1) {
    if (up_neighbor == -1) tauz_local[0][0] = 0.0;
    if (down_neighbor == -1) tauz_local[im-1][0] = 0.0;
    Matrix_set_col(double, 0.0, tauz_local, firstrow, lastrow, 0);
  }

  if (right_neighbor == -1) {
    double curlt = factor* sin(pi*((double) jm-1 + j_off)*res/ysca1);
    if (up_neighbor == -1) { tauz_local[0][jm-1] = curlt; }
    if (down_neighbor == -1) { tauz_local[im-1][jm-1] = curlt; }
    Matrix_set_col(double, curlt, tauz_local, firstrow, lastrow, jm-1);
  }

  if (up_neighbor == -1) { Vector_set(double, (factor* sin(pi*((double) _i + j_off)*res/ysca1)), tauz_local[0], firstcol, lastcol); }
  if (down_neighbor == -1) { Vector_set(double, (factor* sin(pi*((double) _i + j_off)*res/ysca1)), tauz_local[im-1], firstcol, lastcol); }
  
  Matrix_set(double, (factor * sin(pi*((double) _j + j_off)*res/ysca1)), tauz_local, firstrow, lastrow + 1, firstcol, lastcol + 1);

  BARRIER(bars->sl_onetime,nprocs);

  /***************************************************************
  one-time stuff over at this point
  ***************************************************************/

  long iday;

  bool dayflag = false;
  bool dhourflag = false;
  bool endflag = false;
  while (!endflag) {
    while ((!dayflag) || (!dhourflag)) {
      dayflag = false;
      dhourflag = false;
      if (nstep == 1) {
        /* POSSIBLE ENHANCEMENT:  Here is where one might reset the
        statistics that one is measuring about the parallel execution */
      }

      slave2(procid,firstrow,lastrow,numrows,firstcol,lastcol,numcols);

      /* update time and step number
      note that these time and step variables are private i.e. every
      process has its own copy and keeps track of its own time  */

      ttime = ttime + dtau;
      nstep = nstep + 1;
      day = ttime/86400.0;

      if (day > ((double) outday0)) {
        dayflag = 1;
        iday = (long) day;
        dhour = dhour+dtau;
        if (dhour >= 86400.0) {
	        dhourflag = true;
        }
      }
    }
    dhour = 0.0;

    addValues(psium[procid], psium[procid], psim[procid][0], gp[procid].neighbors, firstrow, lastrow, firstcol, lastcol);
    addValues(psilm[procid], psium[procid], psim[procid][1], gp[procid].neighbors, firstrow, lastrow, firstcol, lastcol);
    
    if (iday >= (long) outday3) {
      endflag = true;
    }
  }
}

