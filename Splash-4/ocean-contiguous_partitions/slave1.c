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
   desired.

   One way to do this is as follows.  The function allocate(START,SIZE,I)
   is assumed to place all addresses x such that
   (START <= x < START+SIZE) on node I.

   long d_size;
   unsigned long g_size;
   unsigned long mg_size;

   if (procid == MASTER) {
     g_size = ((jmx[numlev-1]-2)/xprocs+2)*((imx[numlev-1]-2)/yprocs+2)*sizeof(double) +
              ((imx[numlev-1]-2)/yprocs+2)*sizeof(double *);

     mg_size = numlev*sizeof(double **);
     for (i=0;i<numlev;i++) {
       mg_size+=((imx[i]-2)/yprocs+2)*((jmx[i]-2)/xprocs+2)*sizeof(double)+
                ((imx[i]-2)/yprocs+2)*sizeof(double *);
     }
     for (i= 0;i<nprocs;i++) {
       d_size = 2*sizeof(double **);
       allocate((unsigned long) psi[i],d_size,i);
       allocate((unsigned long) psim[i],d_size,i);
       allocate((unsigned long) work1[i],d_size,i);
       allocate((unsigned long) work4[i],d_size,i);
       allocate((unsigned long) work5[i],d_size,i);
       allocate((unsigned long) work7[i],d_size,i);
       allocate((unsigned long) temparray[i],d_size,i);
       allocate((unsigned long) psi[i][0],g_size,i);
       allocate((unsigned long) psi[i][1],g_size,i);
       allocate((unsigned long) psim[i][0],g_size,i);
       allocate((unsigned long) psim[i][1],g_size,i);
       allocate((unsigned long) psium[i],g_size,i);
       allocate((unsigned long) psilm[i],g_size,i);
       allocate((unsigned long) psib[i],g_size,i);
       allocate((unsigned long) ga[i],g_size,i);
       allocate((unsigned long) gb[i],g_size,i);
       allocate((unsigned long) work1[i][0],g_size,i);
       allocate((unsigned long) work1[i][1],g_size,i);
       allocate((unsigned long) work2[i],g_size,i);
       allocate((unsigned long) work3[i],g_size,i);
       allocate((unsigned long) work4[i][0],g_size,i);
       allocate((unsigned long) work4[i][1],g_size,i);
       allocate((unsigned long) work5[i][0],g_size,i);
       allocate((unsigned long) work5[i][1],g_size,i);
       allocate((unsigned long) work6[i],g_size,i);
       allocate((unsigned long) work7[i][0],g_size,i);
       allocate((unsigned long) work7[i][1],g_size,i);
       allocate((unsigned long) temparray[i][0],g_size,i);
       allocate((unsigned long) temparray[i][1],g_size,i);
       allocate((unsigned long) tauz[i],g_size,i);
       allocate((unsigned long) oldga[i],g_size,i);
       allocate((unsigned long) oldgb[i],g_size,i);
       d_size = numlev * sizeof(long);
       allocate((unsigned long) gp[i].rel_num_x,d_size,i);
       allocate((unsigned long) gp[i].rel_num_y,d_size,i);
       allocate((unsigned long) gp[i].eist,d_size,i);
       allocate((unsigned long) gp[i].ejst,d_size,i);
       allocate((unsigned long) gp[i].oist,d_size,i);
       allocate((unsigned long) gp[i].ojst,d_size,i);
       allocate((unsigned long) gp[i].rlist,d_size,i);
       allocate((unsigned long) gp[i].rljst,d_size,i);
       allocate((unsigned long) gp[i].rlien,d_size,i);
       allocate((unsigned long) gp[i].rljen,d_size,i);

       allocate((unsigned long) q_multi[i],mg_size,i);
       allocate((unsigned long) rhs_multi[i],mg_size,i);
       allocate((unsigned long) &(gp[i]),sizeof(struct Global_Private),i);
     }
   }

  */

  double** oldga_local = (double **) oldga[procid];
  double** oldgb_local = (double **) oldgb[procid];
  for (long i = 0; i < im; i++) {
    double* oldga_local_i = (double *) oldga_local[i];
    double* oldgb_local_i = (double *) oldgb_local[i];
    for (long j = 0; j < jm; j++) {
       oldga_local_i[j] = 0.0;
       oldgb_local_i[j] = 0.0;
    }
  }

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
    double* t1a = (double *) f;
    for (long iindex = 0; iindex <= jmx[numlev-1]-1; iindex++) {
      double y = ((double) iindex)*res;
      t1a[iindex] = f0+beta*(y-ysca1);
    }
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

  for(long i = istart; i <= iend; i++) {
    double* rhs_multi_local_i = (double *) rhs_multi_local[i];
    double* psib_local_i = (double *) psib_local[i];
    for(long j = jstart; j <= jend; j++) {
      rhs_multi_local_i[j] = psib_local_i[j] * ressqr;
    }
  }

  if (up_neighbor == -1) {
    double* q_multi_local_0 = (double *) q_multi_local[0];
    double* psib_local_0 = (double *) psib_local[0];
    for(long j = jstart; j <= jend; j++) {
      q_multi_local_0[j] = psib_local_0[j];
    }
  } else {
    double* psib_local_0 = (double *) psib_local[0];
    double* psib_neighbor_last = (double *) psib[up_neighbor][im-2];
    for (long i = 1; i < jm-1; i++) {
      psib_local_0[i] = psib_neighbor_last[i];
    }
  }

  if (down_neighbor == -1) {
    double* q_multi_local_last = (double *) q_multi_local[im-1];
    double* psib_local_last = (double *) psib_local[im-1];
    for(long j = jstart; j <= jend; j++) {
      q_multi_local_last[j] = psib_local_last[j];
    }
  } else {
    double* psib_local_last = (double *) psib_local[im-1];
    double* psib_neighbor_1 = (double *) psib[down_neighbor][1];
    for (long i = 1; i < jm-1; i++) {
      psib_local_last[i] = psib_neighbor_1[i];
    }
  }

  if (left_neighbor == -1) {
    for(long i = istart; i <= iend; i++) {
      q_multi_local[i][0] = psib_local[i][0];
    }
  } else {
    double** psib_neighbor = (double **) psib[left_neighbor];
    for (long i = 1; i < im-1; i++) {
      psib_local[i][0] = psib_neighbor[i][jm-2];
    }
  }

  if (right_neighbor == -1) {
    for(long i = istart; i <= iend; i++) {
      q_multi_local[i][jm-1] = psib_local[i][jm-1];
    }
  } else {
    double** psib_neighbor = (double **) psib[right_neighbor];
    for (long i = 1; i < im-1; i++) {
      psib_local[i][jm-1] = psib_neighbor[i][1];
    }
  } 

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

  for(long i = istart; i <= iend; i++) {
    double* q_multi_local_i = (double *) q_multi_local[i];
    double* psib_local_i = (double *) psib_local[i];
    for(long j = jstart; j <= jend; j++) {
      psib_local_i[j] = q_multi_local_i[j];
    }
  }

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

  if ((up_neighbor == -1) && (left_neighbor == -1)) {
    tauz_local[0][0] = 0.0;
  }
  if ((down_neighbor == -1) && (left_neighbor == -1)) {
    tauz_local[im-1][0] = 0.0;
  }
  if ((up_neighbor == -1) && (right_neighbor == -1)) {
    double sintemp = pi*((double) jm-1+j_off)*res/ysca1;
    sintemp = sin(sintemp);
    tauz_local[0][jm-1] = factor*sintemp;
  }
  if ((down_neighbor == -1) && (right_neighbor == -1)) {
    double sintemp = pi*((double) jm-1+j_off)*res/ysca1;
    sintemp = sin(sintemp);
    tauz_local[im-1][jm-1] = factor*sintemp;
  }

  if (up_neighbor == -1) {
    double* tauz_local_0 = (double *) tauz_local[0];
    for(long j=firstcol;j<=lastcol;j++) {
      double sintemp = pi*((double) j+j_off)*res/ysca1;
      sintemp = sin(sintemp);
      double curlt = factor*sintemp;
      tauz_local_0[j] = curlt;
    }
  }
  if (down_neighbor == -1) {
    double* tauz_local_last = (double *) tauz_local[im-1];
    for(long j=firstcol;j<=lastcol;j++) {
      double sintemp = pi*((double) j+j_off)*res/ysca1;
      sintemp = sin(sintemp);
      double curlt = factor*sintemp;
      tauz_local_last[j] = curlt;
    }
  }
  if (left_neighbor == -1) {
    for(long j=firstrow;j<=lastrow;j++) {
      tauz_local[j][0] = 0.0;
    }
  }
  if (right_neighbor == -1) {
    double sintemp = pi*((double) jm-1+j_off)*res/ysca1;
    sintemp = sin(sintemp);
    double curlt = factor*sintemp;
    for(long j=firstrow;j<=lastrow;j++) {
      tauz_local[j][jm-1] = curlt;
    }
  }

  for(long i = firstrow; i <= lastrow; i++) {
    double* tauz_local_i = (double *) tauz_local[i];
    for(long col = firstcol; col <= lastcol; col++) {
      double sintemp = pi*((double) col+j_off)*res/ysca1;
      sintemp = sin(sintemp);
      double curlt = factor*sintemp;
      tauz_local_i[col] = curlt;
    }
  }

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

