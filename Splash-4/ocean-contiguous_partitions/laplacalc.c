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

/* Performs the laplacian calculation for a subblock */
#include "../common/common.h"

EXTERN_ENV();

#include "decs.h"

void laplacalc(long procid, double ****x, double ****z, long psiindex, long firstrow, long lastrow, long firstcol, long lastcol) {
  double** x_local = (double **) x[procid][psiindex];
  double** z_local = (double **) z[procid][psiindex];

  long up_neighbor = gp[procid].neighbors[UP];
  long down_neighbor = gp[procid].neighbors[DOWN];
  long left_neighbor = gp[procid].neighbors[LEFT];
  long right_neighbor = gp[procid].neighbors[RIGHT];

  if (up_neighbor != -1) {
    double* x_local_0 = (double *) x_local[0];
    double* x_neighbor = (double *) x[up_neighbor][psiindex][im-2];
    for (long col = 1; col <= lastcol; col++) {
      x_local_0[col] = x_neighbor[col];
    }
  }
  if (down_neighbor != -1) {
    double* x_local_0 = (double *) x_local[im-1];
    double* x_neighbor = (double *) x[down_neighbor][psiindex][1];
    for (long col = 1; col <= lastcol; col++) {
      x_local_0[col] = x_neighbor[col];
    }
  }
  if (left_neighbor != -1) {
    double** x_neighbor = (double **) x[left_neighbor][psiindex];
    for (long row = 1; row <= lastrow; row++) {
      x_local[row][0] = x_neighbor[row][jm-2];
    }
  }
  if (right_neighbor != -1) {
    double** x_neighbor = (double **) x[right_neighbor][psiindex];
    for (long row = 1; row <= lastrow; row++) {
      x_local[row][jm-1] = x_neighbor[row][1];
    }
  }

  for (long row = firstrow; row <= lastrow; row++) {
    long row_next = row + 1;
    long row_prev = row - 1;
    double* x_local_i = (double *) x_local[row];
    double* z_local_i = (double *) z_local[row];
    double* x_local_next = (double *) x_local[row_next];
    double* x_local_prev = (double *) x_local[row_prev];
    for (long col = firstcol; col <= lastcol; col++) {
      long col_next = col + 1;
      long col_prev = col - 1;
      z_local_i[col] = factlap*(x_local_next[col]+x_local_prev[col]+x_local_i[col_next]+ x_local_i[col_prev]-4.*x_local_i[col]);
    }
  }

  if (up_neighbor == -1) {
    double* z_local_0 = (double *) z_local[0];
    for (long col = 1; col <= lastcol; col++) {
      z_local_0[col] = 0.0;
    }
  }
  if (down_neighbor == -1) {
    double* z_local_last = (double *) z_local[im-1];
    for (long col = 1; col <= lastcol; col++) {
      z_local_last[col] = 0.0;
    }
  }
  if (left_neighbor == -1) {
    for (long row = 1; row <= lastrow; row++) {
      z_local[row][0] = 0.0;
    }
  }
  if (right_neighbor == -1) {
    for (long row = 1; row <= lastrow; row++) {
      z_local[row][jm-1] = 0.0;
    }
  }

}
