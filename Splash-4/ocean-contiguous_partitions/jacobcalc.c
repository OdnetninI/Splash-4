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

/* Does the arakawa jacobian calculation (of the x and y matrices,
  putting the results in the z matrix) for a subblock.  */
#include "../common/common.h"

EXTERN_ENV();

#include "decs.h"

void jacobcalc(double ***x, double ***y, double ***z, long pid, long firstrow, long lastrow, long firstcol, long lastcol) {
  double** x_local = (double **) x[pid];
  double** y_local = (double **) y[pid];
  double** z_local = (double **) z[pid];

  long up_left_neighbor = gp[pid].neighbors[UPLEFT];
  long up_right_neighbor = gp[pid].neighbors[UPRIGHT];
  long down_left_neighbor = gp[pid].neighbors[DOWNLEFT];
  long down_right_neighbor = gp[pid].neighbors[DOWNRIGHT];
  long up_neighbor = gp[pid].neighbors[UP];
  long left_neighbor = gp[pid].neighbors[LEFT];
  long right_neighbor = gp[pid].neighbors[RIGHT];
  long down_neighbor = gp[pid].neighbors[DOWN];

  if ((up_neighbor == -1) && (left_neighbor == -1)) {
    z_local[0][0] = 0.0;
  }
  if ((down_neighbor == -1) && (left_neighbor == -1)) {
    z_local[im-1][0] = 0.0;
  }
  if ((up_neighbor == -1) && (right_neighbor == -1)) {
    z_local[0][jm-1] = 0.0;
  }
  if ((down_neighbor == -1) && (right_neighbor == -1)) {
    z_local[im-1][jm-1] = 0.0;
  }

  if (up_left_neighbor != -1) {
    x_local[0][0] = x[up_left_neighbor][im-2][jm-2];
    y_local[0][0] = y[up_left_neighbor][im-2][jm-2];
  }
  if (up_right_neighbor != -1) {
    x_local[0][jm-1] = x[up_right_neighbor][im-2][1];
    y_local[0][jm-1] = y[up_right_neighbor][im-2][1];
  }
  if (down_left_neighbor != -1) {
    x_local[im-1][0] = x[down_left_neighbor][1][jm-2];
    y_local[im-1][0] = y[down_left_neighbor][1][jm-2];
  }
  if (down_right_neighbor != -1) {
    x_local[im-1][jm-1] = x[down_right_neighbor][1][1];
    y_local[im-1][jm-1] = y[down_right_neighbor][1][1];
  }

  if (up_neighbor == -1) {
    if (left_neighbor != -1) {
      x_local[0][0] = x[left_neighbor][0][jm-2];
      y_local[0][0] = y[left_neighbor][0][jm-2];
    } else if (down_neighbor != -1) {
      x_local[im-1][0] = x[down_neighbor][1][0];
      y_local[im-1][0] = y[down_neighbor][1][0];
    }
    if (right_neighbor != -1) {
      x_local[0][jm-1] = x[right_neighbor][0][1];
      y_local[0][jm-1] = y[right_neighbor][0][1];
    } else if (down_neighbor != -1) {
      x_local[im-1][jm-1] = x[down_neighbor][1][jm-1];
      y_local[im-1][jm-1] = y[down_neighbor][1][jm-1];
    }
  } else if (down_neighbor == -1) {
    if (left_neighbor != -1) {
      x_local[im-1][0] = x[left_neighbor][im-1][jm-2];
      y_local[im-1][0] = y[left_neighbor][im-1][jm-2];
    } else if (up_neighbor != -1) {
      x_local[0][0] = x[up_neighbor][im-2][0];
      y_local[0][0] = y[up_neighbor][im-2][0];
    }
    if (right_neighbor != -1) {
      x_local[im-1][jm-1] = x[right_neighbor][im-1][1];
      y_local[im-1][jm-1] = y[right_neighbor][im-1][1];
    } else if (up_neighbor != -1) {
      x_local[0][jm-1] = x[up_neighbor][im-2][jm-1];
      y_local[0][jm-1] = y[up_neighbor][im-2][jm-1];
    }
  } else if (left_neighbor == -1) {
    if (left_neighbor != -1) {
      x_local[0][0] = x[left_neighbor][im-2][0];
      y_local[0][0] = y[left_neighbor][im-2][0];
    }
    if (down_neighbor != -1) {
      x_local[im-1][0] = x[down_neighbor][1][0];
      y_local[im-1][0] = y[down_neighbor][1][0];
    }
  } else if (right_neighbor == -1) {
    if (up_neighbor != -1) {
      x_local[0][jm-1] = x[up_neighbor][im-2][jm-1];
      y_local[0][jm-1] = y[up_neighbor][im-2][jm-1];
    }
    if (down_neighbor != -1) {
      x_local[im-1][jm-1] = x[down_neighbor][1][jm-1];
      y_local[im-1][jm-1] = y[down_neighbor][1][jm-1];
    }
  }

  if (up_neighbor != -1) {
    double* x_local_0 = (double *) x_local[0];
    double* x_neighbor = (double *) x[up_neighbor][im-2];
    for (long col = 1; col <= lastcol; col++) {
      x_local_0[col] = x_neighbor[col];
    }
    double* y_local_0 = (double *) y_local[0];
    double* y_neighbor = (double *) y[up_neighbor][im-2];
    for (long col = 1; col <= lastcol; col++) {
      y_local_0[col] = y_neighbor[col];
    }
  }
  if (down_neighbor != -1) {
    double* x_local_last = (double *) x_local[im-1];
    double* x_neighbor = (double *) x[down_neighbor][1];
    for (long col = 1; col <= lastcol; col++) {
      x_local_last[col] = x_neighbor[col];
    }
    double* y_local_last = (double *) y_local[im-1];
    double* y_neighbor = (double *) y[down_neighbor][1];
    for (long col = 1; col <= lastcol; col++) {
      y_local_last[col] = y_neighbor[col];
    }
  }
  if (left_neighbor != -1) {
    double** x_neighbor = (double **) x[left_neighbor];
    for (long row = 1; row <= lastrow; row++) {
      x_local[row][0] = x_neighbor[row][jm-2];
    }
    double** y_neighbor = (double **) y[left_neighbor];
    for (long row = 1; row <= lastrow; row++) {
      y_local[row][0] = y_neighbor[row][jm-2];
    }
  }
  if (right_neighbor != -1) {
    double** x_neighbor = (double **) x[right_neighbor];
    for (long row = 1; row <= lastrow; row++) {
      x_local[row][jm-1] = x_neighbor[row][1];
    }
    double** y_neighbor = (double **) y[right_neighbor];
    for (long row = 1; row <= lastrow; row++) {
      y_local[row][jm-1] = y_neighbor[row][1];
    }
  }

  for (long row = firstrow; row <= lastrow; row++) {
    long next_row = row+1;
    long prev_row = row-1;
    double* x_local_i = (double *) x_local[row];
    double* y_local_i = (double *) y_local[row];
    double* z_local_i = (double *) z_local[row];
    double* y_local_next = (double *) y_local[next_row];
    double* y_local_prev = (double *) y_local[prev_row];
    double* x_local_next = (double *) x_local[next_row];
    double* x_local_prev = (double *) x_local[prev_row];
    for (long col = firstcol; col <= lastcol; col++) {
      long col_next = col+1;
      long col_prev = col-1;
      double f1 = (y_local_i[col_prev]    + y_local_next[col_prev] - y_local_i[col_next]    - y_local_next[col_next]) * (x_local_next[col]      - x_local_i[col]        );
      double f2 = (y_local_prev[col_prev] + y_local_i[col_prev]    - y_local_prev[col_next] - y_local_i[col_next]   ) * (x_local_i[col]         - x_local_prev[col]     );
      double f3 = (y_local_next[col]      + y_local_next[col_next] - y_local_prev[col]      - y_local_prev[col_next]) * (x_local_i[col_next]    - x_local_i[col]        );
      double f4 = (y_local_next[col_prev] + y_local_next[col]      - y_local_prev[col_prev] - y_local_prev[col]     ) * (x_local_i[col]         - x_local_i[col_prev]   );
      double f5 = (y_local_next[col]      - y_local_i[col_next]                                                     ) * (x_local_next[col_next] - x_local_i[col]        );
      double f6 = (y_local_i[col_prev]    - y_local_prev[col]                                                       ) * (x_local_i[col]         - x_local_prev[col_prev]);
      double f7 = (y_local_i[col_next]    - y_local_prev[col]                                                       ) * (x_local_prev[col_next] - x_local_i[col]        );
      double f8 = (y_local_next[col]      - y_local_i[col_prev]                                                     ) * (x_local_i[col]         - x_local_next[col_prev]);

      z_local_i[col] = factjacob*(f1+f2+f3+f4+f5+f6+f7+f8);
    }
  }

  if (up_neighbor == -1) {
    double* z_local_0 = (double *) z_local[0];
    for (long col = firstcol; col <= lastcol; col++) {
      z_local_0[col] = 0.0;
    }
  }
  if (down_neighbor == -1) {
    double* z_local_last = (double *) z_local[im-1];
    for (long col = firstcol; col <= lastcol; col++) {
      z_local_last[col] = 0.0;
    }
  }
  if (left_neighbor == -1) {
    for (long row = firstrow; row <= lastrow; row++) {
      z_local[row][0] = 0.0;
    }
  }
  if (right_neighbor == -1) {
    for (long row = firstrow; row <= lastrow; row++) {
      z_local[row][jm-1] = 0.0;
    }
  }
}

