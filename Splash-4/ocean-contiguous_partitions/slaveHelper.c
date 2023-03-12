/****************************************************************************/
/*                                                                          */
/*  Copyright (c) 2023 Eduardo Jose Gomez-Hernandez (University of Murcia)  */       
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

#include "../common/common.h"

EXTERN_ENV();

#include "decs.h"

void addValues(double** A, double** B, double** C, long neighbors[8], long firstrow, long lastrow, long firstcol, long lastcol) {
  if (neighbors[UP] == -1 && neighbors[LEFT] == -1) {
    A[0][0] = B[0][0] + C[0][0];
  }
  if (neighbors[DOWN] == -1 && neighbors[LEFT] == -1) {
    A[im-1][0] = B[im-1][0] + C[im-1][0];
  }
  if (neighbors[UP] == -1 && neighbors[RIGHT] == -1) {
    A[0][jm-1] = B[0][jm-1] + C[0][jm-1];
  }
  if (neighbors[DOWN] == -1 && neighbors[RIGHT] == -1) {
    A[im-1][jm-1] = B[im-1][jm-1] + C[im-1][jm-1];
  }
  
  if (neighbors[UP] == -1) {
    double* a = (double *) A[0];
    double* b = (double *) B[0];
    double* c = (double *) C[0];
    for(long j = firstcol; j <= lastcol; j++) {
      a[j] = b[j] + c[j];
    }
  }
  if (neighbors[DOWN] == -1) {
    double* a = (double *) A[im-1];
    double* b = (double *) B[im-1];
    double* c = (double *) C[im-1];
    for(long j = firstcol; j <= lastcol; j++) {
      a[j] = b[j] + c[j];
    }
  }
  if (neighbors[LEFT] == -1) {
    for(long j = firstrow; j <= lastrow; j++) {
      A[j][0] = B[j][0] + C[j][0];
    }
  }
  if (neighbors[RIGHT] == -1) {
    for(long j = firstrow; j <= lastrow; j++) {
      A[j][jm-1] = B[j][jm-1] + C[j][jm-1];
    }
  }

  for(long i = firstrow; i <= lastrow; i++) {
    double* a = (double *) A[i];
    double* b = (double *) B[i];
    double* c = (double *) C[i];
    for(long iindex = firstcol; iindex <= lastcol; iindex++) {
      a[iindex] = b[iindex] + c[iindex];
    }
  }
}

void addValuesWithWeight(double** A, double** B, double c_weight, double** C, long neighbors[8], long firstrow, long lastrow, long firstcol, long lastcol) {
  if (neighbors[UP] == -1 && neighbors[LEFT] == -1) {
    A[0][0] = B[0][0] + c_weight * C[0][0];
  }
  if (neighbors[DOWN] == -1 && neighbors[LEFT] == -1) {
    A[im-1][0] = B[im-1][0] + c_weight * C[im-1][0];
  }
  if (neighbors[UP] == -1 && neighbors[RIGHT] == -1) {
    A[0][jm-1] = B[0][jm-1] + c_weight * C[0][jm-1];
  }
  if (neighbors[DOWN] == -1 && neighbors[RIGHT] == -1) {
    A[im-1][jm-1] = B[im-1][jm-1] + c_weight * C[im-1][jm-1];
  }
  
  if (neighbors[UP] == -1) {
    double* a = (double *) A[0];
    double* b = (double *) B[0];
    double* c = (double *) C[0];
    for(long j = firstcol; j <= lastcol; j++) {
      a[j] = b[j] + c_weight * c[j];
    }
  }
  if (neighbors[DOWN] == -1) {
    double* a = (double *) A[im-1];
    double* b = (double *) B[im-1];
    double* c = (double *) C[im-1];
    for(long j = firstcol; j <= lastcol; j++) {
      a[j] = b[j] + c_weight * c[j];
    }
  }
  if (neighbors[LEFT] == -1) {
    for(long j = firstrow; j <= lastrow; j++) {
      A[j][0] = B[j][0] + c_weight * C[j][0];
    }
  }
  if (neighbors[RIGHT] == -1) {
    for(long j = firstrow; j <= lastrow; j++) {
      A[j][jm-1] = B[j][jm-1] + c_weight * C[j][jm-1];
    }
  }

  for(long i = firstrow; i <= lastrow; i++) {
    double* a = (double *) A[i];
    double* b = (double *) B[i];
    double* c = (double *) C[i];
    for(long iindex = firstcol; iindex <= lastcol; iindex++) {
      a[iindex] = b[iindex] + c_weight * c[iindex];
    }
  }
}

void addValuesWith2Weight(double** A, double b_weight, double** B, double c_weight, double** C, long neighbors[8], long firstrow, long lastrow, long firstcol, long lastcol) {
  if (neighbors[UP] == -1 && neighbors[LEFT] == -1) {
    A[0][0] = b_weight * B[0][0] + c_weight * C[0][0];
  }
  if (neighbors[DOWN] == -1 && neighbors[LEFT] == -1) {
    A[im-1][0] = b_weight * B[im-1][0] + c_weight * C[im-1][0];
  }
  if (neighbors[UP] == -1 && neighbors[RIGHT] == -1) {
    A[0][jm-1] = b_weight * B[0][jm-1] + c_weight * C[0][jm-1];
  }
  if (neighbors[DOWN] == -1 && neighbors[RIGHT] == -1) {
    A[im-1][jm-1] = b_weight * B[im-1][jm-1] + c_weight * C[im-1][jm-1];
  }
  
  if (neighbors[UP] == -1) {
    double* a = (double *) A[0];
    double* b = (double *) B[0];
    double* c = (double *) C[0];
    for(long j = firstcol; j <= lastcol; j++) {
      a[j] = b_weight * b[j] + c_weight * c[j];
    }
  }
  if (neighbors[DOWN] == -1) {
    double* a = (double *) A[im-1];
    double* b = (double *) B[im-1];
    double* c = (double *) C[im-1];
    for(long j = firstcol; j <= lastcol; j++) {
      a[j] = b_weight * b[j] + c_weight * c[j];
    }
  }
  if (neighbors[LEFT] == -1) {
    for(long j = firstrow; j <= lastrow; j++) {
      A[j][0] = b_weight * B[j][0] + c_weight * C[j][0];
    }
  }
  if (neighbors[RIGHT] == -1) {
    for(long j = firstrow; j <= lastrow; j++) {
      A[j][jm-1] = b_weight * B[j][jm-1] + c_weight * C[j][jm-1];
    }
  }

  for(long i = firstrow; i <= lastrow; i++) {
    double* a = (double *) A[i];
    double* b = (double *) B[i];
    double* c = (double *) C[i];
    for(long iindex = firstcol; iindex <= lastcol; iindex++) {
      a[iindex] = b_weight * b[iindex] + c_weight * c[iindex];
    }
  }
}

void subValues(double** A, double** B, double** C, long neighbors[8], long firstrow, long lastrow, long firstcol, long lastcol) {
  if (neighbors[UP] == -1 && neighbors[LEFT] == -1) {
    A[0][0] = B[0][0] - C[0][0];
  }
  if (neighbors[DOWN] == -1 && neighbors[LEFT] == -1) {
    A[im-1][0] = B[im-1][0] - C[im-1][0];
  }
  if (neighbors[UP] == -1 && neighbors[RIGHT] == -1) {
    A[0][jm-1] = B[0][jm-1] - C[0][jm-1];
  }
  if (neighbors[DOWN] == -1 && neighbors[RIGHT] == -1) {
    A[im-1][jm-1] = B[im-1][jm-1] - C[im-1][jm-1];
  }
  
  if (neighbors[UP] == -1) {
    double* a = (double *) A[0];
    double* b = (double *) B[0];
    double* c = (double *) C[0];
    for(long j = firstcol; j <= lastcol; j++) {
      a[j] = b[j] - c[j];
    }
  }
  if (neighbors[DOWN] == -1) {
    double* a = (double *) A[im-1];
    double* b = (double *) B[im-1];
    double* c = (double *) C[im-1];
    for(long j = firstcol; j <= lastcol; j++) {
      a[j] = b[j] - c[j];
    }
  }
  if (neighbors[LEFT] == -1) {
    for(long j = firstrow; j <= lastrow; j++) {
      A[j][0] = B[j][0] - C[j][0];
    }
  }
  if (neighbors[RIGHT] == -1) {
    for(long j = firstrow; j <= lastrow; j++) {
      A[j][jm-1] = B[j][jm-1] - C[j][jm-1];
    }
  }

  for(long i = firstrow; i <= lastrow; i++) {
    double* a = (double *) A[i];
    double* b = (double *) B[i];
    double* c = (double *) C[i];
    for(long iindex = firstcol; iindex <= lastcol; iindex++) {
      a[iindex] = b[iindex] - c[iindex];
    }
  }
}

void subValuesWithWeight(double** A, double** B, double c_weight, double** C, long neighbors[8], long firstrow, long lastrow, long firstcol, long lastcol) {
  if (neighbors[UP] == -1 && neighbors[LEFT] == -1) {
    A[0][0] = B[0][0] - c_weight * C[0][0];
  }
  if (neighbors[DOWN] == -1 && neighbors[LEFT] == -1) {
    A[im-1][0] = B[im-1][0] - c_weight * C[im-1][0];
  }
  if (neighbors[UP] == -1 && neighbors[RIGHT] == -1) {
    A[0][jm-1] = B[0][jm-1] - c_weight * C[0][jm-1];
  }
  if (neighbors[DOWN] == -1 && neighbors[RIGHT] == -1) {
    A[im-1][jm-1] = B[im-1][jm-1] - c_weight * C[im-1][jm-1];
  }
  
  if (neighbors[UP] == -1) {
    double* a = (double *) A[0];
    double* b = (double *) B[0];
    double* c = (double *) C[0];
    for(long j = firstcol; j <= lastcol; j++) {
      a[j] = b[j] - c_weight * c[j];
    }
  }
  if (neighbors[DOWN] == -1) {
    double* a = (double *) A[im-1];
    double* b = (double *) B[im-1];
    double* c = (double *) C[im-1];
    for(long j = firstcol; j <= lastcol; j++) {
      a[j] = b[j] - c_weight * c[j];
    }
  }
  if (neighbors[LEFT] == -1) {
    for(long j = firstrow; j <= lastrow; j++) {
      A[j][0] = B[j][0] - c_weight * C[j][0];
    }
  }
  if (neighbors[RIGHT] == -1) {
    for(long j = firstrow; j <= lastrow; j++) {
      A[j][jm-1] = B[j][jm-1] - c_weight * C[j][jm-1];
    }
  }

  for(long i = firstrow; i <= lastrow; i++) {
    double* a = (double *) A[i];
    double* b = (double *) B[i];
    double* c = (double *) C[i];
    for(long iindex = firstcol; iindex <= lastcol; iindex++) {
      a[iindex] = b[iindex] - c_weight * c[iindex];
    }
  }
}

void copyValues(double** A, double** B, long neighbors[8], long firstrow, long lastrow, long firstcol, long lastcol) {
  if (neighbors[UP] == -1 && neighbors[LEFT] == -1) {
    A[0][0] = B[0][0];
  }
  if (neighbors[DOWN] == -1 && neighbors[LEFT] == -1) {
    A[im-1][0] = B[im-1][0];
  }
  if (neighbors[UP] == -1 && neighbors[RIGHT] == -1) {
    A[0][jm-1] = B[0][jm-1];
  }
  if (neighbors[DOWN] == -1 && neighbors[RIGHT] == -1) {
    A[im-1][jm-1] = B[im-1][jm-1];
  }
  
  if (neighbors[UP] == -1) {
    double* a = (double *) A[0];
    double* b = (double *) B[0];
    for(long j = firstcol; j <= lastcol; j++) {
      a[j] = b[j];
    }
  }
  if (neighbors[DOWN] == -1) {
    double* a = (double *) A[im-1];
    double* b = (double *) B[im-1];
    for(long j = firstcol; j <= lastcol; j++) {
      a[j] = b[j];
    }
  }
  if (neighbors[LEFT] == -1) {
    for(long j = firstrow; j <= lastrow; j++) {
      A[j][0] = B[j][0];
    }
  }
  if (neighbors[RIGHT] == -1) {
    for(long j = firstrow; j <= lastrow; j++) {
      A[j][jm-1] = B[j][jm-1];
    }
  }

  for(long i = firstrow; i <= lastrow; i++) {
    double* a = (double *) A[i];
    double* b = (double *) B[i];
    for(long iindex = firstcol; iindex <= lastcol; iindex++) {
      a[iindex] = b[iindex];
    }
  }
}

double calculatePsiipriv(double** A, long neighbors[8], long firstrow, long lastrow, long firstcol, long lastcol) {
  double psiaipriv = 0.0;
  if ((neighbors[UP] == -1) && (neighbors[LEFT] == -1)) {
    psiaipriv = psiaipriv + 0.25*(A[0][0]);
  }
  if ((neighbors[UP] == -1) && (neighbors[RIGHT] == -1)) {
    psiaipriv = psiaipriv + 0.25*(A[0][jm-1]);
  }
  if ((neighbors[DOWN] == -1) && (neighbors[LEFT] == -1)) {
    psiaipriv=psiaipriv+0.25*(A[im-1][0]);
  }
  if ((neighbors[DOWN] == -1) && (neighbors[RIGHT] == -1)) {
    psiaipriv=psiaipriv+0.25*(A[im-1][jm-1]);
  }
  if (neighbors[UP] == -1) {
    double* a = (double *) A[0];
    for(long j = firstcol; j <= lastcol; j++) {
      psiaipriv = psiaipriv + 0.5*a[j];
    }
  }
  if (neighbors[DOWN] == -1) {
    double* a = (double *) A[im-1];
    for(long j = firstcol; j <= lastcol; j++) {
      psiaipriv = psiaipriv + 0.5*a[j];
    }
  }
  if (neighbors[LEFT] == -1) {
    for(long j = firstrow; j <= lastrow; j++) {
      psiaipriv = psiaipriv + 0.5*A[j][0];
    }
  }
  if (neighbors[RIGHT] == -1) {
    for(long j = firstrow; j <= lastrow; j++) {
      psiaipriv = psiaipriv + 0.5*A[j][jm-1];
    }
  }
  for(long i = firstrow; i <= lastrow; i++) {
    double* a = (double *) A[i];
    for(long iindex = firstcol; iindex <= lastcol; iindex++) {
      psiaipriv = psiaipriv + a[iindex];
    }
  }
  return psiaipriv;
}

void initializeWithBorders(double** A, long neighbors[8], long firstrow, long lastrow, long firstcol, long lastcol, double initValue) {
  if (neighbors[UP] == -1 && neighbors[LEFT] == -1) {
    A[0][0] = initValue;
  }
  if (neighbors[DOWN] == -1 && neighbors[LEFT] == -1) {
    A[im-1][0] = initValue;
  }
  if (neighbors[UP] == -1 && neighbors[RIGHT] == -1) {
    A[0][jm-1] = initValue;
  }
  if (neighbors[DOWN] == -1 && neighbors[RIGHT] == -1) {
    A[im-1][jm-1] = initValue;
  }

  if (neighbors[UP] == -1) {
    double* a = (double *) A[0];
    for(long j = firstcol; j <= lastcol; j++) {
      a[j] = initValue;
    }
  }
  if (neighbors[DOWN] == -1) {
    double* a = (double *) A[im-1];
    for(long j = firstcol; j <= lastcol; j++) {
      a[j] = initValue;
    }
  }
  if (neighbors[LEFT] == -1) {
    for(long j = firstrow; j <= lastrow; j++) {
      A[j][0] = initValue;
    }
  }
  if (neighbors[RIGHT] == -1) {
    for(long j = firstrow; j <= lastrow; j++) {
      A[j][jm-1] = initValue;
    }
  }

  for(long i = firstrow; i <= lastrow; i++) {
    double* a = (double *) A[i];
    for(long iindex = firstcol; iindex <= lastcol; iindex++) {
      a[iindex] = 0.0;
    }
  }
}
