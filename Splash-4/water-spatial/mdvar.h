/*************************************************************************/
/*                                                                       */
/*  Copyright (c) 1994 Stanford University                               */
/*                                                                       */
/*  All rights reserved.                                                 */
/*                                                                       */
/*  Permission is given to use, copy, and modify this software for any   */
/*  non-commercial purpose as long as this copyright notice is not       */
/*  removed.  All other uses, including redistribution in whole or in    */
/*  part, are forbidden without prior written permission.                */
/*                                                                       */
/*  This software is provided with absolutely no warranty and no         */
/*  support.                                                             */
/*                                                                       */
/*************************************************************************/
#ifndef __MDVAR_H__
#define __MDVAR_H__

#pragma once
/* some variable declarations */

extern double  TEMP,RHO,TSTEP,BOXL,BOXH,CUTOFF,CUT2;
extern double  BOX_LENGTH;
extern long    NMOL,NORDER,NATMO,NATMO3,NMOL1;
extern long    BOX_PER_SIDE, BPS_SQRD;

#endif /* __MDVAR_H__ */