/****************************************************************************/
/*                                                                          */
/*  Copyright (c) 2023 Eduardo Jose Gomez-Hernandez (University of Murcia)  */
/*  - Rewritten from Splash-3 M4 macros                                     */
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

#ifndef __SPLASH_4__REGION_MARKERS_H__
#define __SPLASH_4__REGION_MARKERS_H__

#pragma once

#define _NOTE_START_LOCK() {;}
#define _NOTE_END_LOCK() {;}

#define _NOTE_START_UNLOCK() {;}
#define _NOTE_END_UNLOCK() {;}

#define _NOTE_START_BARRIER() {;}
#define _NOTE_END_BARRIER() {;}

#define _NOTE_START_ATOMIC() {;}
#define _NOTE_END_ATOMIC() {;}

#define _NOTE_START_CMPXCHG() {;}
#define _NOTE_END_CMPXCHG() {;}

#define _NOTE_START_SEM_WAIT() {;}
#define _NOTE_END_SEM_WAIT() {;}

#define _NOTE_START_SEM_POST() {;}
#define _NOTE_END_SEM_POST() {;}

#define _NOTE_START_WAIT() {;}
#define _NOTE_END_WAIT() {;}

#define _NOTE_START_SIGNAL() {;}
#define _NOTE_END_SIGNAL() {;}

#define _NOTE_START_BDCAST() {;}
#define _NOTE_END_BDCAST() {;}

#endif /* __SPLASH_4__REGION_MARKERS_H__ */
