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

#ifndef __SPLASH_4__PAUSE_H__
#define __SPLASH_4__PAUSE_H__

#pragma once

#define PAUSEDEC(sem) sem_t sem;
#define PAUSEINIT(sem) {sem_init(&(sem), 0, 0);}
#define CLEARPAUSE(sem) {;}

#define SETPAUSE(sem) {				\
    _NOTE_START_SEM_POST();			\
    sem_post(&(sem));				\
    _NOTE_END_SEM_POST();			\
  }						\

#define WAITPAUSE(sem) {			\
    _NOTE_START_SEM_WAIT();			\
    sem_wait(&(sem));				\
    _NOTE_END_SEM_WAIT();			\
  }						\

#endif /* __SPLASH_4__PAUSE_H__ */
