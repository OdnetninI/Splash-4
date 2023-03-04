/****************************************************************************/
/*                                                                          */
/*  Copyright (c) 2023 Eduardo Jose Gomez-Hernandez (University of Murcia)  */
/*  - Rewritten from Splash-3 M4 macros
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

#ifndef __SPLASH_4__CONDVAR_H__
#define __SPLASH_4__CONDVAR_H__

#pragma once

#define CONDVARDEC(condvar) pthread_cond_t condvar;
#define CONDVARINIT(condvar) pthread_cond_init(&(condvar) , NULL);

#define CONDVARWAIT(condvar, lock) {		\
    _NOTE_START_WAIT();				\
    pthread_cond_wait(&(condvar), &(lock));	\
    _NOTE_END_WAIT();				\
  }						\

#define CONDVARSIGNAL(condvar) {		\
    _NOTE_START_SIGNAL();			\
    pthread_cond_signal(&(condvar));		\
    _NOTE_END_SIGNAL();				\
  }						\

#define CONDVARBCAST(condvar) {			\
    _NOTE_START_SIGNAL();			\
    pthread_cond_broadcast(&(condvar));		\
    _NOTE_END_SIGNAL();				\
  }						\


#endif /* __SPLASH_4__CONDVAR_H__ */
