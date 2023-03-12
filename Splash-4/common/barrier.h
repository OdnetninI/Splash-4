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

#ifndef __SPLASH_4__BARRIER_H__
#define __SPLASH_4__BARRIER_H__

#pragma once

#ifdef ATOMIC_BARRIERS

#define BARDEF() 				\
    unsigned __count__;				\
    volatile int __sense__ = 1;			\
    __thread int __local_sense__ = 1;		\
  						\

#define BAREXTERN() 				\
    extern unsigned __count__;			\
    extern volatile int __sense__;		\
    extern __thread int __local_sense__;	\
  						\

#define BARRIER(barrier, numThreads) {			\
    _NOTE_START_BARRIER();				\
    __local_sense__ = !__local_sense__;			\
    if (atomic_fetch_sub(&(__count__), 1) == 1) {	\
      __count__ = numThreads;				\
      STORE(__sense__, __local_sense__);		\
    }							\
    else {						\
      do {} while (LOAD(__sense__) != __local_sense__);	\
    }							\
    _NOTE_END_BARRIER();				\
  }							\

#define BARDEC(name) ; 

#define BARINIT(barrier, numThreads) { __count__ = numThreads; }

#else /* !ATOMIC_BARRIERS */

#define BARDEF()  ; 

#define BAREXTERM()  ; 

#define BARRIER(barrier, numThreads) {					\
    _NOTE_START_BARRIER();						\
    pthread_mutex_lock(&((barrier).bar_mutex));				\
    (barrier).bar_teller++;						\
    if ((barrier).bar_teller == numThreads) {				\
      (barrier).bar_teller = 0;						\
      pthread_cond_broadcast(&((barrier).bar_cond));			\
    }									\
    else {								\
      pthread_cond_wait(&((barrier).bar_cond), &((barrier).bar_mutex));	\
    }									\
    pthread_mutex_unlock(&((barrier).bar_mutex));			\
    _NOTE_END_BARRIER();						\
  }									\

#define BARDEC(name) struct {			\
  pthread_mutex_t bar_mutex;			\
  pthread_cond_t bar_cond;			\
  unsigned bar_teller;				\
  } name;					\


#define BARINIT(barrier, numThreads) {			\
    pthread_mutex_init(&((barrier).bar_mutex), NULL);	\
    pthread_cond_init(&((barrier).bar_cond), NULL);	\
    (barrier).bar_teller = 0;				\
  }							\


#endif /* ATOMIC_BARRIERS */

#endif /* __SPLASH_4__BARRIER_H__ */
