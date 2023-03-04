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

#ifndef __SPLASH_4__THREAD_H__
#define __SPLASH_4__THREAD_H__

#pragma once

#ifdef BIND_CORES

#define BIND(cpu, thread) {						\
    cpu_set_t cpuset;							\
    const pthread_t pid = thread;					\
    cpu_set_t  ____cpuset;						\
    CPU_ZERO(&____cpuset);						\
    CPU_SET((cpu), &____cpuset);					\
    const int set_result = pthread_setaffinity_np(pid, sizeof(cpu_set_t), &____cpuset); \
    assert(set_result == 0);						\
  }									\

#elif defined(BIND_THREADS)

#define BIND(cpu, thread) {						\
    cpu_set_t cpuset;							\
    const pthread_t pid = thread;					\
    cpu_set_t  ____cpuset;						\
    CPU_ZERO(&____cpuset);						\
    CPU_SET((cpu)/2+((cpu)%2), &____cpuset);				\
    const int set_result = pthread_setaffinity_np(pid, sizeof(cpu_set_t), &____cpuset); \
    assert(set_result == 0);						\
  }									\

#else /* NO_BIND */

#define BIND(cpu, thread) {;}

#endif /* BIND_CORES BIND_THREADS */

#define CREATE(function, numThreads) {					\
    int i;								\
    assert(__threads__<__MAX_THREADS__);				\
    pthread_mutex_lock(&__intern__);					\
    for (i = 0; i < (numThreads) - 1; i++) {				\
      BIND(i, __tid__[__threads__-1])					\
	const int Error = pthread_create(&__tid__[__threads__++], NULL, (void * (*)(void *))(function), NULL); \
      if (Error != 0) {							\
	printf("Error in pthread_create().\n");				\
	exit(-1);							\
      }									\
    }									\
    pthread_mutex_unlock(&__intern__);					\
    BIND(i, __tid__[__threads__-1]);					\
    function();								\
  }									\

#define WAIT_FOR_END(numThreads) {		\
    int aantal = numThreads;			\
    while (aantal--)				\
      pthread_join(__tid__[aantal], NULL);	\
  }						\

#endif /* __SPLASH_4__THREAD_H__ */
