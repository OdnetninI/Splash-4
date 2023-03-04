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

#ifndef __SPLASH_4__COMMON_H__
#define __SPLASH_4__COMMON_H__

#pragma once

#include <stdlib.h>
#include <semaphore.h>
#include <assert.h>
#if __STDC_VERSION__ >= 201112L
#include <stdatomic.h>
#endif
#include <stdint.h>

#include <pthread.h>
#include <sched.h>
#include <unistd.h>

#define PAGE_SIZE 4096
#define __MAX_THREADS__ 256


#include "roi.h"
#include "region_markers.h"

#include "alloc.h"
#include "atomic.h"
#include "barrier.h"
#include "condvar.h"
#include "fence.h"
#include "lock.h"
#include "pause.h"
#include "thread.h"

#define MAIN_INITENV(u,u2) {			\
    __tid__[__threads__++] = pthread_self();	\
  }						\

#define MAIN_END() {				\
    exit(0);					\
  }						\

#define MAIN_ENV() 				\
  pthread_t __tid__[__MAX_THREADS__];		\
  unsigned __threads__ = 0;			\
  pthread_mutex_t __intern__;			\
  BARDEF();					\
  						\

#define EXTERN_ENV() 				\
  extern pthread_t __tid__[__MAX_THREADS__];	\
  extern unsigned __threads__;			\
  extern pthread_mutex_t __intern__;		\
  BAREXTERN();					\
						\
    
#endif /* __SPLASH_4__COMMON_H__ */
