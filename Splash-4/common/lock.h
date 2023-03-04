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

#ifndef __SPLASH_4__LOCK_H__
#define __SPLASH_4__LOCK_H__

#pragma once

#define LOCKDEC(lock) pthread_mutex_t lock;
#define LOCKINIT(lock) { pthread_mutex_init(&(lock), NULL); }

#define LOCK(lock) {				\
    _NOTE_START_LOCK();				\
    pthread_mutex_lock(&(lock));		\
    _NOTE_END_LOCK();				\
  }						\

#define UNLOCK(lock) {				\
    _NOTE_START_UNLOCK();			\
    pthread_mutex_unlock(&(lock));		\
    _NOTE_END_UNLOCK();				\
  }						\

#define ALOCKDEC(lock, index) pthread_mutex_t (lock)[index];
#define ALOCKINIT(lock, numLocks) {		\
    for (int i = 0; i < (numLocks); ++i)	\
      pthread_mutex_init(&((lock)[i]), NULL);	\
  }						\

#define ALOCK(lock, index) {			\
    _NOTE_START_LOCK();				\
    pthread_mutex_lock(&((lock)[(index)]));	\
    _NOTE_END_LOCK();				\
  }						\

#define AUNLOCK(lock, index) {			\
    _NOTE_START_UNLOCK();			\
    pthread_mutex_unlock(&((lock)[(index)]));	\
    _NOTE_END_UNLOCK();				\
  }						\

#define AGETL(lock, index) ((lock)[index])

#endif /* __SPLASH_4__LOCK_H__ */
