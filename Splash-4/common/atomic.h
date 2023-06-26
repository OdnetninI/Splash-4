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

#ifndef __SPLASH_4__ATOMIC_H__
#define __SPLASH_4__ATOMIC_H__

#pragma once

#define FETCH_ADD(shared, value) ({				\
      _NOTE_START_ATOMIC();					\
      uint64_t ___x = atomic_fetch_add(&(shared), value);	\
      _NOTE_END_ATOMIC();					\
      ___x;							\
    })								\

#define FETCH_SUB(shared, value) ({				\
      _NOTE_START_ATOMIC();					\
      uint64_t ___x = atomic_fetch_sub(&(shared), value);	\
      _NOTE_END_ATOMIC();					\
      ___x;							\
    })								\

#define STORE(shared, value) ({			\
      atomic_store(&(shared), value);		\
    })						\

#define LOAD(shared) ({			\
      atomic_load(&(shared));		\
    })					\
  
#define CAS(shared, oldvalue, newvalue) ({				\
      _NOTE_START_CMPXCHG();						\
      _Bool ___b = atomic_compare_exchange_weak(&(shared), &(oldvalue), newvalue); \
      _NOTE_END_CMPXCHG();						\
      ___b;								\
    })									\

// _CAS is only intended for internal use, so no tracking is being done
#define _CAS(address, oldvalue, newvalue) ({				\
      _Bool ___b = atomic_compare_exchange_weak((address), &(oldvalue), newvalue); \
      ___b;								\
    })									\

#define FETCH_ADD_DOUBLE(address, value) ({			\
      _NOTE_START_CMPXCHG();					\
      double ___oldValue = LOAD(*address);			\
      double ___newValue;					\
      do {							\
	___newValue = ___oldValue + value;			\
      } while(!_CAS(address, ___oldValue, ___newValue));	\
      _NOTE_END_CMPXCHG();					\
      ___oldValue;						\
    })								\


#define ATOMIC_MAX_DOUBLE(address, value) ({	      \
      _NOTE_START_CMPXCHG();					\
      double ___oldValue = LOAD(*address);		\
      double ___newValue = value;	           		\
      do {							      \
        if (___newValue <= ___oldValue) break;        \
      } while(!_CAS(address, ___oldValue, ___newValue));	\
      _NOTE_END_CMPXCHG();					\
      ___oldValue;						\
    })								\


#endif /* __SPLASH_4__ATOMIC_H__ */
