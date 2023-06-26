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

#ifndef __SPLASH_4__ALLOC_H__
#define __SPLASH_4__ALLOC_H__

#pragma once

#define G_MALLOC(size) ({			\
      void* mem = malloc(size);			\
      assert(mem);				\
      mem;					\
    })						\

#define NU_MALLOC(size) ({			\
      void* mem = malloc(size);			\
      assert(mem);				\
      mem;					\
    })						\
  
#define CLOCK(clock) {							     \
    struct timeval FullTime;						     \
    gettimeofday(&FullTime, NULL);					     \
    (clock) = (unsigned long)(FullTime.tv_usec + FullTime.tv_sec * 1000000); \
  }									     \
    
#endif /* __SPLASH_4__ALLOC_H__ */
