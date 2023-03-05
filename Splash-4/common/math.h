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

#ifndef __SPLASH_4__MATH_H__
#define __SPLASH_4__MATH_H__

#pragma once

#if (ARCH == X86_64) || (ARCH == X86)

#define LOG2_INT(type, x) ({					\
      type _x = (x);						\
      type _y;							\
      __asm__ ( "\tbsr %1, %0\n" : "=r"(_y) : "r" (_x) );	\
      _y;							\
    })								\

#define LOG2_64(x) ({ LOG2_INT(uint64_t, x); })
#define LOG2_32(x) ({ LOG2_INT(uint32_t, x); })

#else /* C version */

#define LOG2_INT(type, x) ({			\
      type _x = (x);			\
      type _y = 0;				\
      while(_x >>= 1) ++_y;			\
      _y;					\
    })						\

#define LOG2_64(x) ({ LOG2_INT(uint64_t, x); })
#define LOG2_32(x) ({ LOG2_INT(uint32_t, x); })

#endif /* ARCH */

#define LOG2(x) ({				\
      _Generic((x),				\
	       uint64_t: LOG2_64(x),		\
	       uint32_t: LOG2_32(x),		\
	       int64_t: LOG2_64(x),		\
	       int32_t: LOG2_32(x)		\
	       );				\
    })						\

#endif /* __SPLASH_4__MATH_H__ */
