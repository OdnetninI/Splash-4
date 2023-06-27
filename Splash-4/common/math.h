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

#ifndef __SPLASH_4__MATH_H__
#define __SPLASH_4__MATH_H__

#pragma once

/* Typical min/max/abs ... fuctions */
#define min(a,b) ((a) < (b) ? (a) : (b))
#define max(a,b) ((a) > (b) ? (a) : (b))
#define absf(a)   ((a) > 0.0 ? (a) : -(a))
#define signf(a)   ((a) > 0.0 ? 1.0 : -1.0)

/*
  NOTE that indexes _i, _j, and _k are avaliable for the caller (example in ocean-coniguous_partitions/slace1.c) 
*/

/* Simple Vector operations */
#define Vector_set(data_type, value, ptr, i, iend) { \
  data_type* __local_ptr = (data_type*) ptr; \
  const uint64_t _start_i = i; \
  const uint64_t _end_i = iend; \
  for (uint64_t _i = _start_i; _i < _end_i; ++_i) { \
    __local_ptr[_i] = (data_type) value; \
  } \
} \

#define Vector_copy(data_type, dest, src, i, iend) { \
  data_type* __local_dest = (data_type*) dest; \
  data_type* __local_src = (data_type*) src; \
  const uint64_t _start_i = i; \
  const uint64_t _end_i = iend; \
  for (uint64_t _i = _start_i; _i < _end_i; ++_i) { \
    __local_dest[_i] = (data_type) __local_src[_i]; \
  } \
} \

/* Simple Matrix operations */
#define Matrix_set(data_type, value, ptr, i, iend, j, jend) { \
  data_type** __local_ptr = (data_type**) ptr; \
  const uint64_t _start_i = i; \
  const uint64_t _end_i = iend; \
  const uint64_t _start_j = j; \
  const uint64_t _end_j = jend; \
  for (uint64_t _i = _start_i; _i < _end_i; ++_i) { \
    data_type* __local_inner_ptr = (data_type*) __local_ptr[_i]; \
    for (uint64_t _j = _start_j; _j < _end_j; ++_j) { \
      __local_inner_ptr[_j] = (data_type) (value); \
    } \
  } \
} \

#define Matrix_multi2_set(data_type, value, ptr1, ptr2, i, iend, j, jend) { \
  data_type** __local_ptr1 = (data_type**) ptr1; \
  data_type** __local_ptr2 = (data_type**) ptr2; \
  const uint64_t _start_i = i; \
  const uint64_t _end_i = iend; \
  const uint64_t _start_j = j; \
  const uint64_t _end_j = jend; \
  for (uint64_t _i = _start_i; _i < _end_i; ++_i) { \
    data_type* __local_inner_ptr1 = (data_type*) __local_ptr1[_i]; \
    data_type* __local_inner_ptr2 = (data_type*) __local_ptr2[_i]; \
    for (uint64_t _j = _start_j; _j < _end_j; ++_j) { \
      __local_inner_ptr1[_j] = (data_type) (value); \
      __local_inner_ptr2[_j] = (data_type) (value); \
    } \
  } \
} \

#define Matrix_set_col(data_type, value, ptr, i, iend, col) { \
  data_type** __local_ptr = (data_type**) ptr; \
  const uint64_t _start_i = i; \
  const uint64_t _end_i = iend; \
  for (uint64_t _i = _start_i; _i < _end_i; ++_i) { \
    __local_ptr[_i][col] = (data_type) (value); \
  } \
} \

#define Matrix_copy_col(data_type, dest, src, i, iend, col1, col2) { \
  data_type** __local_dest = (data_type**) dest; \
  data_type** __local_src = (data_type**) src; \
  const uint64_t _start_i = i; \
  const uint64_t _end_i = iend; \
  for (uint64_t _i = _start_i; _i < _end_i; ++_i) { \
    __local_dest[_i][col1] = (data_type) __local_src[_i][col2]; \
  } \
} \

#define Matrix_copy(data_type, dest, src, i, iend, j, jend) { \
  data_type** __local_dest = (data_type**) dest; \
  data_type** __local_src = (data_type**) src; \
  const uint64_t _start_i = i; \
  const uint64_t _end_i = iend; \
  const uint64_t _start_j = j; \
  const uint64_t _end_j = jend; \
  for (uint64_t _i = _start_i; _i < _end_i; ++_i) { \
    data_type* __local_inner_dest = (data_type*) __local_dest[_i]; \
    data_type* __local_inner_src = (data_type*) __local_src[_i]; \
    for (uint64_t _j = _start_j; _j < _end_j; ++_j) { \
      __local_inner_dest[_j] = (data_type) __local_inner_src[_j]; \
    } \
  } \
} \

#define Matrix_copy_to_two_dest(data_type, dest1, dest2, src, i, iend, j, jend) { \
  data_type** __local_dest1 = (data_type**) dest1; \
  data_type** __local_dest2 = (data_type**) dest2; \
  data_type** __local_src = (data_type**) src; \
  const uint64_t _start_i = i; \
  const uint64_t _end_i = iend; \
  const uint64_t _start_j = j; \
  const uint64_t _end_j = jend; \
  for (uint64_t _i = _start_i; _i < _end_i; ++_i) { \
    data_type* __local_inner_dest1 = (data_type*) __local_dest1[_i]; \
    data_type* __local_inner_dest2 = (data_type*) __local_dest2[_i]; \
    data_type* __local_inner_src = (data_type*) __local_src[_i]; \
    for (uint64_t _j = _start_j; _j < _end_j; ++_j) { \
      __local_inner_dest1[_j] = (data_type) __local_inner_src[_j]; \
      __local_inner_dest2[_j] = (data_type) __local_inner_src[_j]; \
    } \
  } \
} \

#define Matrix_copy_coef(data_type, dest, src, coef, i, iend, j, jend) { \
  data_type** __local_dest = (data_type**) dest; \
  data_type** __local_src = (data_type**) src; \
  const uint64_t _start_i = i; \
  const uint64_t _end_i = iend; \
  const uint64_t _start_j = j; \
  const uint64_t _end_j = jend; \
  for (uint64_t _i = _start_i; _i < _end_i; ++_i) { \
    data_type* __local_inner_dest = (data_type*) __local_dest[_i]; \
    data_type* __local_inner_src = (data_type*) __local_src[_i]; \
    for (uint64_t _j = _start_j; _j < _end_j; ++_j) { \
      __local_inner_dest[_j] = (data_type) __local_inner_src[_j] * (coef); \
    } \
  } \
} \


/* LOG2 functions depending on the architecture */
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
