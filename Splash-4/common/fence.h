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

#ifndef __SPLASH_4__FENCE_H__
#define __SPLASH_4__FENCE_H__

#pragma once

#define RELEASE_FENCE() { atomic_thread_fence(memory_order_release); }
#define ACQUIRE_FENCE() { atomic_thread_fence(memory_order_acquire); }
#define FULL_FENCE() { atomic_thread_fence(memory_order_seq_cst); }

#endif /* __SPLASH_4__FENCE_H__ */
