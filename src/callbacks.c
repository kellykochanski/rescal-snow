/* ReSCAL - GTK callbacks
 *
 * Copyright (C) 2011-2012
 *
 * Author: Olivier Rozier <rozier@ipgp.fr>
 *
 * This file is part of ReSCAL.
 *
 * ReSCAL is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * aint64_t with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */



#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <pthread.h>
#include <sys/time.h>
#include <stdint.h>

#include "defs.h"
#include "macros.h"
#include "entry.h"
#include "callbacks.h"
#include "rescal.h"
#include "format.h"
#include "view.h"
#include "cells.h"
#include "doublets.h"
#include "transitions.h"
#include "space.h"
#include "surface.h"
#include "simul.h"
#include "trace.h"
#ifdef LGCA
#include "lgca.h"
#endif

extern int32_t H, L, D;
extern float frame_delay;
extern uint64_t iter, md_iter;
extern uint8_t opt_quit;
extern double csp_time;
extern int32_t abs_cv, prof_cv;
extern uint8_t opt_cv;
extern int32_t vdir_mode; //display mode of the current orientation
extern uint8_t reorient_flag;
extern uint8_t csphpp_flag;
#ifdef LGCA
extern int32_t use_lgca;
extern int32_t col_iter;
extern float meanvel, maxvel;
#endif
#ifdef ROTATIONS
extern int32_t rot_mode;
#endif

int32_t end_of_rescal = 0;
int32_t rescal_paused = 0;
volatile int32_t dump_wait = 0;

pthread_mutex_t mutex_pause;
pthread_cond_t cond_pause;
pthread_mutex_t mutex_display;
#ifdef PARALLEL_AUTOMATA
pthread_barrier_t lgca_barrier;
pthread_t thread_id_lgca;
#endif // LGCA

int32_t elapsed(double *sec) {
  struct timeval t;
  struct timezone tz;

  int32_t stat;
  stat = gettimeofday(&t, &tz);
  *sec = (double)(t.tv_sec + t.tv_usec / 1000000.0);
  return (stat);
}

void *rescal_thread(void *data) {
  (void)data; //SUPPRESS: unused warning
  sleep(1);
  rescal();
  end_of_rescal = 1;
  pthread_exit(0);
  return 0;
}

#ifdef PARALLEL_AUTOMATA
void *lgca_thread(void *data) {
  while (1) {
    simul_lgca();
    pthread_barrier_wait(&lgca_barrier);
  }
  pthread_exit(0);
  return 0;
}

void wait_lgca_thread() {
  static int32_t nRetVal = -1;

  /// wait for the end of lgca cycle
  if (nRetVal == 0) {
    pthread_barrier_wait(&lgca_barrier);
  }

  /// start lgca thread
  if (nRetVal == -1) {
    pthread_barrier_init(&lgca_barrier, NULL, 2);
    nRetVal = pthread_create(&thread_id_lgca, 0, lgca_thread, 0);
    if (nRetVal != 0) {
      ErrPrintf("ERROR: cannot create lgca_thread\n");
      exit(-1);
    }
  }
}

#endif

void do_thread_sched() {
  int32_t unlocked = 0;
  if (dump_wait) {
    unlocked = 1;
    unlock_csp(0);
    while (dump_wait) {
      sched_yield();
    }
  }
  if (rescal_paused) {
    LogPrintf("paused\n");
    if (!unlocked) {
      unlocked = 1;
      unlock_csp(0);
    }
    pthread_cond_wait(&cond_pause, &mutex_pause);
    LogPrintf("resume\n");
  }
  if (unlocked) {
    lock_csp(0);
    unlocked = 0;
  }
}

void log_info() {
  dump_time();
#ifdef INFO_CEL
  log_cell();
#endif
#ifdef INFO_DBL
  dump_db_info();
#endif
#ifdef INFO_TRANS
  dump_trans_info();
#endif
#if defined(TRACE_TRANS) || defined(TRACE3D_CEL) || defined(TRACE_FLUX)
  trace_dump(1);
#endif
#ifdef LGCA
  if (use_lgca) {
    dump_densite();
#ifndef STABILITY_ANALYSIS
    dump_vel();
#endif
#ifdef CGV
    dump_cgv_coef();
#endif
  }
#endif // LGCA
}

void* do_png(void* delay) {
  while (1) {
    dump_image_inter((intptr_t) delay, "png");
    sleep((intptr_t) delay);
  }
  return 0;
}

void* do_jpeg(void* delay) {
  while (1) {
    dump_image_inter((intptr_t) delay, "jpeg");
    sleep((intptr_t) delay);
  }
  return 0;
}

void set_ss_timeout(int32_t delay, const char* type) {
  pthread_t pth;
  if (strcmp(type, "png") == 0) {
    pthread_create(&pth, 0, do_png, (void*)(intptr_t) delay);
  } else if (strcmp(type, "jpeg") == 0) {
    pthread_create(&pth, 0, do_jpeg, (void*)(intptr_t) delay);
  }
}

void* do_stop(void* arg) {
  sleep((intptr_t) arg);
  LogPrintf("time to quit !\n");
  //push_status("stopping ...", 0);
  exit(-1);
  return 0;
}

void* do_quit(void* arg) {
  (void)arg; //SUPPRESS: unused warning
  sleep(1);
  if (opt_quit && end_of_rescal) {
    LogPrintf("quit\n");
    exit(-1);
  }
  return 0;
}
