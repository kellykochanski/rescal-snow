/* ReSCAL - Engine entry
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
#include <locale.h>
#include <stdint.h>

#include "defs.h"
#include "macros.h"
#include "rescal.h"
#include "cells.h"
#include "space.h"
#include "surface.h"
#include "simul.h"
#include "models.h"
#include "doublets.h"
#include "transitions.h"
#include "param.h"
#include "trace.h"
#include "lgca.h"

extern uint8_t opt_info, opt_nv;
extern uint64_t iter;        // nombre d'iterations
extern double csp_time;                    // temps reel simule
extern int32_t use_lgca;

//int32_t stop_simul = 0;                 //semaphore (pour le multithreading)
int32_t fin_simul = 0;                  //semaphore (pour l'arret des callbacks)


void rescal_params() {
  static int32_t done = 0;

  if (done) {
    return;
  }

  init_list_params();
  params_modele();
  params_simul();
#ifdef ALTI
  params_surface();
#endif
#ifdef LGCA
  params_collisions();
#endif
  params_trace();

  done = 1;
}

void rescal_usage() {
  printf("\nusage: rescal PARAMETERS_FILE [OPTIONS]\n");
  printf("       rescal -h\n");
  printf("       rescal -hm\n");
  //rescal_params();
  //param_usage();
}

//initialization of rescal thread
void rescal_init(char *param_filename) {
#ifdef NUM_MODE
  setlocale(LC_NUMERIC, NUM_MODE);
#endif

  rescal_params();
  read_parameters(param_filename);
  //init_signal();

  init_simul();

  cree_terre();
  //calcule_couloir();
  compute_v_walls();
  compute_h_ceil();
  compute_h_floor();

  init_Ncel();
  init_modele();
  init_transitions();
  init_db_pos();

  if (opt_info) {
    log_cell();
    dump_doublets();
  }

#ifdef ALTI
  calcule_alti(ALTI, ALTI_MODE_BAS);
  calcule_normales();
#endif

#ifdef LGCA
  if (use_lgca) {
    init_collisions();
  }
#endif

#if defined(MODEL_DUN) && defined(LGCA)
  //dump_cgv();
#endif

#if MODE_TRACE
  trace_init();
#endif

  LogPrintf("rescal_init: end\n");

}

int32_t rescal() {
  static uint8_t start = 0;
  uint8_t stop = 0;

  if (!start) {
    LogPrintf("debut de la simulation ...\n");
    start = 1;
  }

  while (!stop) {
#if defined(TRACE_SRC) || defined(TRACE_AIRE) || defined(TRACE_PAR)
    stop = trace_init_loop(0);
    if (stop) {
      break;
    }
#endif
    stop = simul_csp();
#if defined(TRACE_AIRE) || defined(TRACE_PAR)
    trace_end_loop();
    stop = 0;
#endif
  }

  if (stop) {
    LogPrintf("fin de la simulation\n");
    fin_simul = 1;
  }

  return 0;
}
