/* ReSCAL - bin2png entry
 *
 * Copyright (C) 2011
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

#define _MAIN_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include "defs.h"
//#include "bin2png.h"
#include "macros.h"
#include "cells.h"
#include "space.h"
#include "surface.h"
#include "view.h"

int32_t prog = PROG_TOOL;
unsigned char opt_nv = 1, opt_info = 0;
char *output_filename = NULL;

extern int32_t H, L, D;
extern char *bin_filename;

void usage() {
  printf("ReSCAL image generator");
  printf(", %s model", MOD_NAME);
#ifdef CELL_COLOR
  printf(", CELL_COLOR data");
#endif
#ifdef CELL_TIME
  printf(", CELL_TIME data");
#endif
  printf("\nversion %s (%s)\n", VER_NUM, VER_DAT);
  printf("size of cells : %lu byte(s)\n", sizeof(Cell));
  printf("\nusage : bin2png <binary file> H L D [GRAPHICAL OPTIONS] [-o <output PNG file>]\n");
  printf("  H : height\n");
  printf("  L : length\n");
  printf("  D : depth\n");
  show_view_options();
  exit(-1);
}


void general_options(int32_t argc, char *argv[]) {
  int32_t i;
  for (i = 1; i < argc; i++) {
    if (!strcmp(argv[i], "-o")) {
      output_filename = argv[++i];
    }
  }
}

int main(int argc, char **argv) {
  int32_t i;

  if (argc < 5) {
    usage();
    exit(-4);
  }

  i = 1;
  bin_filename = argv[i++];
  H = atoi(argv[i++]);
  L = atoi(argv[i++]);
  D = atoi(argv[i++]);

  general_options(argc, argv);

  cree_terre();

  if (!output_filename) {
    output_filename = bin_filename;
    strcpy(output_filename + strlen(output_filename) - 3, "png");
  }

  //calcule_couloir();
  compute_v_walls();
  compute_h_ceil();

#ifdef ALTI
  calcule_alti(ALTI, ALTI_MODE_BAS);
  calcule_normales();
#endif

  view_init(argc, argv);

  view_dump_init();

  LogPrintf("PNG image: %s\n", output_filename);
  //view_dump_png(output_filename);
  dump_image(output_filename, "png");

  return 1;
}


//// empty functions, for linking with simul.c

void do_thread_sched() {
}

void wait_lgca_thread() {
}

void push_status(char* str) {
  (void)str; //SUPPRESS: unused warning
}

void pop_status() {
}

void lock_display(int32_t log_flag) {
  (void)log_flag; //SUPPRESS: unused warning
}

void unlock_display(int32_t log_flag) {
  (void)log_flag; //SUPPRESS: unused warning
}

int32_t elapsed(double *sec) {
  (void)sec; //SUPPRESS: unused warning
  return 0;
}
