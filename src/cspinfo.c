/* ReSCAL - cspinfo entry
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

//numero et date de version
//#define VER_NUM "1.2"
//#define VER_DAT "2011-07"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include "defs.h"
#include "macros.h"
#include "format.h"
#include "cells.h"


int32_t H = 0, L = 0, D = 0, HL = 0, HLD = 0; // les dimensions de la terre
Cell  *TE = NULL;          // la 'terre'
double csp_time = 0.0;
char *csp_filename = NULL; //nom du fichier CSP
unsigned char opt_count = 0;
const char *etats[MAX_CELL] = ETATS;  // les noms des types de cellules
int32_t Ncel[MAX_CELL];        // nombre de cellules par type

void usage() {
  printf("CSP metadata reader");
  //printf("\nversion %s (%s)\n",VER_NUM,VER_DAT);
  printf("\nusage: cspinfo <CSP-file>[.gz] [-c]\n");
  printf("  -c \t count cells\n");
  exit(-1);
}

void general_options(int32_t argc, char *argv[]) {
  int32_t i;
  for (i = 1; i < argc; i++) {
    if (!strcmp(argv[i], "-c")) {
      opt_count = 1;
    }
  }
}

void count_cells() {
  int32_t i;

  HL = H * L;
  HLD = HL * D;
  AllocMemoryPrint("TE", TE, Cell, HLD);
  ResetMemory(TE, Cell, HLD);

  csp_set_bounds(0, 0, 0);

  read_csp(csp_filename);

  LogPrintf("counting number of cells of each type in CSP file ...\n");

  for (i = 0; i < MAX_CELL; i++) {
    Ncel[i] = 0;
  }

  for (i = 0; i < HLD; i++) {
    Ncel[TE[i].celltype]++;
  }

  for (i = 0; i < MAX_CELL; i++) {
    if ((Ncel[i] > 0) || (*etats[i] != 0)) {
      LogPrintf("type %d - %s :\t %d cells\n", i, etats[i], Ncel[i]);
    }
  }
}

int main(int argc, char **argv) {
  int32_t i;

  if (argc < 2) {
    usage();
    exit(-4);
  }

  i = 1;
  csp_filename = argv[i++];

  general_options(argc, argv);

  csp_set_warning(0);

  read_csp_header(csp_filename);

  if (opt_count) {
    count_cells();
  }

  compress(csp_filename, 0);

  return 1;
}




