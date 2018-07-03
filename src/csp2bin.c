/* ReSCAL - csp2bin entry
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
#include "macros.h"
#include "format.h"
#include "cells.h"


int32_t H = 0, L = 0, D = 0, HL = 0, HLD = 0; // les dimensions de la terre
Cell  *TE = NULL;          // la 'terre'
double csp_time = 0.0;
char *output_filename = NULL;

void usage() {
  printf("CSP to BIN conversion tool");
#ifdef CELL_COLOR
  printf(", CELL_COLOR data");
#endif
#ifdef CELL_TIME
  printf(", CELL_TIME data");
#endif  //printf("\nversion %s (%s)\n",VER_NUM,VER_DAT);
  printf("\nusage: csp2bin <CSP-file>[.gz] [-o <output BIN file>]\n");
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
  char *csp_filename;

  if (argc < 2) {
    usage();
    exit(-4);
  }

  i = 1;
  csp_filename = argv[i++];

  general_options(argc, argv);

  read_csp_header(csp_filename);

  HL = H * L;
  HLD = HL * D;

  AllocMemoryPrint("TE", TE, Cell, HLD);

  csp_set_bounds(0, 0, 0);

  read_csp(csp_filename);

  //compress(csp_filename);

  if (!output_filename) {
    output_filename = csp_filename;
    strcpy(output_filename + strlen(output_filename) - 3, "bin");
  }

  write_csp(DUMP_BIN, output_filename);

  return 1;
}




