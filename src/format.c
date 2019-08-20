/* ReSCAL - CSP Format
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

#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <assert.h>
#include <stdint.h>

#include "defs.h"
#include "macros.h"
#include "format.h"
#include "cells.h"

#define READ_INT(_CH) (int)(*(int*)(_CH))
#define READ_DOUBLE(_CH) (double)(*(double*)(_CH))

extern int32_t H, L, D, HL, HLD;       // les dimensions de la terre
extern Cell  *TE;            // la 'terre'
extern double csp_time;

char header[2048];  //metadata header
int32_t csp_cell_size = 0; //size of cells in a CSP file
int32_t csp_header_size = 0; //size of CSP header
int32_t csp_L = 0;
int32_t csp_H = 0;
int32_t csp_D = 0;
int32_t csp_x_bounds = 0; //east_west thickness of the boundary walls
int32_t csp_y_bounds = 0; //north-south thickness of the boundary walls
int32_t csp_z_bounds = 0; //vertical thickness of the boundary walls
uint8_t decomp_flag = 0;
uint8_t check_model = 1;

uint8_t little_endian() {
  int32_t n = 1;

  // return 1 if litlle endian, 0 otherwise
  return *(char*)&n;
}

void csp_set_warning(uint8_t check) {
  check_model = check;
}

void csp_set_bounds(int32_t xb, int32_t yb, int32_t zb) {
  csp_x_bounds = xb;
  csp_y_bounds = yb;
  csp_z_bounds = zb;
  csp_L = L - 2 * xb;
  csp_D = D - 2 * yb;
  csp_H = H - 2 * zb;
}

void decompress(char *filename) {
  decomp_flag = 0;
  if (strstr(filename, ".gz")) {
    char command[100];
    int32_t ret;
    LogPrintf("decompressing CSP file\n");
    sprintf(command, "gunzip %s", filename);
    ret = system(command);
    if (ret) {
      ErrPrintf("ERROR: impossible to uncompress %s (ret=%d)\n", filename, ret);
      exit(-1);
    }
    filename[strlen(filename) - 3] = 0;
    decomp_flag = 1;
  }
}

void compress(char *filename, char force) {
  if (decomp_flag || force) {
    char command[100];
    int32_t ret;
    sprintf(command, "gzip %s &", filename);
    ret = system(command);
    if (ret) {
      ErrPrintf("ERROR: impossible to compress %s (ret=%d)\n", filename, ret);
      exit(-1);
    }
    decomp_flag = 0;
  }
}

int32_t read_csp_header(char *filename) {
  int32_t chunk, cksize, mdsize, offset, n;
  FILE *fp;

  decompress(filename);

  if (!strstr(filename, ".csp")) {
    ErrPrintf("ERROR: %s not a CSP file\n", filename);
    csp_header_size = 0;
    exit(-4);
  }

  assert(sizeof(int) == 4);

  LogPrintf("reading CSP metadata\n");

  if ((fp = fopen(filename, "r")) == NULL) {
    ErrPrintf("ERROR: cannot open file %s\n", filename);
    exit(-4);
  }

  // reset metadata
  memset(header, 0, sizeof(header));
  offset = 0;

  // magic number
  n = fread(header, 4, 1, fp);
  //LogPrintf("magic=%s   n=%d\n", header, n);
  if ((n != 1) || strcmp(header, CSP_MAGIC_NUM)) {
    ErrPrintf("ERROR: %s not in CSP format (bad magic number)\n", filename);
    exit(-1);
  }
  offset += 4;

  // endianness
  n = fread(header + offset, 4, 1, fp);
  if ((n != 1) || (header[offset + 2] != (little_endian() ? '_' : 'b'))) {
    ErrPrintf("ERROR: bad endianness (%c)\n", header[offset + 2]);
    exit(-1);
  }
  offset += 4;

  // chunk: size of header
  mdsize = sizeof(int);
  cksize = 8 + mdsize;
  n = fread(header + offset, cksize, 1, fp);
  if ((n != 1) || (READ_INT(header + offset + 4) != READ_INT("HDSZ"))) {
    ErrPrintf("ERROR: chunk HDSZ invalid\n");
    exit(-1);
  }
  csp_header_size = READ_INT(header + offset + 8);
  LogPrintf("size of header: %d\n", csp_header_size);
  offset += cksize;

  // reading all chunks left
  n = fread(header + offset, csp_header_size - offset, 1, fp);
  if (n != 1) {
    ErrPrintf("ERROR: failed to read header\n");
    exit(-1);
  }
  fclose(fp);

  while (offset < csp_header_size) {
    cksize = READ_INT(header + offset);
    chunk = READ_INT(header + offset + 4);
    offset += 8;
    assert(cksize >= 8);

    if (chunk == READ_INT("MODL")) {
      LogPrintf("chunk MODL: %s\n", header + offset);
      if (check_model && strcmp(header + offset, MOD_NAME)) {
        WarnPrintf("WARNING: wrong model (%s)\n", MOD_NAME);
      }
    } else if (chunk == READ_INT("SIZE")) {
      csp_H = READ_INT(header + offset);
      csp_L = READ_INT(header + offset + 4);
      csp_D = READ_INT(header + offset + 8);
      LogPrintf("chunk SIZE: %d x %d x %d\n", csp_H, csp_L, csp_D);
      if (((H != 0) && (csp_H != H)) || ((L != 0) && (csp_L != L)) || ((D != 0) && (csp_D != D))) {
        ErrPrintf("ERROR: wrong size (%d x %d x %d)\n", H, L, D);
        exit(-1);
      }
      H = csp_H;
      L = csp_L;
      D = csp_D;
    } else if (chunk == READ_INT("CELL")) {
      csp_cell_size = READ_INT(header + offset);
      LogPrintf("chunk CELL: %d\n", csp_cell_size);
      if (csp_cell_size > (int)sizeof(Cell)) {
        WarnPrintf("WARNING: too much data within cells (%d bytes), it will be truncated\n", csp_cell_size);
      }
    } else if (chunk == READ_INT("TIME")) {
      csp_time = READ_DOUBLE(header + offset);
      LogPrintf("chunk TIME: %f\n", csp_time);
    }

    offset += cksize - 8;
  }

  return csp_header_size;
}

void read_csp(char *filename) {
  FILE *fp;
  uint8_t *buf, *aux;
  int32_t i, j, k, n;
  int32_t szcell;
  Cell *pt = TE;

  LogPrintf("reading CSP data : %s\n", filename);

  if (!filename) {
    ErrPrintf("ERROR: read_csp - no filename\n");
    exit(-1);
  }

  if (!csp_cell_size) {
    csp_cell_size = sizeof(Cell);
  }

  AllocMemory(buf, unsigned char, HL * csp_cell_size);

  if ((fp = fopen(filename, "r")) == NULL) {
    ErrPrintf("ERROR: cannot open file %s\n", filename);
    exit(-4);
  }

  if (csp_header_size) {
    n = fread(header, csp_header_size, 1, fp);
    if (n != 1) {
      ErrPrintf("ERROR: read_csp - failed to read header\n");
      exit(-1);
    }
  }

  szcell = csp_cell_size ? Min(csp_cell_size, (int)sizeof(Cell)) : (int)sizeof(Cell);
  LogPrintf("szcell = %d\n", szcell);
  LogPrintf("csp_cell_size = %d\n", csp_cell_size);

  for (k = csp_y_bounds; k < D - csp_y_bounds; k++) { // profondeur
    fread(buf, csp_cell_size, csp_H * csp_L, fp);
    aux = buf;
    pt = TE + k * HL;
    for (j = 0; j < H; j++) { // hauteur
      for (i = 0; i < L; i++, pt++) { //largeur
        if (pt->celltype != BORD) {
          memcpy(pt, aux, szcell);
          if (*aux >= MAX_CELL) {
            ErrPrintf("ERROR : read_csp - incorrect cell type (%d)\n", *aux);
            exit(-1);
          }
          aux += csp_cell_size;
        }
      }
    }
  }
  fclose(fp);

  FreeMemory(buf, unsigned char, HL * csp_cell_size)

  compress(filename, 0);
}

int32_t generate_csp_header() {
  int32_t offset = 0;
  int32_t hdsz_offset = 0; //offset of header size
  int32_t cksize = 0; //chunk size
  int32_t mdsize = 0; //metadata size
  int32_t md4; //4-bytes metadata

  /// reset metadata
  memset(header, 0, sizeof(header));

  /// magic number
  sprintf(header, "\212CSP@1%c\n", little_endian() ? '_' : 'b');
  offset += 8;

  /// chunk: size of header
  mdsize = sizeof(int);
  cksize = 8 + mdsize;
  memcpy(header + offset, &cksize, 4);
  offset += 4;
  memcpy(header + offset, "HDSZ", 4);
  offset += 4;
  hdsz_offset = offset;
  memset(header + offset, 0, mdsize);
  offset += mdsize;

#ifdef MOD_NAME
  /// chunk: model name
  mdsize = strlen(MOD_NAME) + 1;
  Align4(mdsize); //4-bytes alignement
  cksize = 8 + mdsize;
  memcpy(header + offset, &cksize, 4);
  offset += 4;
  memcpy(header + offset, "MODL", 4);
  offset += 4;
  memcpy(header + offset, MOD_NAME, strlen(MOD_NAME) + 1);
  offset += mdsize;
#endif

  /// chunk: size of cellular space
  if (!csp_H && !csp_L && !csp_D) {
    csp_H = H;
    csp_L = L;
    csp_D = D;
  }
  mdsize = 3 * sizeof(int);
  cksize = 8 + mdsize;
  memcpy(header + offset, &cksize, 4);
  offset += 4;
  memcpy(header + offset, "SIZE", 4);
  offset += 4;
  memcpy(header + offset, &csp_H, sizeof(int));
  offset += sizeof(int);
  memcpy(header + offset, &csp_L, sizeof(int));
  offset += sizeof(int);
  memcpy(header + offset, &csp_D, sizeof(int));
  offset += sizeof(int);

  /// chunk: cells metadata
  md4 = sizeof(Cell);
  mdsize = sizeof(int);
  cksize = 8 + mdsize;
  memcpy(header + offset, &cksize, 4);
  offset += 4;
  memcpy(header + offset, "CELL", 4);
  offset += 4;
  memcpy(header + offset, &md4, mdsize);
  offset += mdsize;

  /// chunk: time
  mdsize =  sizeof(double);
  cksize = 8 + mdsize;
  memcpy(header + offset, &cksize, 4);
  offset += 4;
  memcpy(header + offset, "TIME", 4);
  offset += 4;
  memcpy(header + offset, &csp_time, sizeof(double));
  offset += mdsize;

  /// size of header
  memcpy(header + hdsz_offset, &offset, sizeof(int));

  return offset;
}

void write_csp(char dump_type, char *filename) {
  static int32_t hd_size = 0;
  static char start = 1;
  FILE *fp;
  int32_t i, j;
  Cell *pt = TE;
  Cell *buf, *aux;

  LogPrintf("writing CSP data : %s\n", filename);

  if (dump_type == DUMP_CSP) {
    // header with metadata
    hd_size = generate_csp_header();
    if (start) {
      LogPrintf("header size: %d\n", hd_size);
    }
    start = 0;
  }

  fp = fopen(filename, "w");
  if (!fp) {
    ErrPrintf("Erreur ouverture fichier CSP : %s\n", filename);
    exit(-4);
  }

  if (dump_type == DUMP_CSP) {
    fwrite(header, 1, hd_size, fp);
  }

  AllocMemory(buf, Cell, HL);
  ResetMemory(buf, Cell, HL);

  for (j = 0; j < D; j++) {
    aux = buf;
    for (i = 0; i < HL; i++, pt++) {
      if (pt->celltype != BORD) {
        *aux++ = *pt;
      }
    }
    if (aux > buf) {
      fwrite(buf, sizeof(Cell), csp_H * csp_L, fp);
    }
  }
  fclose(fp);
  FreeMemory(buf, Cell, HL);
}
