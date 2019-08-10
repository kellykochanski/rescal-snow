/* ReSCAL - Cells
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

#include <stdint.h>

#define CELL_DATA defined(CELL_COLOR)

#ifndef IN
#define IN BORD
#endif

#ifndef OUT
#define OUT BORD
#endif


typedef struct cellule {
  uint8_t celltype;
#ifdef CELL_COLOR
  uint8_t color;
#endif
#ifdef CELL_TIME
  int32_t celltime; //date du dernier changement d'etat
#endif
} Cell;

#define CELL_TYPE_SIZE sizeof(unsigned char)
#define CELL_DATA_SIZE (sizeof(Cell) - CELL_TYPE_SIZE)

//typedef int8_t CellData[CELL_DATA_SIZE];


void init_Ncel();
void init_cellule(Cell cel, int32_t index);
void modifie_cellule(int32_t type, int32_t index);
//void ajoute_cellule(int32_t type, int32_t index);
void deplace_cellule(int32_t ix, int32_t ix2);
//void elimine_cellule(int32_t index);
//void log_cell_info();
void log_cell();
