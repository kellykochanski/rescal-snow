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
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

#define CELL_DATA defined(CELL_COLOR)

#ifndef IN
#define IN BORD
#endif

#ifndef OUT
#define OUT BORD
#endif


typedef struct cellule{
  unsigned char celltype;
#ifdef CELL_COLOR
  unsigned char color;
#endif
#ifdef CELL_TIME
  int celltime; //date du dernier changement d'etat
#endif
} Cell;

#define CELL_TYPE_SIZE sizeof(unsigned char)
#define CELL_DATA_SIZE (sizeof(Cell) - CELL_TYPE_SIZE)

//typedef char CellData[CELL_DATA_SIZE];


void init_Ncel();
void init_cellule(Cell cel, int index);
void modifie_cellule(int type, int index);
//void ajoute_cellule(int type, int index);
void deplace_cellule(int ix, int ix2);
//void elimine_cellule(int index);
//void log_cell_info();
void log_cell();
