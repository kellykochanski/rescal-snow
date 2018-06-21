/* ReSCAL - Doublets
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

#define DB_INACTIF -1
#define DB_BORD -2
#define DB_TUNNEL -3

typedef struct un_doublet {
  uint8_t classe, one, two, actif;
} Doublet;

//choix du codage des references aux doublets actifs, stockes dans des tableaux
//#define REFDB_PTR

#ifdef REFDB_PTR
//la reference vers les doublets actifs est un pointeur (32 ou 64 bits selon l'architecture de la machine cible)
typedef int32_t *RefDoublets[3];
#else
//la reference vers les doublets actifs est codee sur un caractere (8 bits) pour le type, et un entier (32 bits) pour l'index dans le tableau de positions
typedef char RefDoublets_Type[3];
typedef int32_t RefDoublets_Ind[3];
#endif

int32_t init_doublet(int32_t classe, char etat1, char etat2, char actif);
void split_db_hor(int32_t db_hor, int32_t *db_eo, int32_t *db_ns);
int32_t type_doublet(int32_t index, int32_t dir);
void elimine_doublet(int32_t type, int32_t index, int32_t dir);
void elimine_doublet_est(int32_t index);
void elimine_doublet_ouest(int32_t index);
void elimine_doublet_bas(int32_t index);
void elimine_doublet_haut(int32_t index);
void elimine_doublet_sud(int32_t index);
void elimine_doublet_nord(int32_t index);
void ajoute_doublet(int32_t type, int32_t index, int32_t dir);
void ajoute_doublet_est(int32_t index);
void ajoute_doublet_ouest(int32_t index);
void ajoute_doublet_bas(int32_t index);
void ajoute_doublet_haut(int32_t index);
void ajoute_doublet_sud(int32_t index);
void ajoute_doublet_nord(int32_t index);
void init_db_inv();
void init_db_pos();
void fin_db_pos();
void dump_db_info();
void dump_doublets();


