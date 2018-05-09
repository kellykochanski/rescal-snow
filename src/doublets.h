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
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

#define DB_INACTIF -1
#define DB_BORD -2
#define DB_TUNNEL -3

typedef struct un_doublet{
  unsigned char classe, one, two, actif;
} Doublet;

//choix du codage des references aux doublets actifs, stockes dans des tableaux
//#define REFDB_PTR

#ifdef REFDB_PTR

//la reference vers les doublets actifs est un pointeur (32 ou 64 bits selon l'architecture de la machine cible)
typedef int *RefDoublets[3];

#else

//la reference vers les doublets actifs est codee sur un caractere (8 bits) pour le type, et un entier (32 bits) pour l'index dans le tableau de positions
typedef char RefDoublets_Type[3];
typedef int RefDoublets_Ind[3];

#endif

int init_doublet(int classe, char etat1, char etat2, char actif);
void split_db_hor(int db_hor, int *db_eo, int *db_ns);
int type_doublet(int index, int dir);
void elimine_doublet(int type, int index, int dir);
void elimine_doublet_est(int index);
void elimine_doublet_ouest(int index);
void elimine_doublet_bas(int index);
void elimine_doublet_haut(int index);
void elimine_doublet_sud(int index);
void elimine_doublet_nord(int index);
void ajoute_doublet(int type, int index, int dir);
void ajoute_doublet_est(int index);
void ajoute_doublet_ouest(int index);
void ajoute_doublet_bas(int index);
void ajoute_doublet_haut(int index);
void ajoute_doublet_sud(int index);
void ajoute_doublet_nord(int index);
void init_db_inv();
void init_db_pos();
void fin_db_pos();
void dump_db_info();
void dump_doublets();


