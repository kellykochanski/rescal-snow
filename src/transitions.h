/* ReSCAL - States transitions
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

#include <stdint.h>

enum TRANS_TYPES {TR_NOTYPE, TR_TRANSPORT, TR_INOUT};

enum TIME_MODES {TIME_EVOL, TIME_CORR, TIME_NO_EVOL};

typedef int32_t Callback_regul(void *data);
typedef int32_t Callback_check(int32_t ix, void *data);
typedef struct data_check DataCheck;

typedef struct une_transition_de_doublet{
  int8_t reference;
  uint8_t classe; //horizontale, verticale ou isotrope
  uint8_t depart, arrivee;  //doublets de depart et d'arrivee
  double intensite;
  int8_t type;  //type de transition
  int8_t cel1_flag; //premiere cellule modifiee
  int8_t cel2_flag; //deuxiemme cellule modifiee
  int8_t lien; //indicateur "lien de transitions"
  Callback_regul *regul; //fonction de regulation globale de l'intensite
  int32_t nb_chk; //nombre de fonctions de controle local
  DataCheck *checks; //tableau des donnees de controle local
  //Callback_check **checks;  //tableau des fonctions de controle local (inhibition si l'une d'elles retourne 0)
  int8_t time_mode; //indique si le temps evolue et s'il doit etre corrige
} TransitionDb;

struct data_check{
  //int32_t trans;   //transition
  Callback_check *func;
  int8_t dir; //direction
  int8_t cel; //premiere ou deuxieme cellule
  int8_t inv; //inversion flag
#ifdef CELL_COLOR
  int8_t col;
#endif
  float coef;
  int8_t char_data1; //generic data
  int8_t char_data2; //generic data
  //void *data;
};

typedef struct une_transition_de_cellule{
  int8_t reference;
  uint8_t depart, arrivee;
  double intensite;
  Callback_regul *regul; //fonction de regulation globale de l'intensite
} TransitionCel;

typedef struct un_lien_entre_deux_transitions_de_doublets{
  int32_t trans1, trans2;
  int8_t cel; //indique si le lien s'effectue sur la premiere ou la deuxieme cellule du doublet
  double intensite;
} LienTransDb;


void init_transitions();
void trans_ref(int8_t ref, int8_t classe, int8_t cel_depart_1, int8_t cel_depart_2, int8_t cel_arrivee_1, int8_t cel_arrivee_2, double intensite);
void trans(int8_t classe, int8_t cel_depart_1, int8_t cel_depart_2, int8_t cel_arrivee_1, int8_t cel_arrivee_2, double intensite);
void trans_type(int8_t ref, int8_t type);
void trans_time(int8_t ref, int8_t mode);
void trans_regul(int8_t ref, Callback_regul f);
int32_t trans_check(int8_t ref, Callback_check chk_func, int8_t cel);
int32_t trans_check_inv(int8_t ref, Callback_check chk_func, int8_t cel);
void trans_check_cell(int8_t ref, int8_t cel, int8_t dir, int8_t cel_type);
void trans_check_no_cell(int8_t ref, int8_t cel, int8_t dir, int8_t cel_type);
#ifdef CELL_COLOR
void trans_check_color(int8_t ref, int8_t cel, float lambda_col);
#endif
void trans_link(int8_t ref1, int8_t ref2, int8_t ltr_cel, double intensite);
int32_t split_trans_hor(int32_t db_hor);
void trans_cel_ref(int8_t ref, int8_t cel_depart, int8_t cel_arrivee, double intensite);
void trans_cel(int8_t cel_depart, int8_t cel_arrivee, double intensite);
void trans_cel_regul(int8_t ref, Callback_regul f);
void do_trans_db(int32_t tr, int32_t ix, int32_t dir);
void do_trans_cel(int32_t tr, int32_t ix);
void dump_transitions();
void dump_trans_info();

//for the compatibility with the old syntax
#define lien_trans trans_link
