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
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */



enum TRANS_TYPES {TR_NOTYPE, TR_TRANSPORT, TR_INOUT};

enum TIME_MODES {TIME_EVOL, TIME_CORR, TIME_NO_EVOL};

typedef int Callback_regul(void *data);
typedef int Callback_check(int ix, void *data);
typedef struct data_check DataCheck;

typedef struct une_transition_de_doublet{
  char reference;
  unsigned char classe; //horizontale, verticale ou isotrope
  unsigned char depart, arrivee;  //doublets de depart et d'arrivee
  double intensite;
  char type;  //type de transition
  char cel1_flag; //premiere cellule modifiee
  char cel2_flag; //deuxiemme cellule modifiee
  char lien; //indicateur "lien de transitions"
  Callback_regul *regul; //fonction de regulation globale de l'intensite
  int nb_chk; //nombre de fonctions de controle local
  DataCheck *checks; //tableau des donnees de controle local
  //Callback_check **checks;  //tableau des fonctions de controle local (inhibition si l'une d'elles retourne 0)
  char time_mode; //indique si le temps evolue et s'il doit etre corrige
} TransitionDb;

struct data_check{
  //int trans;   //transition
  Callback_check *func;
  char dir; //direction
  char cel; //premiere ou deuxieme cellule
  char inv; //inversion flag
#ifdef CELL_COLOR
  char col;
#endif
  float coef;
  char char_data1; //generic data
  char char_data2; //generic data
  //void *data;
};

typedef struct une_transition_de_cellule{
  char reference;
  unsigned char depart, arrivee;
  double intensite;
  Callback_regul *regul; //fonction de regulation globale de l'intensite
} TransitionCel;

typedef struct un_lien_entre_deux_transitions_de_doublets{
  int trans1, trans2;
  char cel; //indique si le lien s'effectue sur la premiere ou la deuxieme cellule du doublet
  double intensite;
} LienTransDb;


void init_transitions();
void trans_ref(char ref, char classe, char cel_depart_1, char cel_depart_2, char cel_arrivee_1, char cel_arrivee_2, double intensite);
void trans(char classe, char cel_depart_1, char cel_depart_2, char cel_arrivee_1, char cel_arrivee_2, double intensite);
void trans_type(char ref, char type);
void trans_time(char ref, char mode);
void trans_regul(char ref, Callback_regul f);
int trans_check(char ref, Callback_check chk_func, char cel);
int trans_check_inv(char ref, Callback_check chk_func, char cel);
void trans_check_cell(char ref, char cel, char dir, char cel_type);
void trans_check_no_cell(char ref, char cel, char dir, char cel_type);
#ifdef CELL_COLOR
void trans_check_color(char ref, char cel, float lambda_col);
#endif
void trans_link(char ref1, char ref2, char ltr_cel, double intensite);
int split_trans_hor(int db_hor);
void trans_cel_ref(char ref, char cel_depart, char cel_arrivee, double intensite);
void trans_cel(char cel_depart, char cel_arrivee, double intensite);
void trans_cel_regul(char ref, Callback_regul f);
void do_trans_db(int tr, int ix, int dir);
void do_trans_cel(int tr, int ix);
void dump_transitions();
void dump_trans_info();

//for the compatibility with the old syntax
#define lien_trans trans_link
