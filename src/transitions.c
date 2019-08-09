/* ReSCAL - States transitions
 *
 * Copyright (C) 2011
 *
 * Author: Olivier Rozier <rozier@ipgp.fr>
 *
 * This file is part of ReSCAL.
 *
 * Code based on dissol program,
 * by Eduardo Sepulveda <edo@espci.fr>
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
#include <string.h>
#include <math.h>
#include <stdint.h>

#include "defs.h"
#include "macros.h"
#include "space.h"
#include "surface.h"
#include "cells.h"
#include "doublets.h"
#include "transitions.h"
#include "models.h"
#include "trace.h"
#include "simul.h" //output

extern const char       *etats[];                                        // les noms des types de cellules
extern Doublet          t_doub[];                                        // la description des doublets
extern int32_t          nb_type_db;                                      // nombre de types de doublets
extern Cell             *TE;                                             // la 'terre'
extern int32_t          H, L, D, HL, HLD;                                // les dimensions de la terre
extern Doublet          t_doub[];                                        // la description des doublets
extern const char       *classes_db[];                                   // noms des classes de doublet
extern int32_t          LNS;                                             // largeur nord-sud
extern int32_t          pbc_mode;                                        // periodic boundary conditions
extern char             *rot_map;                                        // periodic mapping of the rotating space
extern Pos2             *rot_map_pos;                                    // periodic mapping of the rotating space
extern int32_t          *db_pos[];                                       // les tableaux contenant la position des doublets actifs

TransitionDb            t_trans[MAX_TRANSITIONS_DB];                     // la description des transitions de 2 cellules
TransitionCel           t_trans_cel[MAX_TRANSITIONS_CEL];                // la description des transitions d'1 cellule
int32_t                 nb_trans_db = 0, nb_trans_cel = 0, nb_trans = 0; // nombre de transitions de 2 cellules, d'1 cellule et au total
LienTransDb             t_lien[MAX_LIENS];                               // la description des liens de correlation entre transitions de 2 cellules
int32_t                 nb_liens = 0;                                    // nombre de liens entre transitions
int32_t                 cpt_lien[MAX_LIENS];                             // compteur des correlations effectives
int32_t                 cpt_trans_db[MAX_TRANSITIONS_DB];                // compteur de transitions effectives de doublet
int32_t                 cpt_trans_cel[MAX_TRANSITIONS_CEL];              // compteur de transitions effectives de cellule
int32_t                 direction[3];                                    // conversion classe de transition -> direction
int32_t                 orientation[6];

void split_link_hor(int32_t tr0, int32_t tr);

void init_transitions() {
  int32_t i;
  // determination de la matrice inverse des doublets (et decouplage des doublets et transitions)
  init_db_inv();
  nb_trans = nb_trans_db + nb_trans_cel;
  LogPrintf("nombre de transitions = %d\n", nb_trans);

  // initialisation des compteurs de transitions
  for (i = 0; i < nb_trans_db; i++) {
    cpt_trans_db[i] = 0;
  }
  for (i = 0; i < nb_trans_cel; i++) {
    cpt_trans_cel[i] = 0;
  }

  // initialisation des compteurs de correlations
  for (i = 0; i < nb_liens; i++) {
    cpt_lien[i] = 0;
  }

  // initialisation tableau direction
  direction[VERTICAL] = BAS;
  direction[EST_OUEST] = EST;
  direction[NORD_SUD] = SUD;

  // initialisation tableau orientation
  orientation[BAS] = VERTICAL;
  orientation[HAUT] = VERTICAL;
  orientation[EST] = EST_OUEST;
  orientation[OUEST] = EST_OUEST;
  orientation[SUD] = NORD_SUD;
  orientation[NORD] = NORD_SUD;
}

int32_t get_trans(char ref) {
  int32_t i;

  i = 0;
  while ((i < nb_trans_db) && (t_trans[i].reference != ref)) {
    i++;
  }
  if (i >= nb_trans_db) {
    ErrPrintf("ERROR: get_trans - reference does not exist !\n");
    exit(-1);
  }
  return i;
}

void trans_ref(char ref, char classe, char cel_depart_1, char cel_depart_2, char cel_arrivee_1, char cel_arrivee_2, double intensite) {
  if (intensite < 0) {
    LogPrintf("WARNING: negative transition rate (%f). Transition removed.\n", intensite);
    return;
  }

  if (classe == ISOTROPE){
    // transition isotrope
    // KK -- this ISOTROPE clause creates a set of vertical and horizontal conditions which have the
    // net effect of an isotropic (every-direction) transition.
    char sym = (cel_depart_1 == cel_depart_2) && (cel_arrivee_1 == cel_arrivee_2);
    if (H > 3) {
      trans_ref(-1, VERTICAL, cel_depart_1, cel_depart_2, cel_arrivee_1, cel_arrivee_2, intensite);
      if (!sym) {
        trans_ref(-1, VERTICAL, cel_depart_2, cel_depart_1, cel_arrivee_2, cel_arrivee_1, intensite);
      }
    }
    trans_ref(-1, HORIZONTAL, cel_depart_1, cel_depart_2, cel_arrivee_1, cel_arrivee_2, intensite);
    if (!sym) {
      trans_ref(-1, HORIZONTAL, cel_depart_2, cel_depart_1, cel_arrivee_2, cel_arrivee_1, intensite);
    }
    return;
  }

  if ((classe == HORIZONTAL) && (LNS == 1)) {
    classe = EST_OUEST;  //2D
  }

  if (intensite || (ref != -1)) {
    char cel1_flag = (char)(cel_depart_1 != cel_arrivee_1);
    char cel2_flag = (char)(cel_depart_2 != cel_arrivee_2);
    char transport_flag = (cel_depart_1 == cel_arrivee_2) && (cel_depart_2 == cel_arrivee_1);
    char inout_flag = (cel_depart_1 == BORD) || (cel_depart_2 == BORD) || (cel_depart_1 == IN) || (cel_depart_2 == IN) || (cel_depart_1 == OUT) || (cel_depart_2 == OUT);
    // definition des transitions de 2 cellules :
    // reference = numero de reference (optionnel) pour identification
    // classe    = la classe de transition: VERTICAL ou HORIZONTAL
    // depart    = type de doublet de depart
    // arrivee   = type de doublet d'arrivee
    // intensite = intensite de la transition
    // type      = type de la transition
    // cel1_flag = flag (0,1) qui indique s'il y a changement d'etat pour la premiere cellule
    // cel2_flag = flag (0,1) qui indique s'il y a changemant d'etat pour la deuxieme cellule
    // lien       = flag (0,1) qui indique si la transition est le depart d'un
    //              lien de correlation vers une autre transition
    // regul     = pointeur callback de regulation de l'intensite
    // nb_chk    = nombre de callbacks de controle
    // checks    = tableau des donnees de controle
    // time_mode = mode d'evolution du TEMPS (TIME_EVOL | TIME_CORR | TIME_NO_EVOL), surtout pour les transitions avec controle
    t_trans[nb_trans_db].reference = ref;
    t_trans[nb_trans_db].classe = classe;
    t_trans[nb_trans_db].depart = init_doublet(classe, cel_depart_1, cel_depart_2, 1);
    t_trans[nb_trans_db].arrivee = init_doublet(classe, cel_arrivee_1, cel_arrivee_2, 0);
    t_trans[nb_trans_db].intensite = intensite;
    t_trans[nb_trans_db].type = (cel1_flag && transport_flag) ? TR_TRANSPORT : (inout_flag) ? TR_INOUT : TR_NOTYPE;
    t_trans[nb_trans_db].cel1_flag = cel1_flag;
    t_trans[nb_trans_db].cel2_flag = cel2_flag;
    t_trans[nb_trans_db].lien = 0;
    t_trans[nb_trans_db].regul = NULL;
    t_trans[nb_trans_db].nb_chk = 0;
    t_trans[nb_trans_db].checks = NULL;
    t_trans[nb_trans_db].time_mode = TIME_EVOL;
    nb_trans_db++;
  }
} //trans_ref

void trans(char classe, char cel_depart_1, char cel_depart_2, char cel_arrivee_1, char cel_arrivee_2, double intensite) {
  trans_ref(-1, classe, cel_depart_1, cel_depart_2, cel_arrivee_1, cel_arrivee_2, intensite);
}

void trans_type(char ref, char type) {
  int32_t i;
  i = get_trans(ref);
  t_trans[i].type = type;
}

void trans_time(char ref, char mode) {
  int32_t i;
  i = get_trans(ref);
  t_trans[i].time_mode = mode;
}

void trans_regul(char ref, Callback_regul reg_func) {
  int32_t i;
  i = get_trans(ref);
  t_trans[i].regul = reg_func;
}

int32_t trans_check_flag(char ref, Callback_check chk_func, char cel, char inv) {
  DataCheck *pdc;
  int32_t i, j;

  i = get_trans(ref);
  if (cel != 1 && cel != 2) {
    ErrPrintf("WARNING : trans_check - incorrect cell parameter: %d (admissible values: 1,2)\n", cel);
  }
  if (!t_trans[i].checks) {
    //memory allocation
    AllocMemory(t_trans[i].checks, DataCheck, MAX_CHK);
  }
  j = t_trans[i].nb_chk++;
  if (j >= MAX_CHK) {
    ErrPrintf("ERROR: trans_check - too many controls added on transition ref=%d\n", ref);
    exit(-1);
  };
  pdc = &t_trans[i].checks[j];
  pdc->func = chk_func;
  pdc->cel = cel;
  pdc->inv = inv;
#ifdef CGV
  t_trans[i].time_mode = ((chk_func == check_grad_vel) && (t_trans[i].nb_chk == 1)) ? TIME_CORR : TIME_NO_EVOL;
#else
  t_trans[i].time_mode = TIME_NO_EVOL;
#endif

  return j;
}

int32_t trans_check(char ref, Callback_check chk_func, char cel) {
  return trans_check_flag(ref, chk_func, cel, 0);
}

int32_t trans_check_inv(char ref, Callback_check chk_func, char cel) {
  return trans_check_flag(ref, chk_func, cel, 1);
}

void trans_check_cell(char ref, char cel, char dir, char cell_type) {
  DataCheck *pdc;
  int32_t i, j;

  /// find transition
  i = get_trans(ref);

  /// declare callback
  j = trans_check(ref, check_cell_dir, cel);

  /// fill data structure attached to the callback reference
  pdc = t_trans[i].checks + j;
  pdc->char_data1 = dir;
  pdc->char_data2 = cell_type;
}

void trans_check_no_cell(char ref, char cel, char dir, char cell_type) {
  DataCheck *pdc;
  int32_t i, j;

  /// find transition
  i = get_trans(ref);

  /// declare callback
  j = trans_check_inv(ref, check_cell_dir, cel);

  /// fill data structure attached to the callback reference
  pdc = t_trans[i].checks + j;
  pdc->char_data1 = dir;
  pdc->char_data2 = cell_type;
}

#ifdef CELL_COLOR
extern Callback_check check_color;

void trans_check_color(char ref, char cel, float lambda_col) {
  DataCheck *pdc;
  float lambda;
  int32_t i, j;

  /// find transition
  i = get_trans(ref);
  lambda = t_trans[i].intensite;

  if (lambda == lambda_col) {
    LogPrintf("WARNING: trans_check_color - identical lambda values (%f,%f), no control added\n", lambda, lambda_col);
    return;
  }

  /// declare callback
  j = trans_check(ref, check_color, cel);

  /// fill data structure attached to the callback reference
  pdc = t_trans[i].checks + j;
  if (lambda > lambda_col) {
    //control applied on colored cells
    pdc->col = 1;
    pdc->coef = lambda_col / lambda;
  } else {
    //control applied on black cells
    pdc->col = 0;
    pdc->coef = lambda / lambda_col;
    //intensity rate must be increased
    LogPrintf("transition rate increased for colored cells: ref=%d   lambda=%f\n", ref, lambda_col);
    t_trans[i].intensite = lambda_col;
  }
}
#endif

int32_t split_trans_hor(int32_t db_hor) {
  int32_t i, depart_eo, depart_ns, arrivee_eo, arrivee_ns, flag = 0;

  // decouplage des transitions contenant le doublet horizontal db_hor
  for (i = 0; i < nb_trans_db; i++) {
    if ((t_trans[i].depart == db_hor) || (t_trans[i].arrivee == db_hor)) {
      LogPrintf("decouplage de la transition horizontale %d [ %s, %s ] -> [ %s, %s ]\n", i,
                etats[t_doub[t_trans[i].depart].one], etats[t_doub[t_trans[i].depart].two],
                etats[t_doub[t_trans[i].arrivee].one], etats[t_doub[t_trans[i].arrivee].two]);

      //decouplage du doublet de depart
      split_db_hor(t_trans[i].depart, &depart_eo, &depart_ns);
      //decouplage du doublet d'arrivee
      split_db_hor(t_trans[i].arrivee, &arrivee_eo, &arrivee_ns);
      // mise-a-jour de la transition HORIZONTAL -> EST_OUEST
      t_trans[i].classe = EST_OUEST;
      t_trans[i].depart = depart_eo;
      t_trans[i].arrivee = arrivee_eo;
      // ajout d'une nouvelle transition NORD_SUD
      t_trans[nb_trans_db] = t_trans[i];
      t_trans[nb_trans_db].classe = NORD_SUD;
      t_trans[nb_trans_db].depart = depart_ns;
      t_trans[nb_trans_db].arrivee = arrivee_ns;
      // decouplage des liens
      split_link_hor(i, nb_trans_db);
      nb_trans_db++;
      flag = 1;
    }
  }
  return flag;
}

void add_link(int32_t tr1, int32_t tr2, char ltr_cel, double intensite) {
  if (nb_liens >= MAX_LIENS) {
    ErrPrintf("ERROR: max number of links (%d) exceeded\n", MAX_LIENS);
    exit(-1);
  }

  if (intensite > 0) {
    t_trans[tr1].lien = 1;
    if ((ltr_cel != 1) && (ltr_cel != 2)) {
      ErrPrintf("ERROR: ltr_cel is incorrect : %d\n", ltr_cel);
      exit(-1);
    }
    t_lien[nb_liens].cel = ltr_cel;
    t_lien[nb_liens].trans1 = tr1;

    t_lien[nb_liens].trans2 = tr2;

    if ((t_doub[t_trans[tr1].arrivee].one != t_doub[t_trans[tr2].depart].one) &&
        (t_doub[t_trans[tr1].arrivee].one != t_doub[t_trans[tr2].depart].two) &&
        (t_doub[t_trans[tr1].arrivee].two != t_doub[t_trans[tr2].depart].one) &&
        (t_doub[t_trans[tr1].arrivee].two != t_doub[t_trans[tr2].depart].two)) {
      ErrPrintf("ERROR: link between transitions %d and %d is impossible\n", tr1, tr2);
      exit(-1);
    }

    t_lien[nb_liens].intensite = intensite;

    nb_liens++;
  }
}

void trans_link(char ref1, char ref2, char ltr_cel, double intensite) {
  // definition des liens entre 2 transitions de 2 cellules :
  // ref1 = reference de la premiere transition
  // ref2 = reference de la deuxieme transition
  // ltr_cel = indique si le lien s'effectue sur la premiere ou la deuxieme cellule du doublet
  // intensite = intensite de la correlation entre les 2 transitions
  add_link(get_trans(ref1), get_trans(ref2), ltr_cel, intensite);
}

void split_link_hor(int32_t tr0, int32_t tr) {
  // pour chaque lien avec tr0, creer un lien identique avec tr a la place de tr0
  int32_t i;
  LogPrintf("decouplage des liens contenant la transition %d\n", tr0);
  for (i = 0; i < nb_liens; i++) {
    if (t_lien[i].trans1 == tr0) {
      add_link(tr, t_lien[i].trans2, t_lien[i].cel, t_lien[i].intensite);
    }
    if (t_lien[i].trans2 == tr0) {
      add_link(t_lien[i].trans1, tr, t_lien[i].cel, t_lien[i].intensite);
    }
  }
}

void trans_cel_ref(char ref, char cel_depart, char cel_arrivee, double intensite) {
  // definition des transitions d'une cellule
  // reference = numero de reference (optionnel) pour identification
  // depart    = etat de depart
  // arrivee   = etat d'arrivee
  // intensite = intensite de la transition
  // func = pointeur callback
  t_trans_cel[nb_trans_cel].reference = ref;
  t_trans_cel[nb_trans_cel].depart = cel_depart;
  t_trans_cel[nb_trans_cel].arrivee = cel_arrivee;
  t_trans_cel[nb_trans_cel].intensite = intensite;
  t_trans_cel[nb_trans_cel].regul = NULL;
  nb_trans_cel++;
}

void trans_cel(char cel_depart, char cel_arrivee, double intensite) {
  trans_cel_ref(-1, cel_depart, cel_arrivee, intensite);
}

void trans_cel_regul(char ref, Callback_regul f) {
  int32_t i;
  i = get_trans(ref);
  t_trans_cel[i].regul = f;
}

void do_trans_db(int32_t tr, int32_t ix, int32_t dir) {
  int32_t db_depart, db_arrivee;
  int32_t ix2, tc, tc2;

  if (tr >= nb_trans_db) {
    ErrPrintf("ERROR: do_trans_db - transition %d does not exist\n", tr);
    exit(-1);
  }

  db_depart = t_trans[tr].depart;
  db_arrivee = t_trans[tr].arrivee;
  tc = t_doub[db_arrivee].one;
  tc2 = t_doub[db_arrivee].two;
  ix2 = 0;

  switch (dir) {

  case BAS:

    elimine_doublet(db_depart, ix, BAS);

    if (t_trans[tr].cel1_flag) {   // partie haute est modifiee

      elimine_doublet_haut(ix);
      elimine_doublet_ouest(ix);
      elimine_doublet_est(ix);
      elimine_doublet_nord(ix);
      elimine_doublet_sud(ix);

      // modifier cellule haute
      modifie_cellule(tc, ix);

      ajoute_doublet_haut(ix);
      ajoute_doublet_ouest(ix);
      ajoute_doublet_est(ix);
      ajoute_doublet_nord(ix);
      ajoute_doublet_sud(ix);

    } // FIN modif partie haute

    if (t_trans[tr].cel2_flag) {   // DEBUT modif partie basse
      ix2 = get_cell_down(ix);

      elimine_doublet_ouest(ix2);
      elimine_doublet_est(ix2);
      elimine_doublet_nord(ix2);
      elimine_doublet_sud(ix2);
      elimine_doublet_bas(ix2);

      // modifier cellule basse
      modifie_cellule(tc2, ix2);

      ajoute_doublet_ouest(ix2);
      ajoute_doublet_est(ix2);
      ajoute_doublet_nord(ix2);
      ajoute_doublet_sud(ix2);
      ajoute_doublet_bas(ix2);
    } // FIN  modif partie basse

    // traitement du doublet
    if (t_doub[db_arrivee].actif) {   // verifier si le doublet d'arrivee est actif
      ajoute_doublet(db_arrivee, ix, BAS);
    }

    break;  // FIN BAS

  case EST:

    elimine_doublet(db_depart, ix, EST);

    if (t_trans[tr].cel1_flag) {   // DEBUT partie ouest

      elimine_doublet_ouest(ix);
      elimine_doublet_haut(ix);
      elimine_doublet_bas(ix);
      elimine_doublet_nord(ix);
      elimine_doublet_sud(ix);

      // modifier cellule ouest
      modifie_cellule(tc, ix);

      ajoute_doublet_ouest(ix);
      ajoute_doublet_haut(ix);
      ajoute_doublet_bas(ix);
      ajoute_doublet_nord(ix);
      ajoute_doublet_sud(ix);

    } // FIN modification partie ouest

    if (t_trans[tr].cel2_flag) {   // DEBUT partie est

      ix2 = get_cell_east(ix);

      elimine_doublet_haut(ix2);
      elimine_doublet_bas(ix2);
      elimine_doublet_nord(ix2);
      elimine_doublet_sud(ix2);
      elimine_doublet_est(ix2);

      // modifier cellule est
      modifie_cellule(tc2, ix2);

      ajoute_doublet_haut(ix2);
      ajoute_doublet_bas(ix2);
      ajoute_doublet_nord(ix2);
      ajoute_doublet_sud(ix2);
      ajoute_doublet_est(ix2);

    } // FIN  partie est

    if (t_doub[db_arrivee].actif) {   // verifier si le doublet d'arrivee est actif
      ajoute_doublet(db_arrivee, ix, EST);
    } else { // en cas de doublet horizontal actif
      ajoute_doublet_est(ix);
    }

    break;  // FIN EST

  case SUD:
    elimine_doublet(db_depart, ix, SUD);

    if (t_trans[tr].cel1_flag) {   // DEBUT partie nord

      elimine_doublet_nord(ix);
      elimine_doublet_ouest(ix);
      elimine_doublet_est(ix);
      elimine_doublet_haut(ix);;
      elimine_doublet_bas(ix);

      // modifier cellule nord
      modifie_cellule(tc, ix);

      ajoute_doublet_nord(ix);
      ajoute_doublet_ouest(ix);
      ajoute_doublet_est(ix);
      ajoute_doublet_haut(ix);;
      ajoute_doublet_bas(ix);

    } // FIN modification partie nord

    if (t_trans[tr].cel2_flag) {   // DEBUT partie sud

      ix2 = get_cell_south(ix);

      elimine_doublet_ouest(ix2);
      elimine_doublet_est(ix2);
      elimine_doublet_haut(ix2);
      elimine_doublet_bas(ix2);
      elimine_doublet_sud(ix2);

      // modifier cellule sud
      modifie_cellule(tc2, ix2);

      ajoute_doublet_ouest(ix2);
      ajoute_doublet_est(ix2);
      ajoute_doublet_haut(ix2);
      ajoute_doublet_bas(ix2);
      ajoute_doublet_sud(ix2);

    } // FIN  partie sud

    if (t_doub[db_arrivee].actif) {   // verifier si le doublet d'arrivee est actif
      ajoute_doublet(db_arrivee, ix, SUD);
    } else { // en cas de doublet horizontal actif
      ajoute_doublet_sud(ix);
    }

    break; // FIN SUD

  default:
    break;
    // FIN DES TRANSITIONS
  } // SWITCH

#if CELL_DATA
  if ((t_trans[tr].type == TR_TRANSPORT) || (t_trans[tr].type == TR_INOUT)) {
    swap_cell_data(ix, ix2);
  }
  if (t_trans[tr].type == TR_INOUT) {
    update_inout_data(ix, ix2);
  }
#endif

#ifdef INFO_TRANS
  cpt_trans_db[tr]++;
#endif

#ifdef TRACE_PAR_COL
  if (t_trans[tr].type == TR_TRANSPORT) {
    trace_par_col(ix, tr);
    if (ix2) {
      trace_par_col(ix2, tr);
    }
  }
#endif
}

void do_trans_cel(int32_t tr, int32_t ix) {
  if (tr >= nb_trans_cel) {
    ErrPrintf("ERROR: cell transition %d does not exist\n", tr);
    exit(-1);
  }

  cellule_terre(t_trans_cel[tr].arrivee, ix);

#ifdef INFO_TRANS
  cpt_trans_cel[tr]++;
#endif
}

void dump_transitions()
{
  int32_t i, j;
  static char flag_trans[200];
  static char type_trans[100];
  char current_output[256];

  sprintf(current_output, "\n# TRANSITIONS\n\nNB_TRANS_DB = %d\n", nb_trans_db);
  output_write("TRANSITIONS", current_output);

  for(i=0; i<nb_trans_db; i++){
    *flag_trans = 0;
    if (t_trans[i].regul){
      strcat(flag_trans,"regul");
    }
    if (t_trans[i].checks) {
      for (j = 0; j < t_trans[i].nb_chk; j++) {
        if (*flag_trans) {
          strcat(flag_trans, ",");
        }
        if (t_trans[i].checks[j].cel == 1) {
          strcat(flag_trans, "chk[1]");
        } else {
          strcat(flag_trans, "chk[2]");
        }
        if (t_trans[i].checks[j].inv) {
          strcat(flag_trans, "-");
        } else {
          strcat(flag_trans, "+");
        }
      }
    }
    if (t_trans[i].time_mode != TIME_EVOL) {
      if (*flag_trans) {
        strcat(flag_trans, ",");
      }
      if (t_trans[i].time_mode == TIME_CORR) {
        strcat(flag_trans, "time_corr");
      } else {
        strcat(flag_trans, "time_no_evol");  //TIME_NO_EVOL
      }
    }

    *type_trans = 0;
    sprintf(current_output,"TR(%2d): %s, [%s, %s] -> [%s, %s], %.9f", i, classes_db[t_trans[i].classe],
                  etats[t_doub[t_trans[i].depart].one], etats[t_doub[t_trans[i].depart].two],
                  etats[t_doub[t_trans[i].arrivee].one], etats[t_doub[t_trans[i].arrivee].two],
                  t_trans[i].intensite);
    if (*flag_trans){
      sprintf(current_output, " {%s}", flag_trans);
      output_write("TRANSITIONS", current_output);
    }
    if (*type_trans){
      sprintf(current_output, " <%s>", type_trans);
      output_write("TRANSITIONS", current_output);
    }
    output_write("TRANSITIONS", "\n");
  }
  sprintf(current_output,"\nNB_TRANS_CEL = %d\n",nb_trans_cel);
  output_write("TRANSITIONS", current_output);
  for(i=0; i<nb_trans_cel; i++){
    sprintf(current_output,"TR1(%2d) : %s -> %s, %.9f\n", i, etats[t_trans_cel[i].depart], etats[t_trans_cel[i].arrivee], t_trans_cel[i].intensite);
    output_write("TRANSITIONS", current_output);
  }
#ifdef LIENS_TRANSITIONS
  sprintf(current_output,"\n# LINKED TRANSITIONS\n");
  output_write("TRANSITIONS", current_output);
  sprintf(current_output, "\nNB_LINKS = %d\n", nb_liens);
  output_write("TRANSITIONS", current_output);
  for(i=0; i<nb_liens; i++)
    sprintf(current_output,"LT(%2d): %d[%d] -> %d, %.9f\n", i, t_lien[i].trans1, t_lien[i].cel, t_lien[i].trans2, t_lien[i].intensite);
  output_write("TRANSITIONS", current_output);
#endif
}

#ifdef INFO_TRANS
extern int32_t cpt_trans_blocked[];

void dump_trans_info_header(){
  char current_output[256];

  // Write some initial information about transitions
  output_write("TRANSITIONS", "\n# NUMBER OF OCCURRENCES\n\n    ");
  for (int32_t i=0; i<nb_trans_db; i++){
    sprintf(current_output, "     TR(%2d)", i);
    if (t_trans[i].nb_chk) strcat(current_output, "       BLK");
    output_write("TRANSITIONS", current_output);
  }
  for (int32_t i=0; i<nb_trans_cel; i++){
    sprintf(current_output,"    TR1(%2d)", i);
    output_write("TRANSITIONS", current_output);
  }
  for (int32_t i=0; i<nb_liens; i++){
    sprintf(current_output, "     LT(%2d)", i);
    output_write("TRANSITIONS", current_output);
  }
}
  

void dump_trans_info(){
  int32_t i;
  char current_output[256];
  int32_t cpt = 1;

  sprintf(current_output,"\n%04d:",cpt);
  output_write("TRANSITIONS", current_output);
  for (i=0; i<nb_trans_db; i++){
    sprintf(current_output," %10d", cpt_trans_db[i]);
    output_write("TRANSITIONS", current_output);
    cpt_trans_db[i] = 0;
    if (t_trans[i].nb_chk){
      sprintf(current_output," %10d",cpt_trans_blocked[i]);
      output_write("TRANSITIONS", current_output);
      cpt_trans_blocked[i] = 0;
    }
  }
  for (i=0; i<nb_trans_cel; i++){
    sprintf(current_output," %10d", cpt_trans_cel[i]);
    cpt_trans_cel[i] = 0;
    output_write("TRANSITIONS", current_output);
  }
  for(i=0; i<nb_liens; i++){
    sprintf(current_output," %10d", cpt_lien[i]);
    output_write("TRANSITIONS", current_output);
    cpt_lien[i] = 0;
  }
  cpt++;
}
#endif //INFO_TRANS
