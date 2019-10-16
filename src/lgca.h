/* ReSCAL - Lattice gas
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

#ifdef LGCA

/// forcing of the flow
//#define MVT_REGUL  //regulation flux
//#define MVT_HFOR  //forcage horizontal du flux (tapis roulant)
//#define MVT_VFOR  //forcage partiel du flux sur le bord vertical gauche
#define MVT_GFOR  //forcage global du flux

/// solid-fluid interface
/// to set no-slip conditions, comment the four lines below
//#define SOLSLOW   //collisions realistes entre particules solides et particules lentes
//#define SOLFAST   //collisions realistes entre particules solides et particules rapides
//#define SF_REAL   //collisions realistes entre particules solides et particules fluides (via les normales)
#define PLAF_REAL   //collisions realistes (avec glissement) au plafond uniquement

/// sliding window for the velocity averaging
//#define VEL_SLIDE

//densite moyenne des particules en mouvement
#define DENSITE 1.5 //valeur <= 2.0

//masques pour les mouvements de fluides sur reseaux (modele multi-speed)

#define MVT_SOLID 1 //phase solide
#define MVT_OUT 2 //sortie
#define MVT_E 4   //est
#define MVT_O 8   //ouest
#define MVT_B 16  //bas
#define MVT_H 32  //haut
#define MVT_EB 64 //est-bas
#define MVT_EH 128  //est-haut
#define MVT_OB 256  //ouest-bas
#define MVT_OH 512  //ouest-haut

#define MVT_SLOW (MVT_E | MVT_O | MVT_B | MVT_H) //vitesses lentes
#define MVT_FAST (MVT_EB | MVT_EH | MVT_OB | MVT_OH) //vitesses rapides
#define MVT_ALLDIR (MVT_SLOW | MVT_FAST) //toute direction

typedef uint16_t MvtField;

#define SIZE_MVT_FIELD 1024   //(1<<(8*sizeof(MvtField)))

//parametres pour le couplage et le moyennage

#define NB_MVT_VER 1 //5  //nombre de noeuds verticalement par cellules
#define NB_MVT_EO 1 //5  //nombre de noeuds est-ouest par cellules
#ifndef STABILITY_ANALYSIS
#define DIST_MVT_NS 5 //5   //distance nord-sud entre 2 plans verticaux de collisions (par defaut)
#else
#define DIST_MVT_NS 2    //distance nord-sud entre 2 plans verticaux de collisions (analyse de stabilite)
#endif

#define VSTEP 5    //pas pour le moyennage spatial (ATTENTION : VSTEP doit etre un multiple de NB_MVT_EO et NB_MVT_VER !!)
//#define VSTEP_H (VSTEP*NB_MVT_VER)       //pas pour le moyennage spatial vertical
//#define VSTEP_L (VSTEP*NB_MVT_HOR)       //pas pour le moyennage spatial horizontal
#define VSTEP_H VSTEP    //pas pour le moyennage spatial vertical
#define VSTEP_L VSTEP    //pas pour le moyennage spatial est-ouest
#define VSTEP_TIME 50 //100  //pas pour le moyennage temporel
#define VMEAN_REF (VSTEP_TIME*VSTEP_TIME)/100 //seuil de vitesse moyenne (carre de la norme)

//macros

#define SetMask(mvt, mask) mvt |= mask
#define UnsetMask(mvt, mask) mvt &= (-1) ^ mask

#define Calcule_cix(index, cix) \
{ \
  int32_t z = (int)index/HL; \
  int32_t xy = index - z*HL; \
  int32_t cz = (int)(z - LN)/DIST_MVT_NS; \
  cix = cz*CHL + xy*NB_MVT_EO; \
}

#define Calcule_cxyz(index, cx, cy, cz) \
{                           \
  int32_t x,y,z;                  \
  Calcule_xyz(index, x, y, z); \
  cx = x*NB_MVT_EO; \
  cy = y; \
  cz = (int)(z-LN)/DIST_MVT_NS; \
}

void params_collisions();
void init_collisions();
void init_mvt();
void out_of_space_mvt(int32_t reset_mvt);
void collisions_mod_terre();
int32_t mvt(int32_t ii);
int32_t mvt_solid(int32_t ii);
int32_t mvt_bas(int32_t ii);
int32_t mvt_haut(int32_t ii);
int32_t mvt_est(int32_t ii);
int32_t mvt_ouest(int32_t ii);
//Callback_check check_mvt_bas;
//Callback_check check_mvt_haut;
//Callback_check check_mvt_ouest;
//Callback_check check_mvt_est;
int32_t check_mvt_solid(int32_t index);
int32_t check_mvt_bas(int32_t index);
int32_t check_mvt_haut(int32_t index);
int32_t check_mvt_est(int32_t index);
int32_t check_mvt_ouest(int32_t index);
int32_t check_mvt_est_et_bas(int32_t index);
int32_t check_mvt_EB(int32_t index);
int32_t check_no_mvt_est_et_bas(int32_t index);
int32_t check_mvt_mean_bas(int32_t index);
int32_t check_mvt_mean_est(int32_t index);
int32_t check_no_mvt_mean_bas(int32_t index);
void collisions_modcell(int32_t type, int32_t index);
void do_collisions();
void do_propagations();
void compute_vel(int8_t flag_interp);
void dump_mvt(int32_t cpt, int32_t unit);
void dump_densite();
void dump_vel();
void dump_signature_mvt();
void dump_collisions();

#endif //LGCA
