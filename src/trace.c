/* ReSCAL - Cells tracing
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

// KK TODO --- several #ifdef statements here say MODELE_DUN, MODELE_AVA, etc
// This looks like a typo which should say MODEL_DUN, MODEL_AVA, etc
// Not yet corrected because not sure if this is compensated for elsewhere.


//mode TRACE_SRC : generation du bassin versant d'un point32_t de sortie
//mode TRACE_AIRE : calcul de l'aire drainee (~debit) en chaque point32_t (2D)
//mode TRACE_PAR : calcul du debit et de la pente en chaque point32_t (2D), et sauvegarde du parcours effectif de chaque cellule
//mode TRACE_TRANS : comptage des transtions en chaque point32_t (2D)
//mode TRACE3D_CEL : comptage des cellules d'un type donne qui transitent en chaque point32_t (3D)

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <memory.h>
#include <stdint.h>

#include "defs.h"
#include "macros.h"
#include "trace.h"
#include "cells.h"
#include "doublets.h"
#include "transitions.h"
#include "space.h"
#include "surface.h"
#include "lgca.h"
#include "view.h"
#include "callbacks.h"
#include "simul.h"
#include "param.h"

extern int32_t H, L, D, HL, HLD;             // les dimensions de la terre
extern Cell *TE;                   // la 'terre'
extern double csp_time;                  // temps reel simule
extern float csp_angle;                  // angle resultant de toutes les rotations
extern int32_t L_bounds;
extern int32_t D_bounds;
extern int32_t rot_mode;
extern uint8_t reorient_flag;
extern int32_t vdir_mode;

int32_t source_x, source_y, source_z, source_ix;   //point32_t source
int32_t src_flag = 0; //indicateur

void trace_trans_init();
void trace3d_cel_init();
void trace3d_time_init();
void trace_flux_init();

void trace_dump_loop();

#define NB_CEL_SRC 5 //nombre de cellules sources sur chaque colonne (injectees une par une)

void trace_init() {

#ifdef TRACE_PAR_COL
  trace_par_col_init();
#endif

#ifdef TRACE_TRANS
  trace_trans_init();
#endif

#ifdef TRACE3D_CEL
  trace3d_cel_init();
#endif

#ifdef TRACE_FLUX
  trace_flux_init();
#endif
}

int32_t trace_init_loop() {
  static int32_t x = 0, y = 1, z = 1, cpt_src = 0;
  int32_t ix = 0;

  cpt_src++;
  if (cpt_src >= NB_CEL_SRC) {
    cpt_src = 0;
  }

  source_x = x;
  source_y = y;
  source_z = z;
  source_ix = ix;

  return 0;
}
void trace_end_loop() {
}

void trace_test(int32_t ix) {
  (void)ix; //SUPPRESS: unused warning
}

int32_t trace_point(int32_t x, int32_t z) {
  (void)x; //SUPPRESS: unused warning
  (void)z; //SUPPRESS: unused warning
  return 0;
}

void trace_dump_loop() {
}

///////////////////////////////////////////////////////////////////////////////



#ifdef TRACE_TRANS
/// On comptabilise les occurences de certaines transitions en fonction du modèle, en chaque colonne.
/// En sortie, on sauve une carte 2d pour chaque type de transition.

int32_t *trace_diffusion; //carte des transitions de diffusion
int32_t *trace_dissolution; //carte des transitions de dissolution
int32_t *trace_cristallisation; //carte des transitions de cristallisation
int32_t *trace_erosion; //carte des transitions d'erosion
int32_t *trace_chk_erosion; //carte des transitions d'erosion bloquees
int32_t *trace_transport_eo; //carte du transport est-ouest
int32_t *trace_transport; //carte du transport
int32_t *trace_gravi; //carte des transitions de gravite

void trace_trans_init() {
  LogPrintf("init trace_trans\n");

#if defined(MODEL_DUN) || defined(MODEL_SNO)
  AllocMemory(trace_transport_eo, int, L*D);
  ResetMemory(trace_transport_eo, int, L*D);
#elif defined(MODEL_AVA)
  AllocMemory(trace_transport, int, L*D);
  AllocMemory(trace_gravi, int, L*D);
  ResetMemory(trace_transport, int, L*D);
  ResetMemory(trace_gravi, int, L*D);
#else
  AllocMemory(trace_diffusion, int, L * D);
  AllocMemory(trace_dissolution, int, L * D);
  AllocMemory(trace_cristallisation, int, L * D);
  ResetMemory(trace_diffusion, int, L * D);
  ResetMemory(trace_dissolution, int, L * D);
  ResetMemory(trace_cristallisation, int, L * D);
#endif
}

void trace_trans(int32_t tr, int32_t ix) {
  int32_t x, y, z, ind;
  Calcule_xyz(ix, x, y, z);
  ind = x + z*L;
#if defined(MODEL_DUN) || defined(MODEL_SNO)
  if (tr == 4)
  {
    trace_transport_eo[ind]++;
  }
#else
  if (tr <= 3) {
    trace_diffusion[ind]++;
  } else if (tr <= 11) {
    trace_dissolution[ind]++;
  } else if (tr <= 19) {
    trace_cristallisation[ind]++;
  }
#endif
}

void trace_trans_blocked(int32_t tr, int32_t ix) {
  int32_t x, y, z, ind;
  Calcule_xyz(ix, x, y, z);
  ind = x + z * L;
#if defined(MODEL_DUN) || defined(MODEL_SNO)
  if (tr == 0)
  {
    trace_chk_erosion[ind]++;
  }
#endif
}

void trace_trans_dump() {
  FILE *fp;
  static int32_t cpt = 0;
  char nom[128];

#ifdef MODEL_DUN
  sprintf(nom, "TRANSPORT%03d.car", cpt);
  fp = fopen(nom, "wb");
  fwrite(trace_transport_eo, sizeof(int), L * D, fp);
  fclose(fp);
  memset(trace_transport_eo, 0, L * D * sizeof(int));
  cpt++;
#else
  fp = fopen("DIFFUS.car", "wb");
  fwrite(trace_diffusion, sizeof(int), L * D, fp);
  fclose(fp);

  fp = fopen("DISSOL.car", "wb");
  fwrite(trace_dissolution, sizeof(int), L * D, fp);
  fclose(fp);

  fp = fopen("CRISTA.car", "wb");
  fwrite(trace_cristallisation, sizeof(int), L * D, fp);
  fclose(fp);

  //remise a zero
  memset(trace_diffusion, 0, L * D * sizeof(int));
  memset(trace_dissolution, 0, L * D * sizeof(int));
  memset(trace_cristallisation, 0, L * D * sizeof(int));
#endif
}

void trace_trans_quit() {
#ifdef MODEL_DUN
  free(trace_transport_eo);
#else
  free(trace_diffusion);
  free(trace_dissolution);
  free(trace_cristallisation);
#endif
}
#endif  //TRACE_TRANS

///////////////////////////////////////////////////////////////////////////////

#ifdef TRACE_PAR_COL
/// On memorise le parcours 3d de toutes les cellules colorees (1 cellule par couleur).
/// La position est recalculee en cas de rotation, et on veille a l'integrite (non-disparition, unicite) de la cellule coloree.
/// En sortie, on sauve les donnees (x,y,z,t,transition) dans un fichier pour chaque parcours.

//constants
#define NB_CEL_PAR_COL_MAX 100 //max number of colors
#define LG_PAR_COL_MAX 10000  //max length of each course
#define BUF_SIZE_COL 10 //LG_PAR_COL_MAX  //size of buffer for each course

int32_t nb_par_col = 0; //number of colors

//structures
typedef struct un_element_de_parcours_color {
  int32_t x, y, z; //position
  int32_t tr; //transition
  float t;    //time
} par_elt_col;

par_elt_col parcours[NB_CEL_PAR_COL_MAX][BUF_SIZE_COL]; //courses of all colored cells

int32_t cur_elt[NB_CEL_PAR_COL_MAX];            //current indices for each course
int32_t lg_par[NB_CEL_PAR_COL_MAX];             //length for each course


void trace_par_col_init() {
  int32_t i, j, ix, cm1;
  par_elt_col *pe;

  LogPrintf("init trace_par_col\n");

  for (cm1 = 0; cm1 < NB_CEL_PAR_COL_MAX; cm1++) {
    cur_elt[cm1] = -1; //default value
  }

  //initialization of all courses
  for (ix = 0; ix < HLD; ix++) {
    if (TE[ix].color) {
      cm1 = TE[ix].color - 1;
      if (cm1 >= NB_CEL_PAR_COL_MAX) {
        ErrPrintf("ERROR: trace_par_col_init - color=%d not allowed\n", TE[ix].color);
        exit(-1);
      }
      if (cm1 >= nb_par_col) {
        nb_par_col = TE[ix].color;
      }
      cur_elt[cm1] = 0;
      //par_elt_col *pe = parcours[cm1];
      if (lg_par[cm1] == 0) {
        trace_par_col(ix, -1);
      } else {
        ErrPrintf("ERROR: trace_par_col_init - two cells with the same color(%d), this is not allowed in TRACE_PAR_COL mode !\n", TE[ix].color);
        exit(-1);
      }
    }
  }

  LogPrintf("%d color(s) detected\n", nb_par_col)
}

void trace_par_col_rotation(float angle) {
  int32_t ix, col, cur;
  par_elt_col *pe;
  int32_t x, y, z;

  //LogPrintf("trace_par_col_rotation - angle=%f\n",angle);

  //remove the colors
  for (ix = 0; ix < HLD; ix++) {
    TE[ix].color = 0;
  }

  //restore the colors
  for (col = 1; col <= NB_CEL_PAR_COL_MAX; col++) {
    cur = cur_elt[col - 1];
    if (cur >= 0) { //the color is present
      //LogPrintf("color %d is present\n", col);
      if (cur == 0) {
        cur = BUF_SIZE_COL;
      }
      //read the last position
      pe = &parcours[col - 1][cur - 1];
      //rotation
      x = pe->x;
      y = pe->y;
      z = pe->z;
      //LogPrintf("x=%d y=%d z=%d\n",x,y,z);
      if (angle) {
        Rotate_xz(x, z, angle);
      }
      //LogPrintf("x=%d y=%d z=%d\n",x,y,z);
      //new position
      ix = z * HL + y * L + x;
      //update the color
      if ((x >= 1) && (x < L - 1) && (y >= 1) && (y < H - 1) && (z >= 1) && (z < D - 1)) {
        if (TE[ix].color == 0) {
#ifdef MODEL_DUN
          if ((TE[ix].celltype == GR) || (TE[ix].celltype == GRJ)) {
            TE[ix].color = col;
          } else if (TE[ix].celltype == EAUC) {
            LogPrintf("WARNING: trace_par_col_rotation - not a grain (%d), GR added for color %d\n", TE[ix].celltype, col); /*exit(-1);*/
            TE[ix].celltype = GR; //adding a grain !
            TE[ix].color = col;
          } else {
            ErrPrintf("WARNING: trace_par_col_rotation - not a grain (%d), trace is stopped for color %d\n", TE[ix].celltype, col); /*exit(-1);*/
            cur_elt[col - 1] = -1; //trace is stopped for this grain
          }
#else
          TE[ix].color = col;
#endif
        } else {
          ErrPrintf("WARNING: trace_par_col_rotation - cell in %d is already colored (%d), cannot assign color %d\n", ix, TE[ix].color, col);
          cur_elt[col - 1] = -1; //trace is stopped for this grain
        }
      } else {
        ErrPrintf("WARNING: trace_par_col_rotation - colored (%d) grain is going outside\n", col);
        cur_elt[col - 1] = -1; //trace is stopped for this grain
      }
    }
  }
}

void dump_par_col(int32_t col)
{
  char name[100];
  int32_t j;
  char current_output[128];

  //save course for one color from current buffer
  sprintf(name,"PAR_COL%03d.data", col);
  for (j=0, pe = parcours[col-1]; j<cur_elt[col-1]; j++, pe++){
    sprintf(current_output, "%05d %05d %05d %03d %e\n", pe->x-1, pe->z-1, H-1-pe->y, pe->tr, pe->t);
    output_write(name, current_output);
  }
}

void trace_par_col(int32_t ix, int32_t tr) {
  if (TE[ix].color) {
    int32_t cm1 = TE[ix].color - 1;
    //LogPrintf("trace_par_col: ix=%d, tr=%d, col=%d, csp_angle=%f\n", ix, tr, TE[ix].color, csp_angle);
    if (cur_elt[cm1] < 0) {
      return;
    }
    if (lg_par[cm1] < LG_PAR_COL_MAX) {
      int32_t x, y, z;
      par_elt_col *pe = parcours[cm1] + cur_elt[cm1];
      Calcule_xyz(ix, x, y, z);
      //LogPrintf("trace_par_col: x=%d, y=%d, z=%d\n", x, y, z);
      if (csp_angle) {
        Rotate_xz(x, z, -csp_angle);
      }
      pe->x = x;
      pe->y = y;
      pe->z = z;
      pe->tr = tr;
      pe->t = (float)csp_time;
      cur_elt[cm1]++;
      lg_par[cm1]++;
      if (cur_elt[cm1] >= BUF_SIZE_COL) {
        //save course
        //ErrPrintf("WARNING: out of buffer par_col (col=%d, size=%d) !\n", TE[ix].color, BUF_SIZE_COL);
        dump_par_col(TE[ix].color);
        cur_elt[cm1] = 0;
      }
    } else {
      //course too long
      ErrPrintf("ERROR: cell course is too int64_t (col=%d, length=%d) !\n", TE[ix].color, lg_par[cm1]);
      //color is removed
      TE[ix].color = 0;
    }
  }
}


#endif  //TRACE_PAR_COL

///////////////////////////////////////////////////////////////////////////////

#ifdef TRACE3D_CEL
/// Comptage des cellules de type CEL_TEST qui transitent en chaque point32_t de l'espace cellulaire (3d).
/// En sortie, on obtient l'intensite du flux.
/// Option TRACE_MVT : comptage des cellules de gaz sur reseau en chaque point32_t et dans toutes les directions.

//constantes
#ifdef MODEL_DUN
//#define CEL_TEST  EAUT
#define CEL_TEST  GRJ
#endif

#ifdef LGCA
//#define TRACE_MVT
#endif

int32_t *cpt_cel; //compteur de cellules (en chaque point)
#ifdef TRACE_MVT
int32_t *cpt_mvt_est; //compteur de cellules mvt_est (en chaque point)
int32_t *cpt_mvt_ouest; //compteur de cellules mvt_ouest (en chaque point)
int32_t *cpt_check_mvt_est; //compteur de tests reussis pour le controle du transport vers l'est (en chaque point)
int32_t *cpt_check_mvt_bas; //compteur de tests reussis pour le controle du transport vers le bas (en chaque point)
#endif //TRACE_MVT
int32_t nb_trace3d_cel = 0; // nombre d'appels de la procedure trace3d_cel()

void trace3d_cel_init() {
  LogPrintf("init trace3d_cel\n");

  AllocMemory(cpt_cel, int, HLD);
  ResetMemory(cpt_cel, int, HLD);

#ifdef TRACE_MVT
  AllocMemory(cpt_mvt_est, int, HLD);
  ResetMemory(cpt_mvt_est, int, HLD);

  AllocMemory(cpt_mvt_ouest, int, HLD);
  ResetMemory(cpt_mvt_ouest, int, HLD);

  AllocMemory(cpt_check_mvt_est, int, HLD);
  ResetMemory(cpt_check_mvt_est, int, HLD);

  AllocMemory(cpt_check_mvt_bas, int, HLD);
  ResetMemory(cpt_check_mvt_bas, int, HLD);
#endif //TRACE_MVT
}

void trace3d_cel() {
  int32_t i;
  for (i = 0; i < HLD; i++) {
    if (TE[i].celltype == CEL_TEST) {
      cpt_cel[i]++;
    }
#ifdef TRACE_MVT
    if (mvt_est(i)) {
      cpt_mvt_est[i]++;
    }
    if (mvt_ouest(i)) {
      cpt_mvt_ouest[i]++;
    }
    if (check_mvt_est(i)) {
      cpt_check_mvt_est[i]++;
    }
    if (check_mvt_bas(i)) {
      cpt_check_mvt_bas[i]++;
    }
#endif //TRACE_MVT
  }

  nb_trace3d_cel++;
}

void trace3d_cel_dump() {
  FILE *fp;
  static int32_t cpt = 0;
  //char nom[128];

  //sprintf(nom, "CEL%03d.bin", cpt++);
  fp = fopen("CEL.bin", "wb");
  fwrite(cpt_cel, sizeof(int), HLD, fp);
  fclose(fp);

#ifdef TRACE_MVT
  fp = fopen("MVT_EST.bin", "wb");
  fwrite(cpt_mvt_est, sizeof(int), HLD, fp);
  fclose(fp);

  fp = fopen("MVT_OUEST.bin", "wb");
  fwrite(cpt_mvt_ouest, sizeof(int), HLD, fp);
  fclose(fp);

  fp = fopen("CHK_MVT_EST.bin", "wb");
  fwrite(cpt_check_mvt_est, sizeof(int), HLD, fp);
  fclose(fp);

  fp = fopen("CHK_MVT_BAS.bin", "wb");
  fwrite(cpt_check_mvt_bas, sizeof(int), HLD, fp);
  fclose(fp);
#endif //TRACE_MVT

  fp = fopen("TRACE3D_CEL.data", "w");
  fprintf(fp, "%d\n", nb_trace3d_cel);
  fclose(fp);

  //memset(cpt_cel, 0, HLD*sizeof(int));
  cpt++;
}

void trace3d_cel_quit() {
  free(cpt_cel);
#ifdef TRACE_MVT
  free(cpt_mvt_est);
  free(cpt_mvt_ouest);
  free(cpt_check_mvt_est);
  free(cpt_check_mvt_bas);
#endif //TRACE_MVT
}

#endif //TRACE3D_CEL

///////////////////////////////////////////////////////////////////////////////

#ifdef TRACE_FLUX
/// On calcule la direction et l'intensite du flux moyen en chaque colonne a partir des transitions de transport.
/// On prend en compte les eventuelles rotations dans le calcul de la direction de transport.
/// En sortie, on sauve les cartes 2d contenant les composantes est-ouest et nord-sud du flux.

extern Doublet t_doub[];
extern TransitionDb t_trans[];

FluxMap *flux_map = NULL;
FluxMap *total_flux_map = NULL;
float trace_flux_delay = 0;
int32_t trace_flux_auto = 0;

FluxMap *init_flux_map(float angle) {
  FluxMap *map;

  AllocMemory(map, FluxMap, 1);

  map->angle = angle;

  AllocMemory(map->flux_ew, float, L * D);
  ResetMemory(map->flux_ew, float, L * D);

  AllocMemory(map->flux_ns, float, L * D);
  ResetMemory(map->flux_ns, float, L * D);

  return map;
}

void reset_flux_map(FluxMap *map) {
  //map->angle = 0;
  ResetMemory(map->flux_ew, float, L * D);
  ResetMemory(map->flux_ns, float, L * D);
}

void trace_flux_init() {
  LogPrintf("init trace_flux: angle=%f\n", csp_angle);

  if (trace_flux_delay && trace_flux_auto) {
    ErrPrintf("ERROR: parameters Trace_flux_delay and Trace_flux_auto cannot be used simultaneously. Please modify your parameter file.\n");
    exit(-1);
  }

  flux_map = init_flux_map(csp_angle);
  total_flux_map = init_flux_map(0);
}


void trace_flux(int32_t tr, int32_t ix) {
  TransitionDb *ptr;
  int32_t x, y, z;
  int32_t ind;

  //LogPrintf("flux dir = %d\n", ptr->classe);

  /// calcul position
  Calcule_xyz(ix, x, y, z);

  /// sens du transport
  int32_t step = 1;
  ptr = &t_trans[tr];
  int32_t type_cel1=t_doub[ptr->depart].one;
#if defined(MODEL_DUN) || defined(MODEL_SNO)
  if (type_cel1 == EAUC){
    step = -1;
  }
#endif

  /// orientation du transport et sommation
  ind = x + z * L;
  if (ptr->classe == EST_OUEST) {
    flux_map->flux_ew[ind] += step;
  } else {
    flux_map->flux_ns[ind] += step;
  }
}

void compute_total_flux() {
  int32_t x, z, ind, ind0;
  float x0, z0, tx, tz, f00, f10, f01, f11;

  //LogPrintf("compute_total_flux: %f\n", flux_map->angle);

  if (!flux_map->angle) { /// no rotation
    for (z = 0, ind = 0; z < D; z++) {
      for (x = 0; x < L; x++, ind++) {
        total_flux_map->flux_ew[ind] += flux_map->flux_ew[ind];
        total_flux_map->flux_ns[ind] += flux_map->flux_ns[ind];
      }
    }
  } else { /// rotated map
    float alpha = flux_map->angle * PI / 180.0;
    float co = cos(alpha);
    float si = sin(alpha);
    float cx = L / 2.0 - 0.5;
    float cz = D / 2.0 - 0.5;
    float dx, dz, flx0, flz0;

    for (z = 0; z < D; z++) {
      for (x = 0; x < L; x++) {
        dx = x - cx;
        dz = z - cz;
        /// inverse rotation
        x0 = cx + co * dx - si * dz;
        if ((x0 < 0) || (x0 > L - 2)) {
          continue;
        }
        z0 = cz + si * dx + co * dz;
        if ((z0 < 0) || (z0 > D - 2)) {
          continue;
        }
        ind0 = floorf(z0) * L + floorf(x0);
        /// perform 2d linear interpolation
        tx = ceilf(x0) - x0;
        tz = ceilf(z0) - z0;
        f00 = tx * tz;
        f10 = tz - f00;
        f01 = tx - f00;
        f11 = 1 - f10 - f01 - f00;
        ind = z * L + x;
        if ((ind < 0) || (ind >= L * D)) {
          ErrPrintf("ERROR: compute_total_flux - ind=%d\n", ind);
          exit(-1);
        }
        flx0 = f00 * flux_map->flux_ew[ind0] + f10 * flux_map->flux_ew[ind0 + 1] + f01 * flux_map->flux_ew[ind0 + L] + f11 * flux_map->flux_ew[ind0 + L + 1];
        flz0 = f00 * flux_map->flux_ns[ind0] + f10 * flux_map->flux_ns[ind0 + 1] + f01 * flux_map->flux_ns[ind0 + L] + f11 * flux_map->flux_ns[ind0 + L + 1];
        /// rotation and integration into the resultant flux
        total_flux_map->flux_ew[ind] += co * flx0 + si * flz0;
        total_flux_map->flux_ns[ind] += -si * flx0 + co * flz0;
      }
    }
  }
}

void trace_flux_rotate(float angle) {
  //LogPrintf("trace_flux_rotate: angle=%f, time=%f\n", angle, csp_time);

  if (trace_flux_auto && (csp_time > 0)) {
    trace_flux_dump();
  } else {
    compute_total_flux();
  }

  reset_flux_map(flux_map);

  flux_map->angle = angle;
}

void trace_flux_dump()
{
  static int32_t cpt = 1;
  char name[128];
  int32_t x, z;
  char current_output[128];

  LogPrintf("trace_flux_dump: cpt = %d, csp_time = %f\n", cpt, csp_time);

  compute_total_flux();

  /// save east-west total flux
  sprintf(name, "FLUX_EW_%04d", cpt);
  for(z=D_bounds; z<D-D_bounds; z++){
    for(x=L_bounds; x<L-L_bounds; x++){
      sprintf(current_output,"%.3f ",total_flux_map->flux_ew[x+z*L]);
      output_write(name, current_output);
    }
    output_write(name, "\n");
  }

  /// save north-south total flux
  sprintf(name, "FLUX_NS_%04d.data", cpt);
  for(z=D_bounds; z<D-D_bounds; z++){
    for(x=L_bounds; x<L-L_bounds; x++){
      sprintf(current_output,"%.3f ",total_flux_map->flux_ns[x+z*L]);
      output_write(name, current_output);
    }
    output_write(name, "\n");
  }

#ifdef ALTI
  /// reorient and save topography
  if (csp_angle) {
    lock_display(0);
    reorient_flag = 1;
    rotation(0, rot_mode, ROT_REORIENT_TEMP);
    //display wind direction
    vdir_mode = VDIR_WIND;
  }

  dump_surface("FLUX_AL_", cpt, UNIT_COMP);

  if (reorient_flag) {
    rotation(0, rot_mode, ROT_REORIENT_UNDO);
    vdir_mode = VDIR_NONE;
    reorient_flag = 0;
    unlock_display(0);
  }
#endif

  //reset all maps
  reset_flux_map(flux_map);
  reset_flux_map(total_flux_map);

  cpt++;
}
#endif  //TRACE_FLUX

///////////////////////////////////////////////////////////////////////////////

/// Parameters declaration for the trace module

void params_trace() {
#ifdef TRACE_FLUX
  parameter("Trace_flux_delay", "delay between outputs of flux data (t0 unit) - optional", &trace_flux_delay, PARAM_FLOAT, "");
  parameter("Trace_flux_auto", "automatic outputs of flux data for each wind (0|1) - optional", &trace_flux_auto, PARAM_INT, "");
#endif
}


/// Data outputs during current simulation

void trace_dump(char flag_info) {
  (void)flag_info; //SUPPRESS: unused warning
#ifdef TRACE_TRANS
  if (flag_info) {
    trace_trans_dump();
  }
#endif

#ifdef TRACE3D_CEL
  if (flag_info) {
    trace3d_cel_dump();
  }
#endif
}


