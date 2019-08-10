/* ReSCAL - Surface processes
 *
 * Copyright (C) 2011-2013
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
#include <math.h>
#include <assert.h>
#include <stdint.h>
#include "defs.h"
#include "macros.h"
#include "param.h"
#include "space.h"
#include "doublets.h"
#include "cells.h"
#include "surface.h"
#include "lgca.h"
#include "transitions.h"
#include "simul.h"
extern double csp_time;                  // temps reel simule
extern int32_t H, L, D, HL, HLD;       // les dimensions de la terre
extern int32_t LN, LS, LEO, LNS, HLN;    //couloir est-ouest (limite nord, limite sud, largeur nord-sud, ...)
extern Cell  *TE;            // la 'terre'
extern int32_t pbc_mode;  //periodic boundary conditions
extern int32_t ava_mode;  //mode d'avalanches
extern int32_t ava_h_lim;
extern float ava_delay; //delai entre avalanches
extern float ava_angle;
extern float ava_angle_stable;
extern float ava_angle_col;
extern int32_t ava_upwind;
extern int32_t nb_pv;                  // nombre de plans verticaux (non-DUM)
extern char *psol;                    // indicateur des plans solides verticaux est-ouest
extern char *rot_map;        // periodic mapping of the rotating space
extern Pos2 *rot_map_pos;        // periodic mapping of the rotating space
#ifdef PHASES
const uint8_t Phase[MAX_CELL] = PHASES;   //phase (fluide ou solide) des types de cellules
#endif // PHASES
#ifdef ALTI
int16_t *alti = NULL;           // elevations locales de terrain
// elevations locales de terrain, après moyennage
Vec2 *norm2d = NULL;          // normale 2d a la surface definie par alti
Vec3 *norm3d = NULL;          // normale 3d a la surface definie par alti
Vec3 *norm3d_iso = NULL;      // normale 3d a la surface definie par alti (avec tangentes isotropes)
float angle_norm3d_iso = 0.0; // direction variable pour le calcul des normales 3d (avec tangentes isotropes)
float *grdv = NULL;           // gradient de vitesse en surface
float slope_angle = 20.0;     // pente maximale pour la montee des grains
float ava_propag_angle = 30.0; // pente minimale pour les avalanches avec propagation
int32_t nb_ava_propag = 0;        // nombre de cellules deplacees par propagation d'avalanches
int32_t nb_drop_ava_propag = 0;   // nombre de cellules qui tombent par propagation d'avalanches
char *ava_mask = NULL;        // tableau de localisation des avalanches
FifoPos2* ava_fifo = NULL;     // file ordonnee (circulaire) des avalanches
int32_t veg_h_max = 0;       // hauteur max pour la vegetation
void loop_ava_propag(int32_t i, int32_t j);

void output_write(char* output_filename, char* output_content);
void output_initialize();

void params_surface() {
}
// creation d'une file FIFO
FifoPos2* fifo_create(int32_t len) {
  FifoPos2* fifo;
  AllocMemory(fifo, FifoPos2, 1);
  ResetMemory(fifo, FifoPos2, 1);
  fifo->length = len;
  AllocMemoryPrint("fifo->array", fifo->array, Pos2, len);
  ResetMemory(fifo->array, Pos2, len);
  return fifo;
}
// ajout d'un element dans la file
void fifo_put(FifoPos2 *fifo, Pos2 elt) {
  assert(fifo->nb < fifo->length);
  fifo->array[fifo->tail++] = elt;
  if (fifo->tail >= fifo->length) {
    fifo->tail = 0;
  }
  fifo->nb++;
}
// enlevement d'un element dans la file
Pos2 fifo_get(FifoPos2 *fifo) {
  Pos2 elt;
  assert(fifo->nb > 0);
  elt = fifo->array[fifo->head++];
  if (fifo->head >= fifo->length) {
    fifo->head = 0;
  }
  fifo->nb--;
  return elt;
}
// modification de l'altitude par ajout d'une cellule de type 'typ', en ix
// (code optimise)
void modif_alti_cel(int32_t ix, uint8_t typ) {
  int32_t x, y, z, i, j;
  int32_t al;
  Calcule_xyz(ix, x, y, z);
  i = x - 1;
  j = z - LN;
  al = Alti(i, j);
  if (Phase[typ] == SOLID) {
    if (H - y == al + 1) {
      Cell *aux;
      aux = TE + ix - L;
      al++;
      while ((al < H - 2) && (Phase[aux->celltype] == SOLID)) {
        al++;  //on remonte a la surface
        aux -= L;
      }
      Alti(i, j) = al;
    }
  } else {
    if (H - y <= al) {
      Alti(i, j) = H - y - 1;  //disparition d'une cellule dans la pile
    }
  }
}
void calcule_alti(uint8_t typ, char alti_mode) {
  (void)typ; //SUPPRESS unused warning
  Cell *aux;
  short *pt;
  int i, j, k;

  if (!alti) {
    AllocMemoryPrint("alti", alti, short, L * D/*LEO*LNS*/);
  }

  //calcul des elevations locales
  pt = alti;
  if (alti_mode == ALTI_MODE_BAS) { //mode bas : on part d'en bas et on remonte jusqu'a trouver une cellule non-solide
    for (j = 0; j < LNS; j++)
      for (i = 0; i < LEO; i++, pt++) {
        k = H - 2;
        aux = TE + (1 + i) + k * L + (LN + j) * HL;
        while ((k > 0) && (Phase[aux->celltype] == SOLID)) {
          k--;
          aux -= L;
        }
        *pt = H - (k + 1);
      }
  } else { //mode haut : on part d'en haut et on descend jusqu'a trouver une cellule solide
    for (j = 0; j < LNS; j++)
      for (i = 0; i < LEO; i++, pt++) {
        k = 1;
        aux = TE + (1 + i) + k * L + (LN + j) * HL;
        while ((k < H) && (Phase[aux->celltype] != SOLID)) {
          k++;
          aux += L;
        }
        *pt = H - k;
      }
  }
}
int32_t calcule_alti_max(uint8_t typ) {
  Cell *aux;
  int16_t max = 0;
  int32_t i, j, k;
  for (j = 0; j < LNS; j++)
    for (i = 0; i < LEO; i++) {
      k = 1;
      aux = TE + (1 + i) + k * L + (LN + j) * HL;
      while ((k < H - 1) && (aux->celltype != typ)) {
        k++;
        aux += L;
      }
      if (max < H - k) {
        max = H - k;
      }
    }
  return max;
}
#ifndef LGCA
#define NB_MVT_EO 1
#define NB_MVT_VER 1
#endif
// apply boundary conditions on one cell
void apply_bc_cell(int32_t *pt_x, int32_t *pt_y) {
  if (!pbc_mode) {
    /// open or closed boundaries
    if (*pt_x < 0) {
      *pt_x = 0;
    }
    if (*pt_x >= LEO) {
      *pt_x = LEO - 1;
    }
    if (*pt_y < 0) {
      *pt_y = 0;
    }
    if (*pt_y >= LNS) {
      *pt_y = LNS - 1;
    }
  } else {
    /// periodic boundary conditions
    if (!rot_map) {
      if (*pt_x < 0) {
        *pt_x += LEO;
      }
      if (*pt_x >= LEO) {
        *pt_x -= LEO;
      }
      if (*pt_y < 0) {
        *pt_y += LNS;
      }
      if (*pt_y >= LNS) {
        *pt_y -= LNS;
      }
    } else {
      /// rotating space
      int32_t x = *pt_x + 1;
      int32_t y = *pt_y + LN;
      assert((x >= 0) && (y >= 0));
      if (OutOfSpace(x, y)) {
        *pt_x = RotMapPosX(x, y) - 1;
        *pt_y = RotMapPosY(x, y) - LN;
      }
    }
  }
}
// apply boundary conditions on interpolated cell with float coordinates
void apply_bc_cel_interp(float *pt_x, float *pt_y) {
  if (!pbc_mode) {
    /// open or closed boundaries
    if (*pt_x < 0) {
      *pt_x = 0;
    }
    if (*pt_x >= LEO) {
      *pt_x = LEO - 1;
    }
    if (*pt_y < 0) {
      *pt_y = 0;
    }
    if (*pt_y >= LNS) {
      *pt_y = LNS - 1;
    }
  } else {
    /// periodic boundary conditions
    if (!rot_map) {
      if (*pt_x < 0) {
        *pt_x += LEO;
      }
      if (*pt_x >= LEO) {
        *pt_x -= LEO;
      }
      if (*pt_y < 0) {
        *pt_y += LNS;
      }
      if (*pt_y >= LNS) {
        *pt_y -= LNS;
      }
    } else {
      /// rotating space
      int32_t x = floor(*pt_x) + 1;
      int32_t y = floor(*pt_y) + LN;
      if (OutOfSpace(x, y)) {
        *pt_x += RotMapPosX(x, y) - x;
        *pt_y += RotMapPosY(x, y) - y;
      }
    }
  }
}
//interpolation de l'altitude en un point32_t intermediaire
float interpole_alti_cell(float x, float z) {
  int32_t x0, x1, z0, z1;
  int32_t y00, y01, y10, y11;
  float tx0, tx1, tz0, tz1;
  x0 = floor(x);
  z0 = floor(z);
  x1 = x0 + 1;
  z1 = z0 + 1;
  tx1 = x - x0;
  tx0 = 1.0 - tx1;
  tz1 = z - z0;
  tz0 = 1.0 - tz1;
  y00 = Alti(x0, z0);
  apply_bc_cell(&x0, &z1);
  y01 = Alti(x0, z1);
  apply_bc_cell(&x1, &z0);
  y10 = Alti(x1, z0);
  apply_bc_cell(&x1, &z1);
  y11 = Alti(x1, z1);
  return (float)(y00 * tx0 * tz0 + y01 * tx0 * tz1 + y10 * tx1 * tz0 + y11 * tx1 * tz1);
}
// calcule de la normale 3d au point32_t (i,j)
// avec pente calculee suivant les altitudes interpolees dans la direction de la plus grande pente
void calcule_norm3d_cel_pente_max_interp(int32_t i, int32_t j, float r) {
  static char start = 1;
  int32_t r0;
  Vec3 *n3d;
  int32_t x1 = 0, y1 = 0, z1 = 0,
          x2 = 0, y2 = 0, z2 = 0,
          x3 = 0, y3 = 0, z3 = 0,
          x4 = 0, y4 = 0, z4 = 0,
          ux = 0, uy = 0, uz = 0, vy = 0, vz = 0,
          wx = 0, wy = 0, wz = 0;
  float d, co, si, angle, r1, r2;
  float fx1, fx2, fy1, fy2, fz1, fz2, fux, fuy, fuz, fwx, fwy, fwz;
  if (!alti) {
    ErrPrintf("ERROR: calcule_norm3d_cel - alti = NULL\n");
    exit(-1);
  }
  if (!norm3d) {
    LogPrintf("calcul des normales 3d\n");
    AllocMemoryPrint("norm3d", norm3d, Vec3, LEO * LNS);
    ResetMemory(norm3d, Vec3, LEO * LNS);
  }
  /// calcul de la direction de plus grande pente
  r0 = 2;
  if (start) {
    LogPrintf("calcule_norm3d_cel_pente_max_interp : r = %f   r0 = %d\n", r, r0);
    start = 0;
  }
  n3d = norm3d + (i + j * LEO);
  /// U = tangente dans le plan XY
  x1 = i - r0;
  z1 = j;
  x2 = i + r0;
  z2 = j;
  apply_bc_cell(&x1, &z1);
  apply_bc_cell(&x2, &z2);
  ux = (pbc_mode) ? 2 * r0 : x2 - x1;
  y1 = Alti(x1, z1);
  y2 = Alti(x2, z2);
  uy = y2 - y1;
  /// V = tangente dans le plan YZ
  if (LNS >= 5) {
    x3 = i;
    z3 = j - r0;
    x4 = i;
    z4 = j + r0;
    apply_bc_cell(&x3, &z3);
    apply_bc_cell(&x4, &z4);
    y3 = Alti(x3, z3);
    y4 = Alti(x4, z4);
    vy = y3 - y4;
    vz = (pbc_mode) ? -2 * r0 : z3 - z4;
  } else {
    vz = 0;
  }
  /// normale W = U ^ V
  if (vz) {
    n3d->x = uy * vz / NB_MVT_EO;
    n3d->z = ux * vy;
  } else {
    n3d->x = -uy / NB_MVT_EO;
    n3d->z = 0.0;
  }
  /// direction de plus grande pente
  angle = (n3d->x) ? atan(n3d->z / n3d->x) : PI / 2;
  /// recalcul de la normale dans cette direction avec interpolation de l'altitude
  co = cos(angle);
  si = sin(angle);
  r1 = co * r;
  r2 = si * r;
  if (LNS < 2 * r2) { // tangente U
    r1 = r;
    r2 = 0;
  }
  /// coordonnees des points tests
  fx1 = i - r1;
  fz1 = j - r2;
  fx2 = i + r1;
  fz2 = j + r2;
  apply_bc_cel_interp(&fx1, &fz1);
  apply_bc_cel_interp(&fx2, &fz2);
  fux = (pbc_mode) ? 2 * r1 : fx2 - fx1;
  fuz = (pbc_mode) ? 2 * r2 : fz2 - fz1;
  fy1 = interpole_alti_cell(fx1, fz1);
  fy2 = interpole_alti_cell(fx2, fz2);
  fuy = fy2 - fy1;
  /// normale W
  fwx = -fuy * co;
  fwy = fuz * si + fux * co;
  fwz = -fuy * si;
  // normalisation
  d = sqrt(fwx * fwx + fwy * fwy + fwz * fwz);
  if (fwy < 0) {
    d = -d;
  }
  if (!d) {
    LogPrintf("ux = %d   uy = %d   uz = %d   wx = %d   wy = %d   wz = %d   d = %f\n", ux, uy, uz, wx, wy, wz, d);
  }
  n3d->x = fwx / d;
  n3d->y = fwy / d;
  n3d->z = fwz / d;
}

void calcule_normales() {
  int32_t x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4;
  int32_t ux, uy, vy, vz;
  int32_t wx, wy, wz;
  float d;
  int32_t i, j, ii;
  if (!norm2d) {
    LogPrintf("calcul des normales 2d et 3d\n");
    AllocMemory(norm2d, Vec2, L * D);
    ResetMemory(norm2d, Vec2, L * D);
    AllocMemory(norm3d, Vec3, L * D);
    ResetMemory(norm3d, Vec3, L * D);
  }
  if (!alti) {
    ErrPrintf("ERROR: calcule_normales - alti=NULL\n");
    exit(-1);
  }
  for (ii = 0, j = 0; j < LNS; j++)
    for (i = 0; i < LEO; i++, ii++) {
      if (rot_map && OutOfSpace(1 + i, LN + j)) {
        continue;
      }
      // calcul de la normale 2d (dans un plan vertical XY)
      x1 = i - 2;
      z1 = j;
      x2 = i + 2;
      z2 = j;
      apply_bc_cell(&x1, &z1);
      apply_bc_cell(&x2, &z2);
      y1 = Alti(x1, z1);
      y2 = Alti(x2, z2);
      ux = (pbc_mode) ? 4 : x2 - x1;
      uy = y2 - y1;
      Check_min_max(uy, -5, 5);
      d = sqrt(ux * ux + uy * uy);
      norm2d[ii].x = -uy / d;
      norm2d[ii].y = ux / d;
      // calcul de la normale 3d
      if (LNS >= 5) {
        x3 = i;
        z3 = j - 2;
        x4 = i;
        z4 = j + 2;
        apply_bc_cell(&x3, &z3);
        apply_bc_cell(&x4, &z4);
        vz = (pbc_mode) ? -4 : z3 - z4;
        y3 = Alti(x3, z3);
        y4 = Alti(x4, z4);
        vy = y3 - y4;
        Check_min_max(vy, -5, 5);
      } else {
        vz = 0;
      }
      // W = U ^ V
      if (vz) {
        wx = uy * vz;
        wy = - ux * vz;
        wz = ux * vy;
        d = sqrt(wx * wx + wy * wy + wz * wz);
        norm3d[ii].x = wx / d / NB_MVT_EO;
        norm3d[ii].y = wy / d;
        norm3d[ii].z = wz / d;
      } else {
        norm3d[ii].x = norm2d[ii].x;
        norm3d[ii].y = norm2d[ii].y;
        norm3d[ii].z = 0.0;
      }
//#endif
    }
}
#ifdef AVALANCHES
// callback de controle des transitions d'avalanches en fonction de la pente locale en surface
int32_t check_ava(int32_t ix, void *data) {
  static char start = 1;
  int32_t x, y, z, i, j, ij, i2, j2, al, al2, delta_x, delta_y, delta_z, chk_ok = 0;
  float ny;
  static float ny_lim = 0;
  float chk_angle = 0, r;
  Vec3 *pt_n3d = NULL;
  DataCheck *pdc = data;
  uint8_t typ;
  if (start) {
    LogPrintf("check_ava : ava_h_lim = %d\n", ava_h_lim);
    LogPrintf("check_ava : ava_angle = %f\n", ava_angle);
    if (!ava_upwind) {
      LogPrintf("check_ava: no avalanche upwind\n");
    }
    start = 0;
  }
  //calcul de la position (x,y,z)
  Calcule_xyz(ix, x, y, z);
  i = x - 1;
  j = z - LN;
  ij = i + j * LEO;
  //type de cellule
  typ = TE[ix].celltype;
  //altitude
  al = alti[ij];
  if (al + y < H) { //grain qui tombe
    return 0;
  }
  //direction d'avalanche
  if (pdc->dir == EST) { //est-ouest
    delta_x = (pdc->cel == 1) ? 1 : -1;
    delta_z = 0;
  } else { //nord-sud
    delta_x = 0;
    delta_z = (pdc->cel == 1) ? 1 : -1;
  }
  i2 = i + delta_x;
  j2 = j + delta_z;
#ifdef CYCLAGE_HOR
  apply_bc_cell(&i2, &j2);
#endif
  //calcul de la difference de hauteur
  al2 = Alti(i2, j2);
  delta_y = al - al2;
  if (delta_y < 0) {
    int32_t typ2 = TE[ix + delta_x + delta_z * HL].celltype;
    LogPrintf("check_ava : x = %04d   z = %04d   al = %04d   delta_x = %04d   delta_y = %04d   delta_z = %04d   typ = %d   typ2 = %d   dir = %d   cel = %d\n", x, z, al, delta_x, delta_y, delta_z, typ, typ2, pdc->dir, pdc->cel);
  }
  chk_angle = ava_angle;
  if (chk_angle) {
    r = (chk_angle <= 33.0) ? 3.5 : (chk_angle <= 45.0) ? 3.0 : 1.5;
    calcule_norm3d_cel_pente_max_interp(i, j, r);
    pt_n3d = norm3d + ij;
  }
  if (!ava_upwind && chk_angle && (pt_n3d->x <= 0)) {
    return 0;  //transition blocked when slope is windward
  }
  if (delta_y > ava_h_lim) {
    return 1;  //check the height difference
  }
  if (chk_angle) {
    //vertical component of the normal vector
    ny = pt_n3d->y;
    //max value
    ny_lim = 0.999 * cos(PI * ava_angle / 180.0); //coefficient 0.999 a cause des imprecisions de calcul
    chk_ok = (ny < ny_lim);
  }
  if (ava_mode == AVA_PROPAG) {
    // delai fixe
    if (chk_ok) {
      static float time_threshold = 0.0;
      //modele AVA
      float delai = ava_delay;
      if (!time_threshold) {
        time_threshold = delai;
      }
      if (csp_time > time_threshold) {
        //on declenche une grosse avalanche avec propagation
        LogPrintf("appel loop_ava_propag() : i=%d   j=%d   csp_time=%g\n", i, j, csp_time);
        LogPrintf("delai = %f\n", delai);
        loop_ava_propag(i, j);
        if (nb_ava_propag) {
          chk_ok = 0;
        }
        time_threshold = csp_time + delai;
      }
    }
  }
  return chk_ok;
}
#endif //AVALANCHES
// callback de controle des transitions de montee de grains en fonction de la pente locale en surface
// (on suppose que le lien de transitions se fera en est-ouest)
// TODO? : on pourrait aussi verifier la difference de hauteur en est-ouest avant de calculer la normale ...
int32_t check_slope(int32_t ix, void *data) {
  (void)data; //SUPPRESS unused warning
  static char start = 1;
  int32_t x, y, z, i, j, ij, chk_ok = 0;
  float nx;
  static float nx_lim = 0;
  Vec3 *pt_n3d;
  if (start) {
    LogPrintf("check_slope : slope_angle = %f\n", slope_angle);
    start = 0;
  }
  //calcul de la position (x,y,z)
  Calcule_xyz(ix, x, y, z);
  i = x - 1;
  j = z - LN;
  i += 2; // lire la normale en aval
  if (!pbc_mode) {
    if (i >= LEO) {
      i = LEO - 1;
    }
  }
#ifdef CYCLAGE_HOR
  else if (!rot_map) {
    if (i >= LEO) {
      i -= LEO;
    }
  } else {
    apply_bc_cell(&i, &j);
  }
#endif
  //altitude
  ij = i + j * LEO;
  //calcul de la normale 3d
  calcule_norm3d_cel_pente_max_interp(i, j, 3.0);
  pt_n3d = norm3d + ij;
  nx = pt_n3d->x;
  //critere
  if (!nx_lim) {
    nx_lim = -sin(slope_angle * PI / 180.0);
  }
  chk_ok = (nx > nx_lim);
  return chk_ok;
}
#ifdef CGV
float dist_grdv = 10.0; //10.0  //distance caracteristique pour le calcul du gradient de la vitesse
float dist0_grdv = 0.0; //0.0  //distance de reference pour le calcul du gradient de la vitesse
float grdvc = 0.0; //constante critique pour le controle de l'erosion en fonction du gradient de vitesse
float grdvc_max = 0.0; //seuil max pour le controle de l'erosion en fonction du gradient de vitesse
float grdvc_min = 0.0; //seuil min pour le controle de l'erosion en fonction du gradient de vitesse
float cgv_coef = 0.0; //valeur moyenne de probabilite de blocage d'une transition
extern float *Velx_interp, *Vely_interp;
extern float maxvel, meanvel;
void reset_grad_vel() {
  if (grdv) {
    ResetMemory(grdv, float, L * D);
  }
}
void calcule_grad_vel() {
  int32_t i, j, ii, x0, y0, z0, x2, y2, z2, ix2;
  float vx2, vy2, vel2, tgx, tgy;
  if (!grdv) {
    AllocMemory(grdv, float, L * D);
    ResetMemory(grdv, float, L * D);
    LogPrintf("calcule_grad_vel : dist0_grdv = %f   dist_grdv = %f\n", dist0_grdv, dist_grdv);
  }
  if (!alti || !norm2d || !Velx_interp) {
    WarnPrintf("WARNING: cannot compute grdv (%lu,%lu,%lu)\n", (long)alti, (long)norm2d, (long)Velx_interp);
    return;
  }
  for (ii = j = 0; j < LNS; j++)
    for (i = 0; i < LEO; i++, ii++) {
      if (rot_map && OutOfSpace(1 + i, LN + j)) {
        continue;
      }
      x0 = i + 1;
      y0 = H - alti[ii];
      z0 = j + LN;
      //vecteur tangent a la surface en 2d
      tgx = norm2d[ii].y;
      tgy = -norm2d[ii].x;
      x2 = x0 + roundf(norm3d[ii].x * dist_grdv);
      y2 = y0 - roundf(norm3d[ii].y * dist_grdv);
      z2 = z0 + roundf(norm3d[ii].z * dist_grdv);
      if (x2 < 0) {
        x2 = 0;
      }
      if (x2 >= L) {
        x2 = L - 1;
      }
      if (y2 < 0) {
        y2 = 0;
      }
      //plafond
      if (y2 >= H) {
        y2 = H - 1;
      }
      if (z2 < 1) {
        z2 = 1;
      }
      if (z2 >= D - 1) {
        z2 = D - 2;
      }
      ix2 = z2 * HL + y2 * L + x2;
      if (rot_map && pbc_mode && OutOfSpace(x2, z2)) {
        ix2 = hcycle(ix2);
      }
      if ((ix2 < 0) || (ix2 >= HLD)) {
        ErrPrintf("WARNING: calcule_grad_vel, i=%d, j=%d, ix2=%d\n", i, j, ix2);
        continue;
      }
      vx2 = Velx_interp[ix2];
      vy2 = Vely_interp[ix2];
      //norme du vecteur vitesse
      vel2 = vx2 * tgx + vy2 * tgy; //projection du vecteur vitesse sur la tangente a la surface
      grdv[ii] = vel2;
      //gradient max quand la normale traverse le plafond
    }


}
//fonction de probabilite pour le controle en fonction du gradient de la vitesse
float prob_cgv(float gv, float gvmin, float gvmax) {
  static char start = 1;
  if (start) {
    if (!grdvc) {
      grdvc = maxvel / 4;  //grdv_max/4
    }
    //grdv_max/4
    if (!grdvc_min) {
      grdvc_min = 0.0;
    }
    if (!grdvc_max) {
      grdvc_max = maxvel;
    }
    LogPrintf("grdvc_min = %f\n", grdvc_min);
    LogPrintf("grdvc = %f\n", grdvc);
    LogPrintf("grdvc_max = %f\n", grdvc_max);
    gvmin = grdvc_min;
    gvmax = grdvc_max;
    start = 0;
  }



  float prob = ((gv < gvmin) || (gv < 0)) ? 0 : (gv > gvmax) ? 1 : (gv - gvmin) / (gvmax - gvmin);
  return prob;
}
// callback de controle en fonction du gradient de la vitesse selon la normale a la surface
int32_t check_grad_vel(int32_t ix, void *data) {
  (void)data; //SUPPRESS unused warning
  int32_t x, y, z, chk_ok;
  float gv, alea, prob;
  //calcul de la position (x,y,z)
  Calcule_xyz(ix, x, y, z);
  gv = grdv[(x - 1) + (z - LN) * LEO];
  alea = drand48();
  prob = prob_cgv(gv, grdvc_min, grdvc_max);
  chk_ok = (prob > alea);
  //test
  return chk_ok;
}
void compute_coef_cgv() {
  int32_t i, k, ik;
  int32_t n0, nn, np;
  float sum_gv, gv_mean;
  static char start = 1;
  if (!grdv) {
    return;
  }
  if (start) {
    LogPrintf("calcul coef_cgv pour correction du pas de temps\n");
    start = 0;
  }
  //estimation du coefficient correctif sur le pas de temps (dans le cas d'une loi de probabilite lineaire a seuils)
  sum_gv = 0.0;
  n0 = np = nn = 0;
  for (k = 0; k < LNS; k++) {
    if (!psol[LN + k]) {
      ik = k * LEO;
      for (i = 0; i < LEO; i++, ik++) {
        if (rot_map && OutOfSpace(1 + i, LN + k)) {
          continue;
        }
        if (grdv[ik] > grdvc_max) {
          np++;
        } else if (grdv[ik] > grdvc_min) {
          n0++;
          sum_gv += grdv[ik];
        }
        nn++;
      }
    }
  }
  gv_mean = n0 ? (sum_gv / n0) : 0.0;
  cgv_coef = nb_pv ? (prob_cgv(gv_mean, grdvc_min, grdvc_max) * n0 + np) / nn : 0.0;
  assert(cgv_coef <= 1.0);
}
#endif //CGV
#ifdef AVALANCHES
//avalanches dans 8 directions equiprobables
void avalanches(uint8_t typ, int16_t h_lim, int16_t nb_cel_max, char alti_mode) {
  static char first = 1;
  Cell *aux, cel;
  short *pt, *pt2;
  int i, j, k, jmin, jmax, ix, ix2;
  unsigned char alea, cpt, flag_ava;
  int step[8], step_alti[8];
  short i_cel, h_lim_diag;
  int nb_ava, nb_ava_iter;
  char output[128]; // things to be written to AVA log

  h_lim_diag = h_lim;

  if (first) {
    LogPrintf("avalanches : h_lim = %d\n", h_lim);
    LogPrintf("avalanches : h_lim_diag = %d\n", h_lim_diag);
    if (pbc_mode && rot_map) {
      ErrPrintf("ERROR: sorry, rotating space is not compatible with synchronous avalanches\n");
      exit(-1);
    }
    first = 0;
  }

  //avalanches
  nb_ava = nb_ava_iter = 0;
  jmin = 0;
  jmax = LNS - 1;
  do {
    nb_ava_iter++;
    flag_ava = 0;
    pt = alti + jmin * LEO;
    for (j = jmin; j <= jmax; j++) {
      for (i = 0; i < LEO; i++, pt++) {
        // on compare l'elevation locale dans les 8 directions
        cpt = 0;
        if ((i > 0) && (*pt - * (pt - 1) > h_lim)) {
          step_alti[cpt] = -1;    //gauche
          step[cpt++] = -1;
        }
        if ((i < LEO - 1) && (*pt - * (pt + 1) > h_lim)) {
          step_alti[cpt] = 1;    //droite
          step[cpt++] = 1;
        }
        if ((j > jmin) && (*pt - * (pt - LEO) > h_lim)) {
          step_alti[cpt] = -LEO;    //derriere
          step[cpt++] = -HL;
        }
        if ((j < jmax) && (*pt - * (pt + LEO) > h_lim)) {
          step_alti[cpt] = LEO;    //devant
          step[cpt++] = HL;
        }
        if ((i > 0) && (j > jmin) && (*pt - * (pt - LEO - 1) > h_lim_diag)) {
          step_alti[cpt] = -LEO - 1;  //diagonale gauche+derriere
          step[cpt++] = -HL - 1;
        }
        if ((i < LEO - 1) && (j > jmin) && (*pt - * (pt - LEO + 1) > h_lim_diag)) {
          step_alti[cpt] = -LEO + 1;  //diagonale droite+derriere
          step[cpt++] = -HL + 1;
        }
        if ((i > 0) && (j < jmax) && (*pt - * (pt + LEO - 1) > h_lim_diag)) {
          step_alti[cpt] = LEO - 1;  //diagonale gauche+devant
          step[cpt++] = HL - 1;
        }
        if ((i < LEO - 1) && (j < jmax) && (*pt - * (pt + LEO + 1) > h_lim_diag)) {
          step_alti[cpt] = LEO + 1;  //diagonale droite+devant
          step[cpt++] = HL + 1;
        }
#ifdef CYCLAGE_HOR
        if (pbc_mode) {
          if ((i == 0) && (*pt - * (pt + LEO - 1) > h_lim)) {
            step_alti[cpt] = LEO - 1;  //gauche
            step[cpt++] = L - 3;
          }
          if ((i == LEO - 1) && (*pt - * (pt - LEO + 1) > h_lim)) {
            step_alti[cpt] = -(LEO - 1);  //droite
            step[cpt++] = -(L - 3);
          }
          if ((LN + j == 0) && (*pt - * (pt + (LNS - 1)*LEO) > h_lim)) {
            step_alti[cpt] = (LNS - 1) * LEO; //derriere
            step[cpt++] = (D - 3) * HL;
          }
          if ((LN + j == D - 2) && (*pt - * (pt - (LNS - 1)*LEO) > h_lim)) {
            step_alti[cpt] = -(LNS - 1) * LEO; //devant
            step[cpt++] = -(D - 3) * HL;
          }
        }
#endif
        if (cpt) {
          // choix d'une direction pour l'avalanche
          for (i_cel = 0; i_cel < nb_cel_max; i_cel++) {
            alea = floor(cpt * drand48());
            pt2 = pt + step_alti[alea];
            if ((pt2 < alti) || (pt2 > alti + LEO * LNS)) {
              ErrPrintf("ERROR: avalanche, array overflow in alti\n");    //antibug
              exit(-1);
            }
            k = H - (*pt);
            ix = 1 + i + k * L + (LN + j) * HL; //cellule du sommet
            aux = TE + ix;
            cel = *aux;
            if (aux->celltype == DUM) {
              break;
            }
            if (aux->celltype != typ) {
              ErrPrintf("WARNING: avalanche, incorrect type : %d (i = %d  j = %d  k = %d  deniv = %d  step_alti = %d  nb_ava_iter = %d)\n", aux->celltype, i, j, k, *pt - *pt2, step_alti[alea], nb_ava_iter);
              break;
            }
            if ((*pt) > (*pt2)) {
              ix2 = ix + step[alea];  //cellule voisine
              deplace_cellule(ix2, ix);
              //on fait remonter la colonne de cellules
              for (k = *pt; k > (*pt2) + 1; k--) {
                ix = ix2;
                ix2 += L;
                deplace_cellule(ix2, ix);
              }
            } else {
              continue;
            }
            // on fait tomber la cellule haute dans la direction choisie
            init_cellule(cel, ix2);
            //mise-a-jour des elevations locales
            (*pt)--;
            (*pt2)++;
            if (alti_mode == ALTI_MODE_HAUT) { //recalcul de l'altitude
              aux += L;
              while (Phase[aux->celltype] != SOLID) {
                (*pt)--;
                aux += L;
              }
            }
            flag_ava = 1;
            nb_ava++;
          }
        }
      }
    }
  } while (flag_ava);

  //dump avalanches
  sprintf(output, "avalanches : %d   passes : %d   csp_time : %f   h_lim : %d\n", nb_ava, nb_ava_iter, csp_time, h_lim); //fflush(stdout);
  output_write("AVA", output);

  if (nb_ava) {
    //recalcul des tableaux de doublets actifs
    init_db_pos();
    //recalcul des cellules
    init_Ncel();
  }

}
//avalanches dans 8 directions selon les normales
//void avalanches_norm(uint8_t typ, float angle, int16_t nb_cel_max, char alti_mode)
void avalanches_norm(uint8_t typ, int16_t nb_cel_max, char alti_mode) {
  static char first = 1;
  Cell *aux, cel;
  int16_t *pt, *pt2;
  float nx, ny, nz, ny_lim, anh, r;
  int32_t i, j, k, jmin, jmax, ix, ix2, ij;
  uint8_t flag_drop;
  int32_t step, step_alti;
  int16_t i_cel;
  int32_t nb_ava, nb_ava_iter, nb_ava_drop;
  Vec3 *pt_n3d;
  float angle;
  char output[256];

  angle = ava_angle;
  ny_lim = cos(angle * PI / 180);
  r = (angle <= 33.0) ? 3.5 : 3.0;
  if (first) {
    LogPrintf("avalanches_norm : ava_angle = %f    r = %f\n", angle, r);
    if (pbc_mode && rot_map) {
      ErrPrintf("ERROR: sorry, rotating space is not compatible with synchronous avalanches\n");
      exit(-1);
    }
    first = 0;
  }
  //avalanches
  nb_ava = nb_ava_iter = nb_ava_drop = 0;
  jmin = 0;
  jmax = LNS - 1;
  do {
    nb_ava_iter++;
    ij = 0;
    for (j = jmin; j <= jmax; j++) {
      for (i = 0; i < LEO; i++, ij++) {
        calcule_norm3d_cel_pente_max_interp(i, j, r);
        pt_n3d = norm3d + ij;
        ny = pt_n3d->y;
        //if (0){
        if (ny < ny_lim) {
          pt = alti + ij;
          nx = pt_n3d->x; //composante EO
          nz = pt_n3d->z;  //composante NS
          //calcul de l'angle au sol, pour determiner la direction de la plus forte pente
          if (nx) {
            anh = atan(-nz / nx);
            if (nx < 0) {
              anh += PI;
            }
            if (anh < 0) {
              anh += 2 * PI;
            }
          } else {
            anh = (nz > 0) ? 3 * PI / 2 : PI / 2;
          }
          step_alti = 0;
          if ((anh < PI / 8) || (anh > 15 * PI / 8)) {
            if (i < LEO - 1) {
              step_alti = 1;  //droite
              step = 1;
            }
          } else if (anh < 3 * PI / 8) {
            if ((i < LEO - 1) && (j > jmin)) {
              step_alti = -LEO + 1;  //diagonale droite+derriere
              step = -HL + 1;
            }
          } else if (anh < 5 * PI / 8) {
            if (j > jmin) {
              step_alti = -LEO;  //derriere
              step = -HL;
            }
          } else if (anh < 7 * PI / 8) {
            if ((i > 0) && (j > jmin)) {
              step_alti = -LEO - 1;  //diagonale gauche+derriere
              step = -HL - 1;
            }
          } else if (anh < 9 * PI / 8) {
            if (i > 0) {
              step_alti = -1;  //gauche
              step = -1;
            }
          } else if (anh < 11 * PI / 8) {
            if ((i > 0) && (j < jmax)) {
              step_alti = LEO - 1;  //diagonale gauche+devant
              step = HL - 1;
            }
          } else if (anh < 13 * PI / 8) {
            if (j < jmax) {
              step_alti = LEO;  //devant
              step = HL;
            }
          } else if (anh < 15 * PI / 8) {
            if ((i < LEO - 1) && (j < jmax)) {
              step_alti = LEO + 1;  //diagonale droite+devant
              step = HL + 1;
            }
          }
          if (step_alti) {
            for (i_cel = 0; i_cel < nb_cel_max; i_cel++) {
              pt2 = pt + step_alti;
              if ((pt2 < alti) || (pt2 > alti + LEO * LNS)) {
                ErrPrintf("WARNING: avalanches_norm - array overflow in alti\n");  //antibug
                break;
              }
              k = H - (*pt);
              ix = 1 + i + k * L + (LN + j) * HL; //cellule du sommet
              aux = TE + ix;
              cel = *aux;
              if (aux->celltype == DUM) {
                break;
              }
              if (aux->celltype != typ) {
                ErrPrintf("WARNING: avalanches_norm - type incorrect : %d (i = %d  j = %d  k = %d  nb_ava_iter = %d)\n", aux->celltype, i, j, k, nb_ava_iter);
                break;
              }
              if ((*pt) - (*pt2) > 0) {
                if ((*pt) - (*pt2) > 1) {
                  flag_drop = 1; //indicateur de chute verticale, pour forcer l'arret des avalanches en l'absence de chute
                  nb_ava_drop++;
                }
                ix2 = ix + step;  //cellule voisine
                deplace_cellule(ix2, ix);
                //on fait remonter la colonne de cellules
                for (k = *pt; k > (*pt2) + 1; k--) {
                  ix = ix2;
                  ix2 += L;
                  deplace_cellule(ix2, ix);
                }
              } else {
                continue;
              }
              // on fait tomber la cellule haute dans la direction choisie
              init_cellule(cel, ix2);
              //mise-a-jour des elevations locales
              (*pt)--;
              (*pt2)++;
              if (alti_mode == ALTI_MODE_HAUT) { //recalcul de l'altitude
                aux += L;
                while (Phase[aux->celltype] != SOLID) {
                  (*pt)--;
                  aux += L;
                }
              }
              nb_ava++;
            }
          }
        }
      }
    }
  } while (flag_drop);
  //dump avalanches
  sprintf(output,"avalanches_norm : %d (mov) %d (drop)   passes : %d   csp_time : %f   angle : %f\n", nb_ava, nb_ava_drop, nb_ava_iter, csp_time, angle); //fflush(stdout);
  output_write("AVA", output);

  if (nb_ava){
    //recalcul des tableaux de doublets actifs
    init_db_pos();
    //recalcul des cellules
    init_Ncel();
  }
}
// avalanches localisees avec propagation vers les 4 premiers voisins (fonction recursive)
// TODO : mode CELL_COLOR
// TODO : mode CYCLAGE_HOR
void ava_propag(int32_t i, int32_t j) {
  int16_t *pt, *pt2;
  uint8_t flag_ava, flag_drop;
  int32_t step, step_alti;
  int32_t jmin, jmax, ix, ix2, ij;
  int16_t h, k;
  if (pbc_mode && rot_map) {
    ErrPrintf("ERROR: sorry, rotating space is not compatible with synchronous avalanches\n");
    exit(-1);
  }
  flag_ava = flag_drop = 0;
  jmin = 0;
  jmax = LNS - 1;
  step_alti = 0;
  //Verification des differences d'altitude avec les 4 premiers voisins
  ij = i + j * LEO;
  ava_mask[ij] = 0;
  pt = alti + ij;
  if (i + 1 < LEO) {
    h = (*pt) - (*(pt + 1));
    if (h > ava_h_lim) {
      step_alti = 1;  //droite
      step = 1;
    }
  }
  if (!step_alti && (i > 0)) {
    h = (*pt) - (*(pt - 1));
    if (h > ava_h_lim) {
      step_alti = -1;  //gauche
      step = -1;
    }
  }
  if (!step_alti && (j + 1 < LNS)) {
    h = (*pt) - (*(pt + LEO));
    if (h > ava_h_lim) {
      step_alti = LEO;  //devant
      step = HL;
    }
  }
  if (!step_alti && (j > 0)) {
    h = (*pt) - (*(pt - LEO));
    if (h > ava_h_lim) {
      step_alti = -LEO;  //derriere
      step = -HL;
    }
  }
#ifdef CYCLAGE_HOR
  if (pbc_mode) {
    if (!step_alti && (i == 0)) {
      h = (*pt) - (*(pt + LEO - 1));
      if (h > ava_h_lim) {
        step_alti = LEO - 1;  //gauche
        step = L - 3;
      }
    }
    if (!step_alti && (i == LEO - 1)) {
      h = (*pt) - (*(pt - LEO + 1));
      if (h > ava_h_lim) {
        step_alti = -(LEO - 1);  //droite
        step = -(L - 3);
      }
    }
    if (!step_alti && (j == 0)) {
      h = (*pt) - (*(pt + (LNS - 1) * LEO));
      if (h > ava_h_lim) {
        step_alti = (LNS - 1) * LEO;  //derriere
        step = (L - 3) * HL;
      }
    }
    if (!step_alti && (j == LNS - 1)) {
      h = (*pt) - (*(pt - (LNS - 1) * LEO));
      if (h > ava_h_lim) {
        step_alti = -(LNS - 1) * LEO;  //devant
        step = -(L - 3) * HL;
      }
    }
  }
#endif
  if (!step_alti) {
    //Chute dans la direction de plus grande pente (4 ou 8 directions possibles)
    float r, nx, ny, nz, ny_lim, anh;
    Vec3 *pt_n3d;
    ny_lim = cos(PI * ava_propag_angle / 180);
    r = (ava_propag_angle <= 33.0) ? 3.5 : 3.0;
    calcule_norm3d_cel_pente_max_interp(i, j, r);
    pt_n3d = norm3d + ij;
    ny = pt_n3d->y;
    if (ny < ny_lim) {
      nx = pt_n3d->x; //composante EO
      nz = pt_n3d->z;  //composante NS
      //calcul de l'angle au sol, pour determiner la direction de la plus forte pente
      if (nx) {
        anh = atan(-nz / nx);
        if (nx < 0) {
          anh += PI;
        }
        if (anh < 0) {
          anh += 2 * PI;
        }
      } else {
        anh = (nz > 0) ? 3 * PI / 2 : PI / 2;
      }
      if ((anh < PI / 8) || (anh > 15 * PI / 8)) {
        if (i < LEO - 1) {
          step_alti = 1;  //droite
          step = 1;
        }
      } else if (anh < 3 * PI / 8) {
        if ((i < LEO - 1) && (j > jmin)) {
          step_alti = -LEO + 1;  //diagonale droite+derriere
          step = -HL + 1;
        }
      } else if (anh < 5 * PI / 8) {
        if (j > jmin) {
          step_alti = -LEO;  //derriere
          step = -HL;
        }
      } else if (anh < 7 * PI / 8) {
        if ((i > 0) && (j > jmin)) {
          step_alti = -LEO - 1;  //diagonale gauche+derriere
          step = -HL - 1;
        }
      } else if (anh < 9 * PI / 8) {
        if (i > 0) {
          step_alti = -1;  //gauche
          step = -1;
        }
      } else if (anh < 11 * PI / 8) {
        if ((i > 0) && (j < jmax)) {
          step_alti = LEO - 1;  //diagonale gauche+devant
          step = HL - 1;
        }
      } else if (anh < 13 * PI / 8) {
        if (j < jmax) {
          step_alti = LEO;  //devant
          step = HL;
        }
      } else if (anh < 15 * PI / 8) {
        if ((i < LEO - 1) && (j < jmax)) {
          step_alti = LEO + 1;  //diagonale droite+devant
          step = HL + 1;
        }
      }
    }
  }
#ifdef CYCLAGE_HOR
  //TODO
#endif
  if (step_alti) { //la cellule se deplace dans la direction step_alti
    Cell *aux, cel;
    pt2 = pt + step_alti;
    if ((pt2 < alti) || (pt2 > alti + LEO * LNS)) {
      ErrPrintf("WARNING: ava_propag - array overflow in alti\n");  //antibug
      return;
    }
    k = H - (*pt);
    ix = 1 + i + k * L + (LN + j) * HL; //cellule du sommet
    aux = TE + ix;
    cel = *aux;
    if (aux->celltype == DUM) {
      return;
    }
    if (aux->celltype != ALTI) {
      ErrPrintf("WARNING: ava_propag - incorrect type %d (i = %d  j = %d  k = %d )\n", aux->celltype, i, j, k);
      return;
    }
    if ((*pt) - (*pt2) > 0) {
      if ((*pt) - (*pt2) > 1) {
        flag_drop = 1; //indicateur de chute verticale, pour forcer l'arret des avalanches en l'absence de chute
        nb_drop_ava_propag++;
      }
      ix2 = ix + step;  //cellule voisine
      deplace_cellule(ix2, ix);
      for (k = *pt; k > (*pt2) + 1; k--) {
        ix = ix2;
        ix2 += L;
        deplace_cellule(ix2, ix);
      }
    } else {
      return;
    }
    init_cellule(cel, ix2);
    //mise-a-jour des elevations locales
    (*pt)--;
    (*pt2)++;
    flag_ava = 1;
    nb_ava_propag++;
  }
  if (flag_ava) { //Propagation des avalanches
    Pos2 p2;
    if ((i + 1 < LEO) && !ava_mask[ij + 1]) { //propagation a droite
      p2.x = i + 1;
      p2.y = j;
      fifo_put(ava_fifo, p2);
      ava_mask[ij + 1] = 1;
    }
    if ((i > 0) && !ava_mask[ij - 1]) { //propagation a gauche
      p2.x = i - 1;
      p2.y = j;
      fifo_put(ava_fifo, p2);
      ava_mask[ij - 1] = 1;
    }
    if ((j + 1 < LNS) && !ava_mask[ij + LEO]) { //propagation devant
      p2.x = i;
      p2.y = j + 1;
      fifo_put(ava_fifo, p2);
      ava_mask[ij + LEO] = 1;
    }
    if ((j > 0) && !ava_mask[ij - LEO]) { //propagation derriere
      p2.x = i;
      p2.y = j - 1;
      fifo_put(ava_fifo, p2);
      ava_mask[ij - LEO] = 1;
    }
  }
}
void loop_ava_propag(int32_t i, int32_t j) {
  Pos2 p2;
  int32_t nb_iter, nb_fifo_max;
  if (!ava_mask) {
    AllocMemoryPrint("ava_mask", ava_mask, char, LEO * LNS);
    ResetMemory(ava_mask, char, LEO * LNS);
  }
  if (!ava_fifo) {
    ava_fifo = fifo_create(LEO * LNS);
  }
  nb_ava_propag = nb_drop_ava_propag = 0;
  nb_iter = nb_fifo_max = 0;
  ava_propag_angle = ava_angle_stable;
  LogPrintf("ava_propag_angle = %f\n", ava_propag_angle);
  ava_propag(i, j);
  nb_iter++;
  while (ava_fifo->nb) {
    if (ava_fifo->nb > nb_fifo_max) {
      nb_fifo_max = ava_fifo->nb;
    }
    p2 = fifo_get(ava_fifo);
    ava_propag(p2.x, p2.y);
    nb_iter++;
  }
  LogPrintf("nb_ava_propag = %d\n", nb_ava_propag);
  LogPrintf("nb_drop_ava_propag = %d\n", nb_drop_ava_propag);
  LogPrintf("nb_iter = %d\n", nb_iter);
  LogPrintf("nb_fifo_max = %d\n", nb_fifo_max);
  if (nb_ava_propag) {
    //recalcul des tableaux de doublets actifs
    init_db_pos();
    //recalcul des cellules
    init_Ncel();
  }
}
#endif //AVALANCHES
////////////////////////////////////////////////////////////////////////////////////:

void dump_surface(char* name, int32_t cpt, int32_t unit)
{
  char filename[100];
  int32_t i,j,n;
  char current_output[128];
  if (unit == UNIT_COMP)
    sprintf(filename,"%s%04d", name, cpt);
  else
    sprintf(filename,"%s%05d_t0", name, cpt);

  if (alti){
    int16_t *pt_al = alti;
    for(j=0; j<LNS; j++){
      n=0;
      for(i=0; i<LEO; i++, pt_al++){
        if (rot_map && OutOfSpace(1+i,LN+j)) continue;
        sprintf(current_output, "%d ", *pt_al);
        output_write(filename, current_output);
        n++;
      }
      if (n>0) output_write(filename, "\n");
    }

  }
}


#ifdef CGV
void dump_grad_vel(int32_t cpt, int32_t unit) {
  static char filename[100];
  FILE *fp;
  float *pt_gv;
  int32_t i, j, n;
  if (grdv) {
    if (unit == UNIT_COMP) {
      sprintf(filename, "GRAD_VEL%04d.data", cpt);
    } else {
      sprintf(filename, "GRAD_VEL%05d_t0.data", cpt);
    }
    fp = fopen(filename, "w");
    if (!fp) {
      ErrPrintf("erreur ouverture fichier vitesse\n");
      exit(-4);
    }
    pt_gv = grdv;
    for (j = 0; j < LNS; j++) {
      n = 0;
      for (i = 0; i < LEO; i++, pt_gv++) {
        fprintf(fp, "%f ", *pt_gv);
        n++;
      }
      if (n > 0) {
        fprintf(fp, "\n");
      }
    }
    fclose(fp);
  }
}

void dump_cgv()
{
  char current_output[128];
  static float grdv_max = VSTEP_H*VSTEP_L*VSTEP_TIME*2;
  float gv, prob;

  for(gv=Min(0,grdvc_min); gv<grdv_max; gv++){
    prob=prob_cgv(gv, grdvc_min, grdvc_max);
    sprintf(current_output, "%f   %f\n", gv, prob);
    output_write("PROB_CGV", current_output);
  }
}


void dump_cgv_coef()
{
  static int32_t cpt=0;
  char output_string[256];

  sprintf(output_string, "\n%04d:   %.4f", cpt++, cgv_coef);
  output_write("CGV_COEF", output_string);
}

#endif //CGV
#endif //ALTI
