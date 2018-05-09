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
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */



#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <math.h>
#include <assert.h>

#include "defs.h"
#include "macros.h"
#include "param.h"
#include "space.h"
#include "doublets.h"
#include "cells.h"
#include "surface.h"
#include "lgca.h"
#ifdef PARALLEL
#include "synchro.h"
#endif
#include "transitions.h"
#include "simul.h"

extern double csp_time;                  // temps reel simule
extern int H, L, D, HL, HLD;       // les dimensions de la terre
extern int LN, LS, LEO, LNS, HLN;    //couloir est-ouest (limite nord, limite sud, largeur nord-sud, ...)
extern Cell  *TE;	           // la 'terre'
extern int pbc_mode;  //periodic boundary conditions
extern int ava_mode;  //mode d'avalanches
extern int ava_h_lim;
extern float ava_delay; //delai entre avalanches
extern float ava_angle;
extern float ava_angle_stable;
extern float ava_angle_col;
extern int ava_upwind;
#ifdef PARALLEL
extern int mode_par;    //mode parallele
#endif
extern int nb_pv;                  // nombre de plans verticaux (non-DUM)
extern char *psol;                    // indicateur des plans solides verticaux est-ouest
extern char *rot_map;        // periodic mapping of the rotating space
extern Pos2 *rot_map_pos;        // periodic mapping of the rotating space

#ifdef PHASES
const unsigned char Phase[MAX_CELL] = PHASES;   //phase (fluide ou solide) des types de cellules
#endif // PHASES

#ifdef ALTI
short *alti = NULL;           // elevations locales de terrain
//float *alti_mean = NULL;      // elevations locales de terrain, après moyennage
Vec2 *norm2d = NULL;          // normale 2d a la surface definie par alti
Vec3 *norm3d = NULL;          // normale 3d a la surface definie par alti
Vec3 *norm3d_iso = NULL;      // normale 3d a la surface definie par alti (avec tangentes isotropes)
float angle_norm3d_iso=0.0;   // direction variable pour le calcul des normales 3d (avec tangentes isotropes)
float *grdv = NULL;           // gradient de vitesse en surface
float slope_angle=20.0;       // pente maximale pour la montee des grains
float ava_propag_angle=30.0;  // pente minimale pour les avalanches avec propagation
int nb_ava_propag=0;          // nombre de cellules deplacees par propagation d'avalanches
int nb_drop_ava_propag=0;     // nombre de cellules qui tombent par propagation d'avalanches
char *ava_mask = NULL;        // tableau de localisation des avalanches
FifoPos2* ava_fifo = NULL;     // file ordonnee (circulaire) des avalanches
int veg_h_max=0;         // hauteur max pour la vegetation

void loop_ava_propag(int i, int j);

void params_surface()
{
#ifdef USE_VEGETATION
  parameter("Veg_h_max", "maximal height for vegetation", &veg_h_max, PARAM_INT, "VEGETATION");
#endif
}

// creation d'une file FIFO
FifoPos2* fifo_create(int len)
{
  FifoPos2* fifo;

  AllocMemory(fifo, FifoPos2, 1);
  ResetMemory(fifo, FifoPos2, 1);

  fifo->length = len;

  AllocMemoryPrint("fifo->array", fifo->array, Pos2, len);
  ResetMemory(fifo->array, Pos2, len);

  return fifo;
}

// ajout d'un element dans la file
void fifo_put(FifoPos2 *fifo, Pos2 elt)
{
  assert(fifo->nb < fifo->length);
  fifo->array[fifo->tail++] = elt;
  if (fifo->tail >= fifo->length) fifo->tail = 0;
  fifo->nb++;
}

// enlevement d'un element dans la file
Pos2 fifo_get(FifoPos2 *fifo)
{
  Pos2 elt;
  assert(fifo->nb > 0);
  elt = fifo->array[fifo->head++];
  if (fifo->head >= fifo->length) fifo->head = 0;
  fifo->nb--;
  return elt;
}


// modification de l'altitude par ajout d'une cellule de type 'typ', en ix
// (code optimise)
void modif_alti_cel(int ix, unsigned char typ)
{
  int x, y, z, i, j;
  int al;

  //LogPrintf("modif_alti_cel - debut : ix = %d\n", ix);

  Calcule_xyz(ix, x, y, z);

  i = x-1;
  j = z-LN;

  al = Alti(i, j);

  if (Phase[typ] == SOLID){
    if (H-y == al+1){
      Cell *aux;
      aux = TE+ix-L;
      al++;
      while ((al < H-2) && (Phase[aux->celltype] == SOLID)) {al++; aux -= L;} //on remonte a la surface
      Alti(i, j) = al;
    }
  }
  else{
    if (H-y <= al) Alti(i, j) = H-y-1; //disparition d'une cellule dans la pile
  }

  //LogPrintf("modif_alti_cel - fin\n");
}


void calcule_alti(unsigned char typ, char alti_mode)
{
  Cell *aux;
  short *pt;
  int i, j, k;

  //LogPrintf("calcule_alti : typ = %d\n", typ);

  if (!alti){
    AllocMemoryPrint("alti", alti, short, L*D/*LEO*LNS*/);
    //alti = (short *)malloc((L-2)*LNS*sizeof(short));
    //if (!alti) {ErrPrintf("ERREUR : calcule_alti - allocation impossible\n");  exit(-1);}
  }

  //calcul des elevations locales
  pt = alti;
  if (alti_mode == ALTI_MODE_BAS){  //mode bas : on part d'en bas et on remonte jusqu'a trouver une cellule non-solide
    for(j=0; j < LNS; j++)
      for(i=0; i < LEO; i++, pt++){
        k = H-2;
        aux = TE+(1+i)+k*L+(LN+j)*HL;
        //while ((k > 0) && ((aux->celltype == typ) || (aux->celltype == DUM) || (aux->celltype == OUT) || (aux->celltype == IN))) {k--; aux -= L;}
        while ((k > 0) && (Phase[aux->celltype] == SOLID)) {k--; aux -= L;}
        *pt = H-(k+1);
      }
  }
  else{ //mode haut : on part d'en haut et on descend jusqu'a trouver une cellule solide
    for(j=0; j < LNS; j++)
      for(i=0; i < LEO; i++, pt++){
        k = 1;
        aux = TE+(1+i)+k*L+(LN+j)*HL;
        //while ((k < H) && (aux->celltype != typ) && (aux->celltype != DUM) && (aux->celltype != OUT) && (aux->celltype != BORD)) {k++; aux += L;}
        while ((k < H) && (Phase[aux->celltype] != SOLID)) {k++; aux += L;}
        *pt = H-k;
      }
  }
}


int calcule_alti_max(unsigned char typ)
{
  Cell *aux;
  short max=0;
  short *pt;
  int i, j, k;

  for(j=0; j < LNS; j++)
    for(i=0; i < LEO; i++, pt++){
      k = 1;
      aux = TE+(1+i)+k*L+(LN+j)*HL;
      while ((k < H-1) && (aux->celltype != typ)) {k++; aux += L;}
      if (max < H-k) max = H-k;
    }

  //LogPrintf("calcule_alti_max: max=%d\n",max);
  return max;
}

#ifndef LGCA
#define NB_MVT_EO 1
#define NB_MVT_VER 1
#endif

// apply boundary conditions on one cell
void apply_bc_cell(int *pt_x, int *pt_y)
{
  if (!pbc_mode){
    /// open or closed boundaries
    if (*pt_x<0) *pt_x = 0;
    if (*pt_x>=LEO) *pt_x = LEO-1;
    if (*pt_y<0) *pt_y = 0;
    if (*pt_y>=LNS) *pt_y = LNS-1;
  }
  else{
    /// periodic boundary conditions
    if (!rot_map){
      if (*pt_x<0) *pt_x+=LEO;
      if (*pt_x>=LEO) *pt_x-=LEO;
      if (*pt_y<0) *pt_y+=LNS;
      if (*pt_y>=LNS) *pt_y-=LNS;
    }
    else{
      /// rotating space
      int x = *pt_x + 1;
      int y = *pt_y + LN;
      assert((x >= 0) && (y >= 0));
      if (OutOfSpace(x, y)){
        *pt_x = RotMapPosX(x, y) - 1;
        *pt_y = RotMapPosY(x, y) - LN;
      }
    }
  }
}

// apply boundary conditions on interpolated cell with float coordinates
void apply_bc_cel_interp(float *pt_x, float *pt_y)
{
  if (!pbc_mode){
    /// open or closed boundaries
    if (*pt_x<0) *pt_x = 0;
    if (*pt_x>=LEO) *pt_x = LEO-1;
    if (*pt_y<0) *pt_y = 0;
    if (*pt_y>=LNS) *pt_y = LNS-1;
  }
  else{
    /// periodic boundary conditions
    if (!rot_map){
      if (*pt_x<0) *pt_x+=LEO;
      if (*pt_x>=LEO) *pt_x-=LEO;
      if (*pt_y<0) *pt_y+=LNS;
      if (*pt_y>=LNS) *pt_y-=LNS;
    }
    else{
      /// rotating space
      int x = floor(*pt_x) + 1;
      int y = floor(*pt_y) + LN;
      if (OutOfSpace(x, y)){
        *pt_x += RotMapPosX(x, y) - x;
        *pt_y += RotMapPosY(x, y) - y;
      }
    }
  }
}

//interpolation de l'altitude en un point intermediaire
float interpole_alti_cell(float x, float z)
{
  int x0, x1, z0, z1;
  int y00, y01, y10, y11;
  float tx0, tx1, tz0, tz1;

  x0 = floor(x);
  z0 = floor(z);
  x1 = x0 + 1;
  z1 = z0 + 1;

  //assert(x1 < LEO);
  tx1 = x - x0;
  tx0 = 1.0 - tx1;

  //assert(z1 < LNS);
  tz1 = z - z0;
  tz0 = 1.0 - tz1;

  y00 = Alti(x0,z0);

  apply_bc_cell(&x0, &z1);
  y01 = Alti(x0,z1);

  apply_bc_cell(&x1, &z0);
  y10 = Alti(x1,z0);

  apply_bc_cell(&x1, &z1);
  y11 = Alti(x1,z1);

  return (float)(y00*tx0*tz0 + y01*tx0*tz1 + y10*tx1*tz0 + y11*tx1*tz1);
}


// calcule de la normale 3d au point (i,j)
// avec pente calculee suivant les altitudes interpolees dans la direction de la plus grande pente
void calcule_norm3d_cel_pente_max_interp(int i, int j, float r)
{
  static char start=1;
  //static float r;
  int r0;
  Vec3 *n3d;
  int x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4;
  int ux, uy, uz, vy, vz;
  int wx, wy, wz;
  float d, co, si, angle, r1, r2;
  float fx1, fx2, fy1, fy2, fz1, fz2, fux, fuy, fuz, fwx, fwy, fwz;

  if (!alti){
    ErrPrintf("ERROR: calcule_norm3d_cel - alti = NULL\n");
    exit(-1);
  }

  if (!norm3d){
    LogPrintf("calcul des normales 3d\n");
    AllocMemoryPrint("norm3d", norm3d, Vec3, LEO*LNS);
    ResetMemory(norm3d, Vec3, LEO*LNS);
  }

  /// calcul de la direction de plus grande pente
  //r = 3; //3.5;//+1.5*drand48(); //3;
  r0 = 2;//floor(r);
  if (start){
    //r = 3.0/tan(PI*ava_angle/180);
    //r = 3.5;
    LogPrintf("calcule_norm3d_cel_pente_max_interp : r = %f   r0 = %d\n", r, r0);

    start = 0;
  }

  n3d = norm3d + (i+j*LEO);

  /// U = tangente dans le plan XY
  x1 = i-r0;
  z1 = j;
  x2 = i+r0;
  z2 = j;
  apply_bc_cell(&x1, &z1);
  apply_bc_cell(&x2, &z2);

  ux = (pbc_mode) ? 2*r0 : x2-x1;
  y1 = Alti(x1,z1);
  y2 = Alti(x2,z2);
  uy = y2-y1;

#ifndef STABILITY_ANALYSIS
  /// V = tangente dans le plan YZ
  if (LNS>=5){
    x3 = i;
    z3 = j-r0;
    x4 = i;
    z4 = j+r0;
    apply_bc_cell(&x3, &z3);
    apply_bc_cell(&x4, &z4);

    y3 = Alti(x3,z3);
    y4 = Alti(x4,z4);
    vy = y3-y4;
    vz = (pbc_mode) ? -2*r0 : z3-z4;
  }
  else {
    vz = 0;
  }
#else
  vz = 0;
#endif

  /// normale W = U ^ V
  if (vz){
    n3d->x = uy*vz/NB_MVT_EO;
    n3d->z = ux*vy;
  }
  else{
    n3d->x = -uy/NB_MVT_EO;
    n3d->z = 0.0;
  }

  /// direction de plus grande pente
  angle = (n3d->x)? atan(n3d->z/n3d->x) : PI/2;

  /// recalcul de la normale dans cette direction avec interpolation de l'altitude
  co = cos(angle);
  si = sin(angle);
  r1 = co*r;
  r2 = si*r;
  if (LNS<2*r2){// tangente U
    r1 = r;
    r2 = 0;
  }
  //assert((r1 >= 0) && (r2 >= 0));
  //LogPrintf("r1 = %d   r2 = %d\n", r1, r2);

  /// coordonnees des points tests
  fx1 = i-r1;
  fz1 = j-r2;
  fx2 = i+r1;
  fz2 = j+r2;
  apply_bc_cel_interp(&fx1, &fz1);
  apply_bc_cel_interp(&fx2, &fz2);

  fux = (pbc_mode) ? 2*r1 : fx2-fx1;
  fuz = (pbc_mode) ? 2*r2 : fz2-fz1;

  fy1 = interpole_alti_cell(fx1, fz1);
  fy2 = interpole_alti_cell(fx2, fz2);
  fuy = fy2-fy1;

  /// normale W
  fwx = -fuy*co;
  fwy = fuz*si + fux*co;
  fwz = -fuy*si;
  // normalisation
  d = sqrt(fwx*fwx+fwy*fwy+fwz*fwz);
  if (fwy<0) d = -d;
  if (!d) LogPrintf("ux = %d   uy = %d   uz = %d   wx = %d   wy = %d   wz = %d   d = %f\n",ux,uy,uz,wx,wy,wz,d);
  n3d->x = fwx/d;
  n3d->y = fwy/d;
  n3d->z = fwz/d;
  //if (!(n3d->y > 0.0)) {LogPrintf("n3d->y = %f\n",n3d->y);}
  assert(n3d->y > 0.0);
}


void calcule_normales()
{
  int x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4;
  int ux, uy, vy, vz;
  int wx, wy, wz;
  float d;
  int i, j, ii;

  if (!norm2d){
    LogPrintf("calcul des normales 2d et 3d\n");
#ifdef NORM2D
    LogPrintf("mode NORM2D : CGV with 2d normals\n");
#endif
    AllocMemory(norm2d, Vec2, L*D);
    ResetMemory(norm2d, Vec2, L*D);
    AllocMemory(norm3d, Vec3, L*D);
    ResetMemory(norm3d, Vec3, L*D);
  }

  if (!alti){
    ErrPrintf("ERROR: calcule_normales - alti=NULL\n");
    exit(-1);
  }

  for(ii=0, j=0; j < LNS; j++)
    for(i=0; i < LEO; i++, ii++){
      if (rot_map && OutOfSpace(1+i,LN+j)) continue;
      // calcul de la normale 2d (dans un plan vertical XY)
      x1 = i-2;
      z1 = j;
      x2 = i+2;
      z2 = j;
      apply_bc_cell(&x1, &z1);
      apply_bc_cell(&x2, &z2);

      y1 = Alti(x1, z1);
      y2 = Alti(x2, z2);
//#ifdef NORM2D
      ux = (pbc_mode) ? 4 : x2-x1;
      uy = y2-y1;
      Check_min_max(uy,-5,5);
      d = sqrt(ux*ux+uy*uy);
      norm2d[ii].x = -uy/d;
      norm2d[ii].y = ux/d;
//#else
      // calcul de la normale 3d
      if (LNS>=5){
        x3 = i;
        z3 = j-2;
        x4 = i;
        z4 = j+2;
        apply_bc_cell(&x3, &z3);
        apply_bc_cell(&x4, &z4);
        vz = (pbc_mode) ? -4 : z3-z4;
        y3 = Alti(x3, z3);
        y4 = Alti(x4, z4);
        vy = y3-y4;
        Check_min_max(vy,-5,5);
      }
      else {
        vz = 0;
      }
      // W = U ^ V
      if (vz){
        wx = uy*vz;
        wy = - ux*vz;
        wz = ux*vy;
        d = sqrt(wx*wx+wy*wy+wz*wz);
        norm3d[ii].x = wx/d/NB_MVT_EO;
        norm3d[ii].y = wy/d;
        norm3d[ii].z = wz/d;
      }
      else{
        norm3d[ii].x = norm2d[ii].x;
        norm3d[ii].y = norm2d[ii].y;
        norm3d[ii].z = 0.0;
      }
      //LogPrintf("i=%d   j=%d   nx=%f   ny=%f   nz=%f   tgn_xy=%f\n", i, j, norm3d[ii].x, norm3d[ii].y, norm3d[ii].z, sqrt(norm3d[ii].x*norm3d[ii].x+norm3d[ii].y*norm3d[ii].y));
//#endif
    }
}

#ifdef AVALANCHES
// callback de controle des transitions d'avalanches en fonction de la pente locale en surface
int check_ava(int ix, void *data)
{
  static char start=1;
  int x, y, z, i, j, ij, i2, j2, al, al2, delta_x, delta_y, delta_z, chk_ok=0;
  float nx, ny, nz;
  static float ny_lim=0;
  //static float nxz_lim=0.7071; //sqrt(2)/2
  float chk_angle=0, r, alea, prob;
  Vec3 *pt_n3d=NULL;
  DataCheck *pdc = data;
  unsigned char typ;

  if (start){
    LogPrintf("check_ava : ava_h_lim = %d\n", ava_h_lim);
    LogPrintf("check_ava : ava_angle = %f\n", ava_angle);
#ifdef CELL_COLOR
    LogPrintf("check_ava : ava_angle_col = %f\n", ava_angle_col);
#endif
    if (!ava_upwind) LogPrintf("check_ava: no avalanche upwind\n");
    start=0;
  }

  //calcul de la position (x,y,z)
  Calcule_xyz(ix, x, y, z);
  i = x-1;
  j = z-LN;
  ij = i + j*LEO;

  //type de cellule
  typ = TE[ix].celltype;

  //altitude
  al = alti[ij];
  if (al+y < H){ //grain qui tombe
    //LogPrintf("GRAIN QUI TOMBE : al = %d   y = %d   al+y = %d\n", al, y, al+y);
    return 0;
  }

  //direction d'avalanche
  if (pdc->dir == EST){ //est-ouest
    delta_x = (pdc->cel == 1) ? 1 : -1;
    delta_z = 0;
  }
  else{ //nord-sud
    delta_x = 0;
    delta_z = (pdc->cel == 1) ? 1 : -1;
  }
  i2 = i+delta_x;
  j2 = j+delta_z;
#ifdef CYCLAGE_HOR
  apply_bc_cell(&i2, &j2);
#endif

  //calcul de la difference de hauteur
  al2 = Alti(i2, j2);
  delta_y = al - al2;

  //if (delta_y <= 0) return 0;
  if (delta_y < 0) {
    int typ2 = TE[ix+delta_x+delta_z*HL].celltype;
    LogPrintf("check_ava : x = %04d   z = %04d   al = %04d   delta_x = %04d   delta_y = %04d   delta_z = %04d   typ = %d   typ2 = %d   dir = %d   cel = %d\n", x, z, al, delta_x, delta_y, delta_z, typ, typ2, pdc->dir, pdc->cel);
  }

#ifdef CELL_COLOR
  chk_angle = (TE[ix].color) ? ava_angle_col : ava_angle;
#else
  chk_angle = ava_angle;
#endif

  if (chk_angle){
    //r = (chk_angle > 33.0) ? 3.0 : 3.5;
    r = (chk_angle <= 33.0) ? 3.5 : (chk_angle <= 45.0) ? 3.0 : 1.5;
    calcule_norm3d_cel_pente_max_interp(i, j, r);
    //calcule_norm3d_cel_iso(x, z);
    pt_n3d = norm3d + ij;
    //calcule_norm3d_cel_iso(i-1, j);
    //pt_n3d = norm3d_iso + ij;
  }

  if (!ava_upwind && chk_angle && (pt_n3d->x <= 0)) return 0; //transition blocked when slope is windward

  if (delta_y > ava_h_lim) return 1; //check the height difference

  if (chk_angle){
    //vertical component of the normal vector
    ny = pt_n3d->y;

    //max value
    ny_lim = 0.999*cos(PI*ava_angle/180.0); //coefficient 0.999 a cause des imprecisions de calcul

    chk_ok = (ny < ny_lim);
  #ifdef CELL_COLOR
    float ny_lim2 = 0.999*cos(PI*ava_angle_col/180.0); //coefficient 0.999 a cause des imprecisions de calcul
    if (TE[ix].color) chk_ok = (ny < ny_lim2);
  #endif

    // loi de proba lineaire avec seuil
    /*alea = drand48();
    prob = (ny < ny_lim) ? (ny_lim - ny)/ny_lim : 0.0;
    chk_ok = (prob > alea);*/

    //if (chk_ok){
    if (0){ // disabled condition (for test only)
      //float ps, nxz;
      nx = pt_n3d->x; //composante EO
      nz = pt_n3d->z;  //composante NS

      //on verfie que la chute a lieu dans la direction de plus grande pente ( +-PI/4)
      //chk_ok = (nx*delta_x + nz*delta_z >= nxz_lim*sqrt(nx*nx+nz*nz));

      //loi de proba fonction de la direction
      //nxz = sqrt(nx*nx+nz*nz);
      float nxz2 = nx*nx+nz*nz;
      if (nxz2){
        alea = drand48();
        //prob = (pdc->dir == EST) ? abs(nx)/nxz : abs(nz)/nxz; //bug
        //prob = (pdc->dir == EST) ? fabs(nx)/nxz : fabs(nz)/nxz;
        prob = (pdc->dir == EST) ? nx*nx/nxz2 : nz*nz/nxz2;
        //prob = (pdc->dir == EST) ? 1 : 0.5;
        //prob = (1 + prob)/2.0; // mix entre hasard et determinisme ( prob entre 0.5 et 1 )
        prob = (1 + 3*prob)/4.0; // mix entre hasard et determinisme ( prob entre 0.25 et 1 )
        //assert(prob <= 1.0);
        chk_ok = (prob > alea);
      }
      else
        chk_ok = 1;

      /*
      // loi de proba fonction de la pente lineaire avec seuil
      alea = drand48();
      //prob = (ny < ny_lim) ? (ny_lim - ny)/ny_lim : 0.0;
      ps = nx*delta_x + nz*delta_z;
      nxz = sqrt(nx*nx+nz*nz);
      //ps_lim = nxz_lim*nxz;
      //prob = (ps > 0) ? ps/nxz : 0.0;
      prob = (ps > 0) ? ps : 0.0;
      chk_ok = (prob > alea);
      */
      //LogPrintf("chk_ok = %d   nxz = %f\n", chk_ok, (nx*delta_x + nz*delta_z)/sqrt(nx*nx+nz*nz));

      //if (!chk_ok) LogPrintf("blocage avalanche (prob=%f)\n", prob);
    }
  }

  //if (chk_ok) {LogPrintf("i=%d   ny=%.12f   ny_lim=%.12f\n", i, ny, ny_lim); }

  if (ava_mode == AVA_PROPAG){
#if 0
    // delai variable
    if (chk_ok){
      //float ava_propag_prob = 0.00001;
      float ava_propag_prob = 0.00005;
      alea = drand48();
      if (alea < ava_propag_prob){
        //on declenche une grosse avalanche avec propagation
        LogPrintf("appel loop_ava_propag() : i=%d   j=%d   alea=%f   csp_time=%g\n", i, j, alea, csp_time);
        loop_ava_propag(i, j);
        if (nb_ava_propag) chk_ok = 0;
      }
    }
#endif
#if 1
    // delai fixe
    if (chk_ok){
      static float time_threshold = 0.0;
      //float delai = 200; //modele AVA
      float delai = ava_delay;
      if (!time_threshold) time_threshold = delai;
      if (csp_time > time_threshold){
        //on declenche une grosse avalanche avec propagation
        LogPrintf("appel loop_ava_propag() : i=%d   j=%d   csp_time=%g\n", i, j, csp_time);
        LogPrintf("delai = %f\n", delai);
        loop_ava_propag(i, j);
        if (nb_ava_propag) chk_ok = 0;
        //time_threshold += delai;
        time_threshold = csp_time + delai;
      }
    }
#endif
  }

#ifdef CELL_COLOR
  //coloration des cellules qui tombent du fait d'une avalanche
  //if (chk_ok) TE[ix].color = 1; //else TE[ix].color = 0;
#endif

  return chk_ok;
}
#endif //AVALANCHES

// callback de controle des transitions de montee de grains en fonction de la pente locale en surface
// (on suppose que le lien de transitions se fera en est-ouest)
// TODO? : on pourrait aussi verifier la difference de hauteur en est-ouest avant de calculer la normale ...
int check_slope(int ix, void *data)
{
  static char start=1;
  int x, y, z, i, j, ij, chk_ok=0;
  float nx;
  static float nx_lim=0;
  Vec3 *pt_n3d;

  if (start){
    LogPrintf("check_slope : slope_angle = %f\n", slope_angle);
    start=0;
  }

  //calcul de la position (x,y,z)
  Calcule_xyz(ix, x, y, z);
  i = x-1;
  j = z-LN;

  i+=2; // lire la normale en aval
  if (!pbc_mode){
    if (i >= LEO) i = LEO-1;
  }
#ifdef CYCLAGE_HOR
  else if (!rot_map){
    if (i >= LEO) i -= LEO;
  }
  else{
    apply_bc_cell(&i, &j);
  }
#endif

  //altitude
  ij = i + j*LEO;
  /*
  al = alti[ij];
  if (al+y < H){ //grain qui tombe
    //LogPrintf("GRAIN QUI TOMBE : al = %d   y = %d   al+y = %d\n", al, y, al+y);
    return 0;
  }*/

  //calcul de la normale 3d
  //calcule_norm3d_cel(i, j);
  calcule_norm3d_cel_pente_max_interp(i, j, 3.0);
  pt_n3d = norm3d + ij;
  nx = pt_n3d->x;
  //ny = pt_n3d->y;

  //critere
  //angle = PI/6;
  //ny_lim = cos(angle);
  if (!nx_lim) nx_lim = -sin(slope_angle*PI/180.0);
  //LogPrintf("nx_lim = %f\n", nx_lim);

  //chk_ok = (ny > ny_lim);
  chk_ok = (nx > nx_lim);
  //chk_ok = (nx < 0) && (nx > nx_lim);
  //chk_ok = (nx < 0);

  return chk_ok;
}




#ifdef CGV
float dist_grdv = 10.0; //10.0  //distance caracteristique pour le calcul du gradient de la vitesse
float dist0_grdv = 0.0; //0.0  //distance de reference pour le calcul du gradient de la vitesse
float grdvc = 0.0; //constante critique pour le controle de l'erosion en fonction du gradient de vitesse
float grdvc_max = 0.0; //seuil max pour le controle de l'erosion en fonction du gradient de vitesse
float grdvc_min = 0.0; //seuil min pour le controle de l'erosion en fonction du gradient de vitesse
#ifdef CELL_COLOR
float grdvc_max_col = -1.0; //seuil max pour le controle de l'erosion des grains colores en fonction du gradient de vitesse
float grdvc_min_col = -1.0; //seuil min pour le controle de l'erosion des grains colores en fonction du gradient de vitesse
#endif
float cgv_coef = 0.0; //valeur moyenne de probabilite de blocage d'une transition

extern float *Velx_interp, *Vely_interp;
extern float maxvel, meanvel;

void reset_grad_vel()
{
  if (grdv) ResetMemory(grdv, float, L*D);
}

void calcule_grad_vel()
{
  //int i, j, x0, y0, x1, y1, x2, y2, ix1, ix2, vx1, vy1, vx2, vy2;
  int i, j, ii, x0, y0, z0, x2, y2, z2, ix2;
  float vx2, vy2, vel2, tgx, tgy;//, sum_gv, gv_mean;
  //float *pt_gv;

  //LogPrintf("calcule_grad_vel : debut\n");

  if (!grdv){
    //grdv = (float *)malloc((L-2)*LNS*sizeof(float));
    AllocMemory(grdv, float, L*D);
    ResetMemory(grdv, float, L*D);
    LogPrintf("calcule_grad_vel : dist0_grdv = %f   dist_grdv = %f\n", dist0_grdv, dist_grdv);
  }

  if (!alti || !norm2d || !Velx_interp){
    WarnPrintf("WARNING: cannot compute grdv (%lu,%lu,%lu)\n", (long)alti, (long)norm2d, (long)Velx_interp);
    return;
  }

  //sum_gv = 0.0;

  for(ii=j=0; j<LNS; j++)
    for(i=0; i<LEO; i++, ii++){
      if (rot_map && OutOfSpace(1+i,LN+j)) continue;
      x0 = i+1;
      y0 = H-alti[ii];
      z0 = j+LN;
      //vecteur tangent a la surface en 2d
      tgx = norm2d[ii].y;
      tgy = -norm2d[ii].x;
#ifdef NORM2D
      x2 = x0 + roundf(norm2d[ii].x*dist_grdv);
      y2 = y0 - roundf(norm2d[ii].y*dist_grdv);
      z2 = z0;
#else
      x2 = x0 + roundf(norm3d[ii].x*dist_grdv);
      y2 = y0 - roundf(norm3d[ii].y*dist_grdv);
      z2 = z0 + roundf(norm3d[ii].z*dist_grdv);
#endif
      /*if ((x2<0) || (x2>=L) || (y2<0) || (y2>=H)){
        *pt_gv=0.0;
        continue;
      }*/
      if (x2<0) x2=0;
      if (x2>=L) x2=L-1;
      if (y2<0) y2=0;
      //if (y2<=h_plaf) y2=h_plaf+1; //plafond
      if (y2>=H) y2=H-1;
      if (z2<1) z2=1;
      if (z2>=D-1) z2=D-2;
      ix2 = z2*HL + y2*L + x2;
      if (rot_map && pbc_mode && OutOfSpace(x2,z2)) ix2 = hcycle(ix2);
      if ((ix2<0) || (ix2>=HLD)){
        ErrPrintf("WARNING: calcule_grad_vel, i=%d, j=%d, ix2=%d\n", i, j, ix2);
        //*pt_gv=0.0;
        continue;
      }
      vx2 = Velx_interp[ix2];
      vy2 = Vely_interp[ix2];
      //vel2 = sqrt(vx2*vx2+vy2*vy2); //norme du vecteur vitesse
      vel2 = vx2*tgx + vy2*tgy; //projection du vecteur vitesse sur la tangente a la surface
      grdv[ii] = vel2;
      //if (y2 <= h_plaf) *pt_gv = grdvc_max; //gradient max quand la normale traverse le plafond
      //sum_gv += vel2;
      //LogPrintf("x0=%d   vx2=%d   vy2=%d   gv=%f\n",x0,vx2,vy2,vel2);
    }

    //gv_mean = sum_gv / ii;
    //cgv_coef = prob_cgv(gv_mean);

  //LogPrintf("calcule_grad_vel : fin\n");
}

//initialisation de la fonction de probabilite
/*
void init_cgv()
{
  int ix, vx, vy, vv;
  maxvel = 0.0;
  meanvel = 0.0;
  for (ix=0; ix<HLD; ix++){
    vx = Velx_interp[ix];
    vy = Vely_interp[ix];
    vv = vx*vx+vy*vy;
    if (vv>maxvel) maxvel=vv;
    meanvel += vv;
  }
  maxvel=sqrt(maxvel);
  meanvel=sqrt(meanvel/HLD);
  LogPrintf("init_prob_cgv : maxvel = %f\n", maxvel);
  LogPrintf("init_prob_cgv : meanvel = %f\n", meanvel);
}
*/

//fonction de probabilite pour le controle en fonction du gradient de la vitesse
float prob_cgv(float gv, float gvmin, float gvmax)
{
  //static float grdv_max = VSTEP_H*VSTEP_L*VSTEP_TIME*2;
  static char start = 1;
  //static float grdvc_min = 0.0;
  //static float grdvc = 0.0;
  //static float grdvc_max = 0.0;

  if (start){
/*#ifdef SOLFAST
    grdvc_min = 0; //0; //grdv_max/10; //grdv_max/10.0;
    grdvc = grdv_max/4; //grdv_max/4
    grdvc_max = grdv_max/2; //grdv_max/2; //grdv_max/10.0;
#else
    grdvc_min = 0; //grdv_max/20.0;
    grdvc = grdv_max/20; //grdv_max/10;
    grdvc_max = grdv_max/4; //grdv_max/10.0;
#endif*/
    //init_cgv();
    //if (!maxvel) maxvel = grdv_max;
    if (!grdvc) grdvc = maxvel/4; //grdv_max/4
    //grdvc = meanvel/10; //grdv_max/4
    if (!grdvc_min) grdvc_min = 0.0; //0; //grdv_max/10; //grdv_max/10.0;
    if (!grdvc_max) grdvc_max = maxvel; //grdv_max/2; //grdv_max/10.0;
    //grdvc_max = meanvel; //grdv_max/2; //grdv_max/10.0;

    //LogPrintf("grdv_max = %f\n", grdv_max);
    LogPrintf("grdvc_min = %f\n", grdvc_min);
    LogPrintf("grdvc = %f\n", grdvc);
    LogPrintf("grdvc_max = %f\n", grdvc_max);
    gvmin = grdvc_min;
    gvmax = grdvc_max;
    //dump_cgv();
    start=0;
  }

  //float prob = (gv>grdvc)? 1 : gv/grdvc;
  //float prob = (gv>grdvc)? 1 : pow((gv/grdvc),1.5);
  //float prob = 1.0-exp(-gv/grdvc);
  //float prob = 1.0-exp(-pow(fabs(gv)/grdvc,1.5));
  //float prob = (gv<grdvc)? 0 : 1.0-exp(-(gv-grdvc)/grdvc);
  //float prob = (gv<gvmin)? 0 : (gv>gvmax)? 1 : (gv-gvmin)/(gvmax-gvmin);
  float prob = ((gv<gvmin) || (gv<0))? 0 : (gv>gvmax)? 1 : (gv-gvmin)/(gvmax-gvmin);
  //float prob = (gv<grdvc_min)? 0 : (gv>grdvc_max)? 1 : pow((gv-grdvc_min)/(grdvc_max-grdvc_min),1.5);
  return prob;
}

// callback de controle en fonction du gradient de la vitesse selon la normale a la surface
int check_grad_vel(int ix, void *data)
{
  int x, y, z, chk_ok;
  float gv, alea, prob;

  //calcul de la position (x,y,z)
  Calcule_xyz(ix, x, y, z);
  //if (rot_map && OutOfSpace(x,z)) WarnPrintf("WARNING: OutOfSpace - x=%d - z=%d\n", x, z);
  gv = grdv[(x-1)+(z-LN)*LEO];
  alea = drand48();
  prob = prob_cgv(gv, grdvc_min, grdvc_max);

  //calcul de la pente
  /*int i = x-1;
  int j = z-LN;
  int ij = i + j*LEO;
  calcule_norm3d_cel_pente_max_interp(i, j, 3.0);
  Vec3 *pt_n3d = norm3d + ij;
  float nx = pt_n3d->x;
  if (nx){
    float slope = fabsf(asinf(nx));
    float slope_min = 1/6.0;
    if (slope < slope_min) slope = slope_min;
    prob = prob*slope_min/slope;
  }*/

  chk_ok = (prob > alea);
  //if (!chk_ok) LogPrintf("check_grad_vel: x=%d   H-y=%d   alti=%d   gv=%f\n", x, H-y, alti[(x-1)+(z-LN)*LEO], gv);
  //chk_ok = (z < L/2) ? 1 : chk_ok; //test
  //if (!chk_ok) {LogPrintf(" check_grad_vel : not ok\n");}
  return chk_ok;
}

#ifdef CELL_COLOR
// callback de controle en fonction du gradient de la vitesse selon la normale a la surface
int check_grad_vel_color(int ix, void *data)
{
  int x, y, z, chk_ok;
  float gv, alea, prob;

  //calcul de la position (x,y,z)
  Calcule_xyz(ix, x, y, z);
  gv = grdv[(x-1)+(z-LN)*LEO];
  alea = drand48();
  prob = (!TE[ix].color) ? prob_cgv(gv, grdvc_min, grdvc_max) : prob_cgv(gv, grdvc_min_col, grdvc_max_col);
  chk_ok = (prob > alea);
  return chk_ok;
}
#endif

void compute_coef_cgv(/*int idb*/)
{
  int i, k, ik;
  int n0, nn, np;
  float sum_gv, gv_mean;
  static char start = 1;
  FILE *fp;

  if (!grdv) return /*0.0*/;
  if (start){
    LogPrintf("calcul coef_cgv pour correction du pas de temps\n");
    start = 0;
  }

  //estimation du coefficient correctif sur le pas de temps (dans le cas d'une loi de probabilite lineaire a seuils)
  sum_gv = 0.0;
  n0 = np = nn = 0;
  //for (i=0; i<(L-2)*LNS; i++){
  for (k=0; k<LNS; k++){
    if (!psol[LN+k]){
      ik = k*LEO;
      for (i=0; i<LEO; i++, ik++){
        if (rot_map && OutOfSpace(1+i,LN+k)) continue;
        if (grdv[ik] > grdvc_max)
          np++;
        else if (grdv[ik] > grdvc_min){
          n0++;
          sum_gv += grdv[ik];
        }
        nn++;
      }
    }
  }

  gv_mean = n0 ? (sum_gv / n0) : 0.0;
  //cgv_coef = nb_pv ? (prob_cgv(gv_mean,grdvc_min,grdvc_max)*n0 + np)/(nb_pv*LEO) : 0.0;
  cgv_coef = nb_pv ? (prob_cgv(gv_mean,grdvc_min,grdvc_max)*n0 + np)/nn : 0.0;
  //cgv_coef = 1.0;

  assert(cgv_coef<=1.0);
}
#endif //CGV


#ifdef AVALANCHES
//avalanches dans 8 directions equiprobables
void avalanches(unsigned char typ, short h_lim, short nb_cel_max, char alti_mode)
{
  //static short *alti = NULL;
  static char first=1;
  Cell *aux, *aux2, cel;
  //CellData data;
  short *pt, *pt2;
  int i, j, k, jmin, jmax, ix, ix2;
  unsigned char alea, cpt, flag_ava;
  int step[8], step_alti[8];
  short i_cel, h_lim_diag;
  int nb_ava, nb_ava_iter;
  FILE *fp;

#ifdef MODEL_AVA
//#if 1
  h_lim_diag = ceil(h_lim*sqrt(2)); //-> octogone
#else
  h_lim_diag = h_lim;
#endif

//#ifdef MODEL_AVA
//#ifdef CELL_COLOR
#if 0
  short h_lim0 = h_lim;
  short h_lim1 = h_lim-1;
  short h_lim_diag0 = h_lim_diag;
  short h_lim_diag1 = h_lim_diag-1;
  if (h_lim1==0) h_lim1 = 1;
  assert(h_lim1 > 0);
#endif
//#endif

  if (first){
//#ifdef CELL_COLOR
#if 0
    LogPrintf("avalanches : h_lim0 = %d\n", h_lim0);
    LogPrintf("avalanches : h_lim1 = %d\n", h_lim1);
#else
    LogPrintf("avalanches : h_lim = %d\n", h_lim);
#endif
    LogPrintf("avalanches : h_lim_diag = %d\n", h_lim_diag);
    if (pbc_mode && rot_map){
      ErrPrintf("ERROR: sorry, rotating space is not compatible with synchronous avalanches\n");
      exit(-1);
    }
    first = 0;
  }

  //calcule_alti(typ, alti_mode);

  //avalanches
  nb_ava = nb_ava_iter = 0;
  jmin=0;
  jmax=LNS-1;
  /*if (LNS<L-2){
    //couloir nord-sud
    jmin=0;
    jmax=LNS-1;
  }
  else{
    //pas d'avalanches sur les bords
    LogPrintf("pas de couloir LNS=%d L=%d\n", LNS, L);
    jmin=1;
    jmax=L-2;
  }*/
  //LogPrintf("h_lim_diag = %d\n", h_lim_diag);
  do{
    nb_ava_iter++;
    flag_ava = 0;
    pt = alti+jmin*LEO;
    for(j=jmin; j <= jmax; j++){
      for(i=0; i < LEO; i++, pt++){
        // on compare l'elevation locale dans les 8 directions
        cpt = 0;
//#ifdef MODEL_AVA
//#ifdef CELL_COLOR
#if 0
        //h_lim depend de la couleur de la cellule la plus haute
        ix = i+(H-(*pt))*L+(LN+j)*HL; //cellule du sommet
        h_lim = (TE[ix].color) ? h_lim1 : h_lim0;
        h_lim_diag = (TE[ix].color) ? h_lim_diag1 : h_lim_diag0;
#endif
//#endif
        if ((i>0) && (*pt - *(pt-1) > h_lim)) {step_alti[cpt] = -1; step[cpt++] = -1;}  //gauche
        if ((i<LEO-1) && (*pt - *(pt+1) > h_lim)) {step_alti[cpt] = 1; step[cpt++] = 1;}  //droite
        if ((j>jmin) && (*pt - *(pt-LEO) > h_lim)) {step_alti[cpt] = -LEO; step[cpt++] = -HL;}  //derriere
        if ((j<jmax) && (*pt - *(pt+LEO) > h_lim)) {step_alti[cpt] = LEO; step[cpt++] = HL;}  //devant
        if ((i>0) && (j>jmin) && (*pt - *(pt-LEO-1) > h_lim_diag)) {step_alti[cpt] = -LEO-1; step[cpt++] = -HL-1;}  //diagonale gauche+derriere
        if ((i<LEO-1) && (j>jmin) && (*pt - *(pt-LEO+1) > h_lim_diag)) {step_alti[cpt] = -LEO+1; step[cpt++] = -HL+1;}  //diagonale droite+derriere
        if ((i>0) && (j<jmax) && (*pt - *(pt+LEO-1) > h_lim_diag)) {step_alti[cpt] = LEO-1; step[cpt++] = HL-1;}  //diagonale gauche+devant
        if ((i<LEO-1) && (j<jmax) && (*pt - *(pt+LEO+1) > h_lim_diag)) {step_alti[cpt] = LEO+1; step[cpt++] = HL+1;}  //diagonale droite+devant
#ifdef CYCLAGE_HOR
        if (pbc_mode){
          if ((i==0) && (*pt - *(pt+LEO-1) > h_lim)) {step_alti[cpt] = LEO-1; step[cpt++] = L-3;}  //gauche
          if ((i==LEO-1) && (*pt - *(pt-LEO+1) > h_lim)) {step_alti[cpt] = -(LEO-1); step[cpt++] = -(L-3);}  //droite
          //if (LNS==L-2){ //cyclage nord-sud
            if ((LN+j==0) && (*pt - *(pt+(LNS-1)*LEO) > h_lim)) {step_alti[cpt] = (LNS-1)*LEO; step[cpt++] = (D-3)*HL;}  //derriere
            if ((LN+j==D-2) && (*pt - *(pt-(LNS-1)*LEO) > h_lim)) {step_alti[cpt] = -(LNS-1)*LEO; step[cpt++] = -(D-3)*HL;}  //devant
            //printf("coucou i= %d   j=%d   cpt=%d\n",i,j,cpt); fflush(stdout);
          //}
        }
#endif
        if (cpt){
          // choix d'une direction pour l'avalanche
          //LogPrintf("avalanche : i = %d   j = %d   cpt = %d\n", i, j, cpt);
          for (i_cel=0; i_cel<nb_cel_max; i_cel++){
            alea = floor(cpt*drand48()); //LogPrintf("alea = %d   step_alti = %d\n", alea, step_alti[alea]);
            pt2 = pt + step_alti[alea];
            if ((pt2<alti) || (pt2>alti+LEO*LNS)){ErrPrintf("ERROR: avalanche, array overflow in alti\n");exit(-1);} //antibug
            k = H-(*pt);
            ix = 1+i+k*L+(LN+j)*HL; //cellule du sommet
            //if ((ix<0) || (ix >= HLD)) {ErrPrintf("ERREUR: avalanche - ix = %d\n", ix);  exit(-1);}
            aux = TE+ix;
            cel = *aux;
//#if CELL_DATA
//            get_cell_data(data, ix);
//#endif
            if (aux->celltype == DUM) break;
            if (aux->celltype != typ) {
              ErrPrintf("WARNING: avalanche, incorrect type : %d (i = %d  j = %d  k = %d  deniv = %d  step_alti = %d  nb_ava_iter = %d)\n", aux->celltype, i, j, k, *pt - *pt2, step_alti[alea], nb_ava_iter);
              //LogPrintf("type  : %d (i = %d  j = %d  k = %d  nb_ava_iter = %d)\n", (aux-L)->celltype, i, j, k-1, nb_ava_iter);
              //LogPrintf("type  : %d (i = %d  j = %d  k = %d  nb_ava_iter = %d)\n", (aux+L)->celltype, i, j, k+1, nb_ava_iter);
              //LogPrintf("nb_ava = %d\n", nb_ava);
              break;
            }
            if ((*pt)>(*pt2)){
            //if (1){
              ix2 = ix + step[alea];  //cellule voisine
              //aux2 = aux + ix2;
              //LogPrintf("i = %d  j = %d  k = %d  deniv = %d  step_alti = %d  cpt = %d  alea = %d\n", i, j, k, *pt - *pt2, step_alti[alea], cpt, alea);
              //on decale la cellule voisine vers le sommet (pour conservation de la matiere)
              //aux->celltype = aux2->celltype;
              //ajoute_cellule(TE[ix2].celltype, ix);
              deplace_cellule(ix2, ix);
              //on fait remonter la colonne de cellules
              //while (((aux2+L)->celltype != typ) && ((aux2+L)->celltype != DUM) && ((aux2+L)->celltype != BORD)){
              for (k=*pt; k>(*pt2)+1; k--){
                //aux2->celltype = (aux2+L)->celltype;
                ix = ix2;
                ix2+=L;
                //ajoute_cellule(TE[ix2].celltype, ix);
                deplace_cellule(ix2, ix);
              }
            }
            else{
              //cas ou la cellule "haute" est plus basse que la cellule "basse" !!!
              //on echange simplement les cellules "haute" et "basse"
              //ix2 = ix + step[alea] + ((*pt)-(*pt2)-1)*L; //cellule "basse"
              //LogPrintf("cellule basse : %d   passe : %d   k1 = %d   k2 = %d\n", TE[ix2].celltype, nb_ava_iter, (*pt), (*pt2));
              //ajoute_cellule(TE[ix2].celltype, ix);
              continue;
            }
            // on fait tomber la cellule haute dans la direction choisie
            //aux2->celltype = typ;
            //ajoute_cellule(typ, ix2);
            init_cellule(cel, ix2);
            //mise-a-jour des elevations locales
            (*pt)--;
            (*pt2)++;
            if (alti_mode == ALTI_MODE_HAUT){ //recalcul de l'altitude
              aux += L;
              //while ((aux->celltype != typ) && (aux->celltype != DUM) && (aux->celltype != BORD)) {(*pt)--; aux += L;}
              while (Phase[aux->celltype] != SOLID) {(*pt)--; aux += L;}
            }
            flag_ava = 1;
            nb_ava++;
          }
        }
      }
    }
  }while (flag_ava);

  //dump avalanches
  fp = fopen("AVA.log","a");
  //fprintf(fp,"avalanches : %d   passes : %d   temps : %f\n", nb_ava, nb_ava_iter, temps); //fflush(stdout);
  fprintf(fp,"avalanches : %d   passes : %d   csp_time : %f   h_lim : %d\n", nb_ava, nb_ava_iter, csp_time, h_lim); //fflush(stdout);
  fclose(fp);

  if (nb_ava){
    //recalcul des tableaux de doublets actifs
    init_db_pos();
    //recalcul des cellules
    init_Ncel();

#ifdef PARALLEL
  //synchronisation des zones de recouvrement
  if (mode_par) synchro_tunnel_zones(1);
#endif
  }

}


//avalanches dans 8 directions selon les normales
//void avalanches_norm(unsigned char typ, float angle, short nb_cel_max, char alti_mode)
void avalanches_norm(unsigned char typ, short nb_cel_max, char alti_mode)
{
  //static short *alti = NULL;
  static char first=1;
  static unsigned char flag_diag=1;
  Cell *aux, *aux2, cel;
  //Vec3 n3d;
  //CellData data;
  short *pt, *pt2;
  float nx, ny, nz, ny_lim, anh, r;
  int i, j, k, jmin, jmax, ix, ix2, ij;
  unsigned char flag_ava, flag_drop;
  int step, step_alti;
  short i_cel;
  int nb_ava, nb_ava_iter, nb_ava_drop;
  FILE *fp;
  Vec3 *pt_n3d;
  float angle;

  //calcule_alti(typ, alti_mode);
  //calcule_alti_mean();

  fp = fopen("AVA.log","a");
  angle = ava_angle;
  ny_lim = cos(angle*PI/180);
  r = (angle <= 33.0) ? 3.5 : 3.0;
  //fprintf(fp, "avalanches_norm : ny_lim = %f\n", ny_lim);

#ifdef CELL_COLOR
//#if 0
  float angle0 = angle;
  float angle1 = ava_angle_col; //angle - 15.0*PI/180.0;
  float ny_lim0 = cos(angle0*PI/180);
  float ny_lim1 = cos(angle1*PI/180);
  assert(angle1 > 0);
#endif

  if (first){
#ifdef CELL_COLOR
//#if 0
    LogPrintf("avalanches_norm : ava_angle0 = %f \n", angle0);
    LogPrintf("avalanches_norm : ava_angle1 = %f \n", angle1);
#else
    LogPrintf("avalanches_norm : ava_angle = %f    r = %f\n", angle,r);
#endif
    if (pbc_mode && rot_map){
      ErrPrintf("ERROR: sorry, rotating space is not compatible with synchronous avalanches\n");
      exit(-1);
    }
    first = 0;
  }

  //avalanches
  nb_ava = nb_ava_iter = nb_ava_drop = 0;
  jmin=0;
  jmax=LNS-1;
  do{
    nb_ava_iter++;
    //fprintf(fp, "avalanches_norm : passe %d\n", nb_ava_iter);
    flag_ava = flag_drop = 0;
    //flag_diag = 1 - flag_diag;
    //angle_norm3d_iso = drand48()*PI/2.0;
    ij = 0;
    //calcule_normales();
    for(j=jmin; j <= jmax; j++){
      for(i=0; i < LEO; i++, ij++){
#ifdef CELL_COLOR
//#if 0
        //angle depend de la couleur de la cellule la plus haute
        pt = alti+ij;
        ix = 1+i+(H-(*pt))*L+(LN+j)*HL; //cellule du sommet
        angle = (TE[ix].color) ? angle1 : angle0;
        ny_lim = (TE[ix].color) ? ny_lim1 : ny_lim0;
#endif
        //calcule_norm3d_cel(i, j);
        //pt_n3d = norm3d + ij;
        //calcule_norm3d_cel_iso(i, j);
        calcule_norm3d_cel_pente_max_interp(i, j, r);
        //pt_n3d = norm3d_iso + ij;
        pt_n3d = norm3d + ij;
        ny = pt_n3d->y;
        //ny = n3d.y;
        //if (0){
        if (ny < ny_lim){
          pt = alti+ij;
          nx = pt_n3d->x; //composante EO
          nz = pt_n3d->z;  //composante NS
          //calcul de l'angle au sol, pour determiner la direction de la plus forte pente
          if (nx){
            anh = atan(-nz/nx);
            if (nx<0) anh += PI;
            if (anh<0) anh += 2*PI;
          }
          else{
            anh = (nz>0) ? 3*PI/2 : PI/2;
          }
          step_alti = 0;
#if 1 // 8 directions
          if ((anh<PI/8) || (anh>15*PI/8)){
            if (i<LEO-1) {step_alti = 1; step = 1;}  //droite
          }
          else if (anh<3*PI/8){
            if ((i<LEO-1) && (j>jmin)) {step_alti = -LEO+1; step = -HL+1;}  //diagonale droite+derriere
          }
          else if (anh<5*PI/8){
            if (j>jmin) {step_alti = -LEO; step = -HL;}  //derriere
          }
          else if (anh<7*PI/8){
            if ((i>0) && (j>jmin)) {step_alti = -LEO-1; step = -HL-1;}  //diagonale gauche+derriere
          }
          else if (anh<9*PI/8){
            if (i>0) {step_alti = -1; step = -1;}  //gauche
          }
          else if (anh<11*PI/8){
            if ((i>0) && (j<jmax)) {step_alti = LEO-1; step = HL-1;}  //diagonale gauche+devant
          }
          else if (anh<13*PI/8){
            if (j<jmax) {step_alti = LEO; step = HL;}  //devant
          }
          else if (anh<15*PI/8){
            if ((i<LEO-1) && (j<jmax)) {step_alti = LEO+1; step = HL+1;}  //diagonale droite+devant
          }
#else // 4 directions
          if ((anh<PI/4) || (anh>7*PI/4)){
            if (i<LEO-1) {step_alti = 1; step = 1;}  //droite
          }
          else if (anh<3*PI/4){
            if (j>jmin) {step_alti = -LEO; step = -HL;}  //derriere
          }
          else if (anh<5*PI/4){
            if (i>0) {step_alti = -1; step = -1;}  //gauche
          }
          else if (anh<7*PI/4){
            if (j<jmax) {step_alti = LEO; step = HL;}  //devant
          }
#endif
#ifdef CYCLAGE_HOR
        //TODO
#endif
          if (step_alti) {
            for (i_cel=0; i_cel<nb_cel_max; i_cel++){
              pt2 = pt + step_alti;
              //if ((pt2<alti) || (pt2>alti+(L-2)*LNS)){ErrPrintf("ERREUR : avalanches_norm, depassement de tableau alti\n");exit(-1);} //antibug
              if ((pt2<alti) || (pt2>alti+LEO*LNS)){ErrPrintf("WARNING: avalanches_norm - array overflow in alti\n");break;} //antibug
              k = H-(*pt);
              ix = 1+i+k*L+(LN+j)*HL; //cellule du sommet
              //if ((ix<0) || (ix >= HLD)) {ErrPrintf("ERREUR: avalanche - ix = %d\n", ix);  exit(-1);}
              aux = TE+ix;
              cel = *aux;
  //#if CELL_DATA
  //            get_cell_data(data, ix);
  //#endif
              if (aux->celltype == DUM) break;
              if (aux->celltype != typ) {
                ErrPrintf("WARNING: avalanches_norm - type incorrect : %d (i = %d  j = %d  k = %d  nb_ava_iter = %d)\n", aux->celltype, i, j, k, nb_ava_iter);
                //ErrPrintf("type  : %d (i = %d  j = %d  k = %d  nb_ava_iter = %d)\n", (aux-L)->celltype, i, j, k-1, nb_ava_iter);
                //ErrPrintf("type  : %d (i = %d  j = %d  k = %d  nb_ava_iter = %d)\n", (aux+L)->celltype, i, j, k+1, nb_ava_iter);
                //ErrPrintf("nb_ava = %d\n", nb_ava);
                break;
                //exit(-1);
              }
              if ((*pt)-(*pt2)>0){
              //if (1){
                if ((*pt)-(*pt2)>1){
                  flag_drop=1; //indicateur de chute verticale, pour forcer l'arret des avalanches en l'absence de chute
                  nb_ava_drop++;
                }
                ix2 = ix + step;  //cellule voisine
              //aux2 = aux + ix2;
              //on decale la cellule voisine vers le sommet (pour conservation de la matiere)
              //aux->celltype = aux2->celltype;
                //ajoute_cellule(TE[ix2].celltype, ix);
                deplace_cellule(ix2, ix);
              //on fait remonter la colonne de cellules
              //while (((aux2+L)->celltype != typ) && ((aux2+L)->celltype != DUM) && ((aux2+L)->celltype != BORD)){
                for (k=*pt; k>(*pt2)+1; k--){
                //aux2->celltype = (aux2+L)->celltype;
                  ix = ix2;
                  ix2+=L;
                  //ajoute_cellule(TE[ix2].celltype, ix);
                  deplace_cellule(ix2, ix);
                }
              }
              else{
                //cas ou la cellule "haute" est plus basse que la cellule "basse" !!!
                //on echange simplement les cellules "haute" et "basse"
                //ix2 = ix + step + ((*pt)-(*pt2)-1)*L; //cellule "basse"
                //LogPrintf("cellule basse : %d   passe : %d   k1 = %d   k2 = %d\n", TE[ix2].celltype, nb_ava_iter, (*pt), (*pt2));
                //ajoute_cellule(TE[ix2].celltype, ix);
                continue;
              }
              // on fait tomber la cellule haute dans la direction choisie
              //aux2->celltype = typ;
              //ajoute_cellule(typ, ix2);
              init_cellule(cel, ix2);
              //mise-a-jour des elevations locales
              (*pt)--;
              (*pt2)++;
              if (alti_mode == ALTI_MODE_HAUT){ //recalcul de l'altitude
                aux += L;
                //while ((aux->celltype != typ) && (aux->celltype != DUM) && (aux->celltype != BORD)) {(*pt)--; aux += L;}
                while (Phase[aux->celltype] != SOLID) {(*pt)--; aux += L;}
              }
              flag_ava = 1;
              nb_ava++;
              //LogPrintf("anh = %f   step_alti = %d   i = %d   j = %d\n", anh, step_alti, i, j);
            }
          }
        }
      }
    }
    //fprintf(fp, "avalanches_norm : nb_ava = %d\n", nb_ava);
  //}while ((flag_ava) && ((flag_drop) || (nb_ava_iter < 10)));
  }while (flag_drop);
  //}while ((flag_ava) && (nb_ava_iter < 100));
  //}while (flag_ava);

  //dump avalanches
  //fprintf(fp,"avalanches_norm : %d   passes : %d   temps : %f   angle : %f\n", nb_ava, nb_ava_iter, temps, angle*180/PI); //fflush(stdout);
  fprintf(fp,"avalanches_norm : %d (mov) %d (drop)   passes : %d   csp_time : %f   angle : %f\n", nb_ava, nb_ava_drop, nb_ava_iter, csp_time, angle); //fflush(stdout);
  fclose(fp);

  if (nb_ava){
    //recalcul des tableaux de doublets actifs
    init_db_pos();
    //recalcul des cellules
    init_Ncel();

#ifdef PARALLEL
  //synchronisation des zones de recouvrement
    if (mode_par) synchro_tunnel_zones(1);
#endif
  }

}


// avalanches localisees avec propagation vers les 4 premiers voisins (fonction recursive)
// TODO : mode CELL_COLOR
// TODO : mode CYCLAGE_HOR
void ava_propag(int i, int j)
{
  short *pt, *pt2;
  unsigned char flag_ava, flag_drop;
  int step, step_alti;
  int jmin, jmax, ix, ix2, ij;
  short h, k;

  if (pbc_mode && rot_map){
    ErrPrintf("ERROR: sorry, rotating space is not compatible with synchronous avalanches\n");
    exit(-1);
  }

  flag_ava = flag_drop = 0;
  jmin = 0;
  jmax = LNS-1;
  step_alti = 0;

  //Verification des differences d'altitude avec les 4 premiers voisins
  ij = i+j*LEO;
  ava_mask[ij] = 0;

  pt = alti+ij;
  if (i+1<LEO){
    h = (*pt) - (*(pt+1));
    if (h>ava_h_lim) {step_alti = 1; step = 1;}  //droite
  }
  if (!step_alti && (i>0)){
    h = (*pt) - (*(pt-1));
    if (h>ava_h_lim) {step_alti = -1; step = -1;}  //gauche
  }
  if (!step_alti && (j+1<LNS)){
    h = (*pt) - (*(pt+LEO));
    if (h>ava_h_lim) {step_alti = LEO; step = HL;}  //devant
  }
  if (!step_alti && (j>0)){
    h = (*pt) - (*(pt-LEO));
    if (h>ava_h_lim) {step_alti = -LEO; step = -HL;}  //derriere
  }
#ifdef CYCLAGE_HOR
  if (pbc_mode){
    if (!step_alti && (i==0)){
      h = (*pt) - (*(pt+LEO-1));
      if (h>ava_h_lim) {step_alti = LEO-1; step = L-3;}  //gauche
    }
    if (!step_alti && (i==LEO-1)){
      h = (*pt) - (*(pt-LEO+1));
      if (h>ava_h_lim) {step_alti = -(LEO-1); step = -(L-3);}  //droite
    }
    if (!step_alti && (j==0)){
      h = (*pt) - (*(pt+(LNS-1)*LEO));
      if (h>ava_h_lim) {step_alti = (LNS-1)*LEO; step = (L-3)*HL;}  //derriere
    }
    if (!step_alti && (j==LNS-1)){
      h = (*pt) - (*(pt-(LNS-1)*LEO));
      if (h>ava_h_lim) {step_alti = -(LNS-1)*LEO; step = -(L-3)*HL;}  //devant
    }
  }
#endif

  //if (step_alti) {LogPrintf("Chute : h=%d   step_alti=%d\n", h, step_alti);}

  if (!step_alti){
    //Chute dans la direction de plus grande pente (4 ou 8 directions possibles)
    float r, nx, ny, nz, ny_lim, anh;
    Vec3 *pt_n3d;

    ny_lim = cos(PI*ava_propag_angle/180);
    r = (ava_propag_angle <= 33.0) ? 3.5 : 3.0;

    calcule_norm3d_cel_pente_max_interp(i, j, r);
    pt_n3d = norm3d + ij;
    ny = pt_n3d->y;
    //if (0){
    if (ny < ny_lim){
      nx = pt_n3d->x; //composante EO
      nz = pt_n3d->z;  //composante NS
      //calcul de l'angle au sol, pour determiner la direction de la plus forte pente
      if (nx){
        anh = atan(-nz/nx);
        if (nx<0) anh += PI;
        if (anh<0) anh += 2*PI;
      }
      else{
        anh = (nz>0) ? 3*PI/2 : PI/2;
      }
  #if 1 // 8 directions
      if ((anh<PI/8) || (anh>15*PI/8)){
        if (i<LEO-1) {step_alti = 1; step = 1;}  //droite
      }
      else if (anh<3*PI/8){
        if ((i<LEO-1) && (j>jmin)) {step_alti = -LEO+1; step = -HL+1;}  //diagonale droite+derriere
      }
      else if (anh<5*PI/8){
        if (j>jmin) {step_alti = -LEO; step = -HL;}  //derriere
      }
      else if (anh<7*PI/8){
        if ((i>0) && (j>jmin)) {step_alti = -LEO-1; step = -HL-1;}  //diagonale gauche+derriere
      }
      else if (anh<9*PI/8){
        if (i>0) {step_alti = -1; step = -1;}  //gauche
      }
      else if (anh<11*PI/8){
        if ((i>0) && (j<jmax)) {step_alti = LEO-1; step = HL-1;}  //diagonale gauche+devant
      }
      else if (anh<13*PI/8){
        if (j<jmax) {step_alti = LEO; step = HL;}  //devant
      }
      else if (anh<15*PI/8){
        if ((i<LEO-1) && (j<jmax)) {step_alti = LEO+1; step = HL+1;}  //diagonale droite+devant
      }
  #else // 4 directions
      if ((anh<PI/4) || (anh>7*PI/4)){
        if (i<LEO-1) {step_alti = 1; step = 1;}  //droite
      }
      else if (anh<3*PI/4){
        if (j>jmin) {step_alti = -LEO; step = -HL;}  //derriere
      }
      else if (anh<5*PI/4){
        if (i>0) {step_alti = -1; step = -1;}  //gauche
      }
      else if (anh<7*PI/4){
        if (j<jmax) {step_alti = LEO; step = HL;}  //devant
      }
#endif
    }
  }
#ifdef CYCLAGE_HOR
  //TODO
#endif

  if (step_alti) { //la cellule se deplace dans la direction step_alti
    Cell *aux, cel;

    pt2 = pt + step_alti;
    //if ((pt2<alti) || (pt2>alti+(L-2)*LNS)){ErrPrintf("ERREUR : avalanches_norm, depassement de tableau alti\n");exit(-1);} //antibug
    if ((pt2<alti) || (pt2>alti+LEO*LNS)){ErrPrintf("WARNING: ava_propag - array overflow in alti\n");return;} //antibug
    k = H-(*pt);
    ix = 1+i+k*L+(LN+j)*HL; //cellule du sommet
    //if ((ix<0) || (ix >= HLD)) {ErrPrintf("ERREUR: avalanche - ix = %d\n", ix);  exit(-1);}
    aux = TE+ix;
    cel = *aux;
//#if CELL_DATA
//            get_cell_data(data, ix);
//#endif
    if (aux->celltype == DUM) return;
    if (aux->celltype != ALTI) {
      ErrPrintf("WARNING: ava_propag - incorrect type %d (i = %d  j = %d  k = %d )\n", aux->celltype, i, j, k);
      //ErrPrintf("type  : %d (i = %d  j = %d  k = %d  nb_ava_iter = %d)\n", (aux-L)->celltype, i, j, k-1, nb_ava_iter);
      //ErrPrintf("type  : %d (i = %d  j = %d  k = %d  nb_ava_iter = %d)\n", (aux+L)->celltype, i, j, k+1, nb_ava_iter);
      //ErrPrintf("nb_ava = %d\n", nb_ava);
      return;
      //exit(-1);
    }
    if ((*pt)-(*pt2)>0){
    //if (1){
      if ((*pt)-(*pt2)>1){
        flag_drop=1; //indicateur de chute verticale, pour forcer l'arret des avalanches en l'absence de chute
        nb_drop_ava_propag++;
      }
      ix2 = ix + step;  //cellule voisine
    //aux2 = aux + ix2;
    //on decale la cellule voisine vers le sommet (pour conservation de la matiere)
    //aux->celltype = aux2->celltype;
      //ajoute_cellule(TE[ix2].celltype, ix);
      deplace_cellule(ix2, ix);
    //on fait remonter la colonne de cellules
    //while (((aux2+L)->celltype != typ) && ((aux2+L)->celltype != DUM) && ((aux2+L)->celltype != BORD)){
      for (k=*pt; k>(*pt2)+1; k--){
      //aux2->celltype = (aux2+L)->celltype;
        ix = ix2;
        ix2+=L;
        //ajoute_cellule(TE[ix2].celltype, ix);
        deplace_cellule(ix2, ix);
      }
    }
    else{
      //cas ou la cellule "haute" est plus basse que la cellule "basse" !!!
      return;
    }
    // on fait tomber la cellule haute dans la direction choisie
    //aux2->celltype = typ;
    //ajoute_cellule(typ, ix2);
    init_cellule(cel, ix2);
    //mise-a-jour des elevations locales
    (*pt)--;
    (*pt2)++;
    flag_ava = 1;
    nb_ava_propag++;
    //if (nb_ava_propag % 100 == 0) init_Ncel(); //pour ralentir !!
    //LogPrintf("anh = %f   step_alti = %d   i = %d   j = %d\n", anh, step_alti, i, j);
  }

  //if (flag_drop){
  if (flag_ava){ //Propagation des avalanches
    Pos2 p2;
    if ((i+1<LEO) && !ava_mask[ij+1]){ //propagation a droite
      p2.x = i+1;
      p2.y = j;
      fifo_put(ava_fifo, p2);
      ava_mask[ij+1] = 1;
    }
    if ((i>0) && !ava_mask[ij-1]){ //propagation a gauche
      p2.x = i-1;
      p2.y = j;
      fifo_put(ava_fifo, p2);
      ava_mask[ij-1] = 1;
    }
    if ((j+1<LNS) && !ava_mask[ij+LEO]){ //propagation devant
      p2.x = i;
      p2.y = j+1;
      fifo_put(ava_fifo, p2);
      ava_mask[ij+LEO] = 1;
    }
    if ((j>0) && !ava_mask[ij-LEO]){ //propagation derriere
      p2.x = i;
      p2.y = j-1;
      fifo_put(ava_fifo, p2);
      ava_mask[ij-LEO] = 1;
    }
    //if (i+1<LEO) ava_propag(i+1,j);
    //if (i>0) ava_propag(i-1,j);
    //if (j+1<LNS) ava_propag(i,j+1);
    //if (j>0) ava_propag(i,j-1);
  }
}


void loop_ava_propag(int i, int j)
{
  Pos2 p2;
  int nb_iter, nb_fifo_max;

  if (!ava_mask){
    AllocMemoryPrint("ava_mask", ava_mask, char, LEO*LNS);
    ResetMemory(ava_mask, char, LEO*LNS);
  }

  if (!ava_fifo){
    ava_fifo = fifo_create(LEO*LNS);
  }

  nb_ava_propag = nb_drop_ava_propag = 0;
  nb_iter = nb_fifo_max = 0;
  ava_propag_angle = ava_angle_stable;
  LogPrintf("ava_propag_angle = %f\n",ava_propag_angle);

  ava_propag(i ,j);
  nb_iter++;

  while (ava_fifo->nb)
  {
    if (ava_fifo->nb > nb_fifo_max) nb_fifo_max = ava_fifo->nb;
    p2 = fifo_get(ava_fifo);
    ava_propag(p2.x, p2.y);
    nb_iter++;
  }

  LogPrintf("nb_ava_propag = %d\n", nb_ava_propag);
  LogPrintf("nb_drop_ava_propag = %d\n", nb_drop_ava_propag);
  LogPrintf("nb_iter = %d\n", nb_iter);
  LogPrintf("nb_fifo_max = %d\n", nb_fifo_max);

  if (nb_ava_propag){
    //recalcul des tableaux de doublets actifs
    init_db_pos();
    //recalcul des cellules
    init_Ncel();
#ifdef PARALLEL
  //synchronisation des zones de recouvrement
    if (mode_par) synchro_tunnel_zones(1);
#endif
    //sleep(3);
  }
}
#endif //AVALANCHES

#ifdef USE_VEGETATION
int check_vegetation(int ix, void *data)
{
  static char start=1;
  int x, y, z;

  if (start){
    LogPrintf("check_vegetation : veg_h_max = %d\n", veg_h_max);
    start=0;
  }

  //calcul de la position (x,y,z)
  Calcule_xyz(ix, x, y, z);

  return (y>=H-1-veg_h_max);
}
#endif

////////////////////////////////////////////////////////////////////////////////////:

void dump_surface(char* name, int cpt, int unit)
{
  char filename[100];
  FILE *fp;
  int i,j,n;
/*
  //dump de la composante verticale des normales
  if (norm3d){
    sprintf(nom,"NORM%04d.data", cpt);

    fp = fopen(nom,"w");
    if ( ! fp ){
      ErrPrintf("erreur ouverture fichier vitesse\n");
      exit(-4);
    }

    Vec3 *pt_n = norm3d;
    for(j=0; j<LNS; j++){
      for(i=0; i<LEO; i++, pt_n++){
        fprintf(fp, "%f ", pt_n->y);
      }
      fprintf(fp, "\n");
    }

    fclose(fp);
  }
*/
  if (unit == UNIT_COMP)
    sprintf(filename,"%s%04d.data", name, cpt);
  else
    sprintf(filename,"%s%05d_t0.data", name, cpt);

  if (alti){
    fp = fopen(filename,"w");
    if ( ! fp ){
      ErrPrintf("ERROR: cannot open file %s\n", filename);
      exit(-4);
    }

    short *pt_al = alti;
    for(j=0; j<LNS; j++){
      n=0;
      for(i=0; i<LEO; i++, pt_al++){
        if (rot_map && OutOfSpace(1+i,LN+j)) continue;
        fprintf(fp, "%d ", *pt_al);
        n++;
      }
      if (n>0) fprintf(fp, "\n");
    }

    fclose(fp);
  }
}


#ifdef DUMP_SIGMA
//extern double temps;
// calcul de l'ecart-type du signal d'altitude
void dump_sigma_alti()
{
  static float *mean;
  static float *sigma;
  long sum;
  float diff, var;
  int i, j;
  FILE *fp;
  int i0 = 0; //L/10; //offset pour reduire l'artefact lie aux effets de bords (non-raccordement des oscillations)

  //LogPrintf("dump_sigma_alti\n");

  //if (!alti) return;
//#ifdef  AVALANCHES
#if 0
  extern int ava_h_lim; // hauteur limite avant avalanche
  extern int ava_nb_cel_max;  // nb de cellules qui tombent simultanement
  avalanches(GR, ava_h_lim, ava_nb_cel_max, ALTI_MODE_BAS);
#else
  //calcule_alti(GR, ALTI_MODE_BAS);
#endif

  if (!mean){
    AllocMemory(mean, float, LNS);
    AllocMemory(sigma, float, LNS);
    LogPrintf("dump_sigma_alti : i0=%d\n", i0);
  }

  //calcul de l'altitude moyenne dans chaque couloir
  for (j=0; j<LNS; j+=2){
    sum = 0;
    for (i=i0; i<LEO; i++){
      sum += Alti(i, j);
    }
    mean[j] = (float)sum/(LEO-i0);
  }
  fp = fopen("MEAN_ALTI.data", "a");
  fprintf(fp, "%f  ", csp_time);
  for (j=0; j<LNS; j+=2) fprintf(fp, "%f  ", mean[j]);
  fprintf(fp, "\n");
  fclose(fp);

  //calcul de l'ecart-type dans chaque couloir
  for (j=0; j<LNS; j+=2){
    var = 0.0;
    for (i=i0; i<LEO; i++){
      diff = Alti(i, j) - mean[j];
      var += diff*diff;
    }
    var /= (LEO-i0);
    sigma[j] = sqrt(var);
  }
  fp = fopen("SIGMA.data", "a");
  fprintf(fp, "%f  ", csp_time);
  for (j=0; j<LNS; j+=2) fprintf(fp, "%f  ", sigma[j]);
  fprintf(fp, "\n");
  fclose(fp);

//#ifdef STABILITY_ANALYSIS
#if defined(STABILITY_ANALYSIS) && defined(LGCA)
  dump_vel();
#endif
}
#endif //DUMP_SIGMA


#ifdef DUMP_AUTOCORREL
//extern double temps;
// calcul de l'auto-correlation spatiale du signal d'altitude
void dump_autocorrel()
{
  static float *mean;
  static float *ac;
  static char nom[32];
  static int cpt=0;
  short *pt_al;
  float *pt_ac;
  long sum;
  float diff, var;
  int i, j, k;
  FILE *fp;
  int i0 = 0; //L/10; //offset pour reduire l'artefact lie aux effets de bords (non-raccordement des oscillations)

  //LogPrintf("dump_autocorrel\n");

  //if (cpt!=10*((int)cpt/10)) {cpt++; return;}

  //if (!alti) return;
//#ifdef  AVALANCHES
#if 0
  extern int ava_h_lim; // hauteur limite avant avalanche
  extern int ava_nb_cel_max;  // nb de cellules qui tombent simultanement
  avalanches(GR, ava_h_lim, ava_nb_cel_max, ALTI_MODE_BAS);
#else
  //calcule_alti(GR, ALTI_MODE_BAS);
#endif

  if (!mean){
    AllocMemory(mean, float, LNS);
    AllocMemory(ac, float, LNS*LEO);
    LogPrintf("dump_autocorrel : i0=%d\n", i0);
  }

  //calcul de l'altitude moyenne dans chaque couloir
  for (j=0; j<LNS; j+=2){
    sum = 0;
    for (i=i0; i<LEO; i++){
      sum += Alti(i, j);
    }
    mean[j] = (float)sum/(LEO-i0);
  }
  fp = fopen("MEAN_ALTI.data", "a");
  fprintf(fp, "%03d   %f  ", cpt, csp_time);
  for (j=0; j<LNS; j+=2) fprintf(fp, "%f  ", mean[j]);
  fprintf(fp, "\n");
  fclose(fp);

  //calcul de l'auto-correlation dans chaque couloir
  pt_ac = ac;
  for (j=0; j<LNS; j+=2){
    pt_al = alti+j*LEO; //altitudes dans le couloir j
    for (k=0; k<LEO/2; k++){
      sum = 0;
      for (i=i0; i<LEO; i++){
        if (i+k<LEO)
          sum += pt_al[i]*pt_al[i+k];
        else
          sum += pt_al[i]*pt_al[i+k-LEO];
      }
      *(pt_ac++) = (float)sum/LEO - mean[j]*mean[j];
    }
  }

  sprintf(nom,"AUTOCORREL%05d.data",cpt++);
  fp = fopen(nom, "w");
  //fprintf(fp, "%f  ", csp_time);
  pt_ac = ac;
  for (j=0; j<LNS; j+=2){
    for (k=0; k<LEO/2; k++){
      fprintf(fp, "%f  ", *(pt_ac++));
    }
    fprintf(fp, "\n");
  }
  fclose(fp);
}
#endif //DUMP_AUTOCORREL

#ifdef CGV
void dump_grad_vel(int cpt, int unit)
{
  static char filename[100];
  FILE *fp;
  float *pt_gv;
  int i, j, n;

  if (grdv){
    //calcule_grad_vel();

    if (unit == UNIT_COMP)
      sprintf(filename,"GRAD_VEL%04d.data", cpt);
    else
      sprintf(filename,"GRAD_VEL%05d_t0.data", cpt);

    fp = fopen(filename,"w");
    if (!fp){
      ErrPrintf("erreur ouverture fichier vitesse\n");
      exit(-4);
    }

    pt_gv=grdv;
    for(j=0; j<LNS; j++){
      n=0;
      for(i=0; i<LEO; i++, pt_gv++){
        //if (rot_map && OutOfSpace(1+i,LN+j)) continue;
        fprintf(fp, "%f ", *pt_gv);
        n++;
      }
      if (n>0) fprintf(fp, "\n");
    }

    fclose(fp);
  }
}

void dump_cgv()
{
  static float grdv_max = VSTEP_H*VSTEP_L*VSTEP_TIME*2;
  FILE *fp;
  float gv, prob;

  fp = fopen("PROB_CGV.data","w");
  if ( ! fp ){
    ErrPrintf("Erreur ouverture fichier PROB_CGV.data\n");
    exit(-4);
  }

  for(gv=Min(0,grdvc_min); gv<grdv_max; gv++){
    prob=prob_cgv(gv, grdvc_min, grdvc_max);
    fprintf(fp, "%f   %f\n", gv, prob);
  }

  fclose(fp);
}


void dump_cgv_coef()
{
  static int start=1, cpt=0;
  FILE *fp;

  fp = fopen("CGV_COEF.log","a");
  if (start){
    start = 0;
    fprintf(fp,"        cgv_coef");
  }
  fprintf(fp,"\n%04d:   %.4f", cpt++, cgv_coef);
  fclose(fp);
}
#endif //CGV
#endif //ALTI
