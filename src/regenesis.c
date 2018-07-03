/* ReSCAL - regenesis
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

#define _MAIN_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <libgen.h>
#include <stdint.h>

#include <math.h>
#include "defs.h"
#include "macros.h"
#include "cells.h"
#include "format.h"
//#include "regenesis.h"


#ifndef MOINS
#define MOINS 2
#endif

#define ResetCell(adrcell,type) { \
  memset(adrcell, 0, sizeof(Cell)); \
  (adrcell)->celltype = type; \
  }

int32_t H0, L0, D0, HL0, HLD0;       // les dimensions de la terre initiale
int32_t H, L, D, HL, HLD;       // les dimensions de la terre finale
double csp_time = 0.0;
int32_t Hd2, Ld2, Pd2, Hd3, Hd4;
Cell  *TE = NULL;           // notre 'terre'
Cell  *TE0 = NULL;           // ancienne 'terre'
char *csp_filename = NULL; //nom du fichier CSP
char out_format = DUMP_CSP; //format en sortie
int32_t graine = 68374;
float alpha = 0.0;  // angle de rotation
char opt_prefix = 0;

double drand48();
void  srand48();


void read_terre() {
  int32_t sz;
  sz = read_csp_header(csp_filename);

  out_format = sz ? DUMP_CSP : DUMP_BIN;

  if (!sz) {
    //binary file
    H = H0;
    L = L0;
    D = D0;
  } else {
    //CSP file
    H0 = H;
    L0 = L;
    D0 = D;
  }

  HL0 = HL = H * L;
  HLD0 = HLD = HL * D;

  AllocMemoryPrint("TE", TE, Cell, HLD);
  ResetMemory(TE, Cell, HLD);

  csp_set_bounds(0, 0, 0);

  read_csp(csp_filename);
}

// on copie TE dans TE0
void cp_TE_TE0() {
  if (!TE0) {
    AllocMemoryPrint("TE0", TE0, Cell, HLD);
  } else if (HLD != HLD0) {
    ReallocMemory(TE0, Cell, HLD);
    H0 = H;
    L0 = L;
    D0 = D;
    HL0 = HL;
    HLD0 = HLD;
  }

  memcpy(TE0, TE, sizeof(Cell) * HLD); //TE0 <- TE
}

void change_cel()
// in : TE
// out : TE
{
//   int32_t i, j, k, n;
//   float di, dj, dk;
//   int32_t Ld2, Pd2, Hd2;
//   float Ldx, Ldy, Ldz;
//   Cell *aux;

//   Ld2 = (int) L / 2;
//   Pd2 = (int) D / 2;
//   Hd2 = (int) H / 2;
//   int32_t lc = L * 0.3; //Ld2; //L*0.75; //Ld2+40; //L/4; //Ld2; //largeur couloir
//   Ldx = Ld2 - 0.5; //Ld3;//Ld4;
//   Ldy = Pd2 - 0.5; //Ld4; //(int) L/2;
//   Ldz = Hd2;
//   int32_t rcyl = Ld2 - 5; //rayon cylindre

  LogPrintf("change_cel\n");

//   aux = TE;
  // TODO: Figure out what is going on here
//   for (k = 0, n = 1; k < D; k++, n = 1 - n) // profondeur
//     for (j = 0; j < H; j++) // hauteur
//       for (i = 0; i < L; i++, aux++) { //largeur
//         di = i - Ldx;
//         dj = j - Ldz;
//         dk = k - Ldy;
        //if ((n==0) || (j > H0) ||  (j==0) || (j-0.08*k < 0))
        //aux->celltype   = DUM;
        //else if ((n==1))
        //aux->celltype   = GR;
        //if (aux->celltype == GR) aux->celltype = DUM; //dune en beton
        //if (aux->celltype == GRJ) aux->celltype = EAUC; //dune en beton
        //if (aux->celltype == IN) aux->celltype = DUM; //dune en beton
        //if ((aux->celltype == EAUT) || (aux->celltype == EAUC)) aux->celltype = (drand48() < 0.08) ? EAUC : EAUT; //antidunes
        //if (aux->celltype == EAUT) aux->celltype = EAUC; //pas de turbulence
        //if ((i<40) && (j==H-8)) aux->celltype=EAUC; //lissage
        //if (k!=Ld2) aux->celltype = DUM; //couloir de largeur 1 -> modele 2d
        //if (aux->celltype == DUM) aux->celltype = EAUC; //pas de DUM
        //if (aux->celltype == IN) aux->celltype = (k==Ld2)? EAUC : DUM; //pas de IN
        //if (aux->celltype == IN) aux->celltype = DUM; //pas de IN
        //if ((j>=1) && (j<=H-8) && (aux->celltype == DUM)) aux->celltype = EAUC; //on enleve le couloir, mais on garde le sol et le plafond
        //if ((k<Ld2-lc/2) || (k>Ld2+lc/2)) aux->celltype = DUM; //couloir de largeur lc
        //if ((i==5) && (j==0) && (k==Ld2)) aux->celltype = IN; //source lineaire
        //if ((i==5) && (j==0) && (k>=Ld2-60) && (k<=Ld2+60)) aux->celltype = IN; //source lineaire dans le couloir
        //if ((i==5) && (j==0) && (k >= Ld2-lc) && (k <= Ld2+lc))//source lineaire dans le couloir
        //if ((i==5) && (j==0)) aux->celltype = IN;
        //if ((j==1) && (aux->celltype == IN)) aux->celltype = EAUC;
        //if ((i==0) && (aux->celltype == GR)) aux->celltype = EAUC; //pas de grain contre le bord, SVP ...
        //if ((j == 0) && ((di*di + dk*dk >= rcyl*rcyl) && ((di+1)*(di+1) + dk*dk < rcyl*rcyl))) ////source de grains en demi-cercle
        //  aux->celltype = IN;
        //if ((j ==H-1) && (i > Ld2) && (di*di + dk*dk >= rcyl*rcyl)) ////sortie de grains en cercle
        //aux->celltype = OUT;
        /*float Ldx1 = Ld2-0.5-(L/3)*cos(PI/4);
        float Ldy1 = Pd2-0.5-(L/3)*sin(PI/4);
        if ((j >= H-1) && ((i-Ldx1)*(i-Ldx1) + (k-Ldy1)*(k-Ldy1) <= 25*16)) ////source de grains en disque au sol
          aux->celltype = IN;*/
//       }
}

//#if defined(MODEL_DUN) || defined(MODEL_SNO)
void rotation()
// in : TE0
// out : TE
{
  int32_t i, j, k;
  int32_t di, dk, i0, k0;
  Cell *aux;

  //alpha = PI;
  float co = cos(alpha);
  float si = sin(alpha);

  LogPrintf("rotation de %f degres\n", alpha * 180 / PI);
  cp_TE_TE0();
  aux = TE;
  for (k = 0; k < D; k++) // profondeur
    for (j = 0; j < H; j++) // hauteur
      for (i = 0; i < L; i++, aux++) { //largeur
        di = i - Ld2;
        dk = k - Pd2;
        if (di * di + dk * dk < Ld2 * Ld2) {
          i0 = Ld2 + co * di - si * dk;
          k0 = Ld2 + si * di + co * dk;
          *aux = *(TE0 + i0 + j * L0 + k0 * HL0);
        } else {
          *aux   = *(TE0 + i + j * L0 + k * HL0);
        }
        //if (aux->celltype == DUM) aux->celltype = EAUC;
      }
}
//#endif

void rotation_cyclic()
// in : TE0
// out : TE
{
  int32_t i, j, k;
  int32_t i0, k0;
  float di, dk;
  float di0, dk0;
  Cell *aux, *aux0;

  //alpha = PI;
  float co = cos(alpha);
  float si = sin(alpha);
  di0 = Ld2 - 0.5;
  dk0 = Pd2 - 0.5;

  LogPrintf("rotation cyclique de %f degres\n", alpha * 180 / PI);
  cp_TE_TE0();
  aux0 = TE0;
  aux = TE;

  for (k = 0; k < D; k++) { // profondeur
    for (j = 0; j < H; j++) { // hauteur
      for (i = 0; i < L; i++, aux++, aux0++) { //largeur
        di = i - di0;
        dk = k - dk0;
        if ((aux0->celltype != BORD) && (aux0->celltype != DUM) && (aux0->celltype != IN)) {
          i0 = roundf(di0 + co * di + si * dk);
          k0 = roundf(dk0 - si * di + co * dk);
          while (i0 < 0) {
            i0 += L;
          }
          while (i0 >= L) {
            i0 -= L;
          }
          while (k0 < 0) {
            k0 += D;
          }
          while (k0 >= D) {
            k0 -= D;
          }
          *aux = TE0[i0 + j * L + k0 * HL];

#if defined(MODEL_DUN) || defined(MODEL_SNO)
          if (aux->celltype == DUM) {
            aux->celltype = EAUC;
          }
#endif
        } else
          //aux2->celltype = aux->celltype;
        {
          *aux = *aux0;
        }
      }
    }
  }
}

void rotation90(int32_t n)
// in : TE0
// out : TE
{
  Cell *aux;
  int32_t i, j, k;
  int32_t i0, k0;
  float a = n * PI / 2;
  float co = cos(a);
  float si = sin(a);
  float di0 = Ld2 - 0.5;
  float dk0 = Ld2 - 0.5;
  float di, dk;

  LogPrintf("rotation complete de %f degres\n", a * 180 / PI);
  cp_TE_TE0();
  aux = TE;
  for (k = 0; k < D; k++) // profondeur
    for (j = 0; j < H; j++) // hauteur
      for (i = 0; i < L; i++, aux++) { //largeur
        di = i - di0;
        dk = k - dk0;
        i0 = di0 + co * (float)di - si * (float)dk;
        k0 = dk0 + si * (float)di + co * (float)dk;
        *aux = *(TE0 + i0 + j * L0 + k0 * HL0);
      }
}

void translation()
// in : TE
// out : TE
{
  int32_t i, j, k;
  Cell *aux;

  //translation est-ouest
  int32_t d_eo = L / 3; //L/3; //L/4; //L/10; //L/5;
  aux = TE;
  LogPrintf("translation : %d\n", d_eo);
  for (k = 0; k < D; k++) // profondeur
    for (j = 0; j < H; j++) // hauteur
      for (i = 0; i < L; i++, aux++) { //largeur
#if defined(MODEL_DUN) || defined(MODEL_SNO)
        //aux->celldata = (i<L-d_eo)? (aux+d_eo)->celldata : (aux->celldata == DUM)? DUM : EAUC;
        //if ((i>=L-d_eo) && (j==H-2)) aux->celldata = (drand48() < 0.2) ? GR : EAUC; //sable
        //if (i>=L-d_eo) aux->celldata = DUM;//(j<=H-3)? EAUC : DUM;
        if (i < L - d_eo) {
          if (aux->celltype != IN) {
            *aux = *(aux + d_eo);
          }
        } else {
          //if (aux->celldata != DUM && (aux+d_eo)->celldata != DUM) aux->celldata = (aux+d_eo-L)->celldata;//EAUC;
          //if (aux->celldata != DUM) aux->celldata = (aux+d_eo)->celldata;
          if (aux->celltype != DUM) {
            ResetCell(aux, EAUC);
          }
        }
#endif
      }
}

void translation_alldir()
// in : TE0
// out : TE
{
  int32_t i, j, k;
  Cell *aux;

  //translation est-ouest
  int32_t d_eo = 100; //L/3; //L/4; //L/10; //L/5;
  //translation nord-sud
  int32_t d_ns = 100;

  LogPrintf("translation : %d , %d\n", d_eo, d_ns);
  cp_TE_TE0();
  aux = TE;
  for (k = 0; k < D; k++) // profondeur
    for (j = 0; j < H; j++) // hauteur
      for (i = 0; i < L; i++, aux++) { //largeur
#if defined(MODEL_DUN) || defined(MODEL_SNO)
        if ((i + d_eo >= 0) && (i + d_eo < L) && (k + d_ns >= 0) && (k + d_ns < D)) {
          if (aux->celltype != IN) {
            *aux = *(TE0 + d_eo + i + j * L + (k + d_ns) * HL);
          }
        } else {
          if (aux->celltype != DUM) {
            ResetCell(aux, EAUC);
          }
        }
#endif
      }
}

void translation_cyclic()
// in : TE0
// out : TE
{
  int32_t i, j, k;
  Cell *aux;

  //translation est-ouest
  int32_t d_eo = 140; //L/3; //L/4; //L/10; //L/5;

  LogPrintf("translation : %d\n", d_eo);
  cp_TE_TE0();
  aux = TE;
  for (k = 0; k < D; k++) // profondeur
    for (j = 0; j < H; j++) // hauteur
      for (i = 0; i < L; i++, aux++) { //largeur
#if defined(MODEL_DUN) || defined(MODEL_SNO)
        if (i < L - d_eo) {
          if (aux->celltype != IN) {
            *aux = *(TE0 + d_eo + i + j * L + k * HL);
          }
        } else {
          if (aux->celltype != DUM) {
            *aux = *(TE0 + d_eo - L + i + j * L + k * HL);
          }
        }
#endif
      }
}

void translation_cyclic_ns()
// in : TE0
// out : TE
{
  int32_t i, j, k;
  Cell *aux;

  //translation nord_sud
  int32_t d_ns = 220; //L/3; //L/4; //L/10; //L/5;

  LogPrintf("translation : %d\n", d_ns);
  cp_TE_TE0();
  aux = TE;
  for (k = 0; k < D; k++) // profondeur
    for (j = 0; j < H; j++) // hauteur
      for (i = 0; i < L; i++, aux++) { //largeur
#if defined(MODEL_DUN) || defined(MODEL_SNO)
        if (k < D - d_ns) {
          if (aux->celltype != IN) {
            *aux = *(TE0 + i + j * L + (k + d_ns) * HL);
          }
        } else {
          if (aux->celltype != DUM) {
            *aux = *(TE0 + i + j * L + (d_ns - D + k) * HL);
          }
        }
#endif
      }
}

void transcale1()
// in : TE0
// out : TE
{
  int32_t i, j, k;
  Cell *aux;

  //translation + homothetie vers l'ouest

  int32_t d_eo = L / 10; //L/8;
  int32_t ii, jj, kk;
  int32_t d_ns = 0;
  float coef1 = 0.4;

  LogPrintf("translation 1 + rescaling\n");
  cp_TE_TE0();
  aux = TE;
  for (k = 0; k < D; k++) // profondeur
    for (j = 0; j < H; j++) // hauteur
      for (i = 0; i < L; i++, aux++) { //largeur
        ii = (i * coef1) + d_eo;
        jj = H - 1 - ((H - 1 - j) * coef1);
        kk = Pd2 + ((k - Pd2) * coef1) + d_ns;
#if defined(MODEL_DUN) || defined(MODEL_SNO)
        if (((ii > 0) && (ii < L) && (jj > 0) && (jj < H) && (kk > 0) && (kk < D))) {
          *aux = *(TE0 + ii + jj * L + kk * HL);
        } else {
          ResetCell(aux, EAUC);
        }
        //if ((j<H-1) && (aux->celltype == DUM)) ResetCell(aux,EAUC);
        if (j == 0) {
          ResetCell(aux, DUM);
        }
        //else if ( j == H-1) ResetCell(aux,EAUC);
#endif
      }
}

void transcale2()
// in : TE0
// out : TE
{
  int32_t i, j, k;
  Cell *aux;

  //translation + homothetie vers l'ouest
  //conservation de l'ancienne topographie

  int32_t d_eo = L / 8; //L/8;
  int32_t ii, jj, kk;
  //int32_t d_ns=2*D/15; //position 1
  //int32_t d_ns=1.7*D/15;
  int32_t d_ns = D / 15; //position 2
  //int32_t d_ns=0; //position 3
  float coef2 = 1.5;

  LogPrintf("translation 2 + rescaling\n");
  LogPrintf("d_ns = %d\n", d_ns);
  aux = TE;
  cp_TE_TE0();
  for (k = 0; k < D; k++) // profondeur
    for (j = 0; j < H; j++) // hauteur
      for (i = 0; i < L; i++, aux++) { //largeur
        ii = (i * coef2) + d_eo;
        jj = H - 1 - ((H - 1 - j) * coef2);
        kk = Pd2 + ((k - Pd2) * coef2) + d_ns;
#if defined(MODEL_DUN) || defined(MODEL_SNO)
        if (aux->celltype != GR) {
          if (((ii > 0) && (ii < L) && (jj > 0) && (kk > 0) && (kk < D))) {
            *aux = *(TE0 + ii + jj * L + kk * HL);
          } else {
            ResetCell(aux, EAUC);
          }
        }
        if ((j < H - 1) && (aux->celltype == DUM))
          ResetCell(aux, EAUC)
          else if (j == H - 1) {
            ResetCell(aux, DUM);
          }
#endif
      }
}

void rescale_height()
// in : TE0
// out : TE
{
  int32_t i, j, k, j0;
  Cell *aux;

  cp_TE_TE0();

  //redimensionnement
  H = 2 * H0;
  HL = H * L;
  HLD = HL * D;
  LogPrintf("dilatation verticale : %d\n", H);
  ReallocMemory(TE, Cell, HLD);
  aux = TE;

  for (k = 0; k < D; k++) { // profondeur
    for (j = 0; j < H; j++) { // hauteur
      for (i = 0; i < L; i++, aux++) { //largeur
        j0 = j * H0 / H;
        *aux = *(TE0 + i + j0 * L0 + k * HL0);
#if defined(MODEL_DUN) || defined(MODEL_SNO)
        if ((H < H0) && ((j == 0) || (j == H - 1))) {
          ResetCell(aux, DUM);
        }
#endif
      }
    }
  }
}

void change_height()
// in : TE0
// out : TE
{
  int32_t i, j, k;
  Cell *aux;

  cp_TE_TE0();

  //redimensionnement
  H = 200; //H0*2/3; //H0*2;
  Hd2 = H / 2;
  HL = H * L;
  HLD = HL * D;

  LogPrintf("redimensionnement vertical : %d x %d x %d \n", H, L, D);
  ReallocMemory(TE, Cell, HLD);
  aux = TE;

  for (k = 0; k < D; k++) // profondeur
    for (j = 0; j < H; j++) // hauteur
      for (i = 0; i < L; i++, aux++) { //largeur
#if defined(MODEL_DUN) || defined(MODEL_SNO)
        if (H0 - H + j > 0) {
          *aux = *(TE0 + i + (H0 - H + j) * L + k * HL0);
        } else if ((TE0 + i + k * HL0 + L)->celltype == DUM)
          ResetCell(aux, DUM)
          else {
            ResetCell(aux, EAUC);
          }
        if (j == 0) {
          ResetCell(aux, DUM);  //plafond
        }
#else
        if (H0 - H + j >= 0) {
          *aux = *(TE0 + i + (H0 - H + j) * L + k * HL0);
        } else {
          ResetCell(aux, MOINS);
        }
#endif
        //if (j==H-1) ResetCell(aux, DUM); //sol plat
      }

}


void add_layers()
// in : TE0
// out : TE
{
  int32_t i, j, k;
  int32_t h1, h2;
  Cell *aux;
  uint8_t typ;

  cp_TE_TE0();

  //redimensionnement : ajout de h1 couches en haut et h2 couches en bas
  h1 = 30; //en haut
  h2 = 10; //en bas
  H = H0 + h1 + h2;
  Hd2 = H / 2;
  HL = H * L;
  HLD = HL * D;

#if defined(MODEL_DUN) || defined(MODEL_SNO)
  typ = EAUC;
#else
  typ = 2;//MOINS;
#endif

  ReallocMemory(TE, Cell, HLD);
  aux = TE;

  printf("redimensionnement vertical : %d x %d x %d \n", H, L, D);

  for (k = 0; k < D; k++) // profondeur
    for (j = 0; j < H; j++) // hauteur
      for (i = 0; i < L; i++, aux++) { //largeur
        if ((j >= h1) && (j < H - h2)) {
          *aux = *(TE0 + i + (j - h1) * L + k * HL0);
        } else {
          ResetCell(aux, typ);
        }
      }
}

void change_width()
// in : TE0
// out : TE
{
  int32_t i, j, k;
  Cell *aux;

  cp_TE_TE0();

  //redimensionnement horizontal est-ouest
  int32_t d;
  L = L0 * 2; //L0*2;
  Ld2 = L / 2;
  HL = H * L;
  HLD = HL * D;
  d = (L - L0) / 2;

  ReallocMemory(TE, Cell, HLD);
  aux = TE;

  printf("redimensionnement horizontal est-ouest : %d x %d x %d \n", H, L, D);

  for (k = 0; k < D; k++) { // profondeur
    for (j = 0; j < H; j++) { // hauteur
      for (i = 0; i < L; i++, aux++) { //largeur
        if (i < d) {
          *aux = *(TE0 + j * L0 + k * HL0);
        } else if (i >= L - d) {
          *aux = *(TE0 + (L0 - 1) + j * L0 + k * HL0);
        } else {
          *aux = *(TE0 + (i - d) + j * L0 + k * HL0);
        }
      }
    }
  }
}


void change_depth()
// in : TE0
// out : TE
{
  int32_t i, j, k, d;
  Cell *aux;

  cp_TE_TE0();

  //redimensionnement horizontal nord-sud centre

  D = 20;
  d = (D - D0) / 2;
  HLD = HL * D;

  ReallocMemory(TE, Cell, HLD);
  aux = TE;

  printf("redimensionnement horizontal nord-sud : %d x %d x %d \n", H, L, D);

  for (k = 0; k < D; k++) { // profondeur
    for (j = 0; j < H; j++) { // hauteur
      for (i = 0; i < L; i++, aux++) { //largeur
        if (k < d) {
          *aux = *(TE0 + i + j * L0);
        } else if (k >= D - d) {
          *aux = *(TE0 + i + j * L0 + (D0 - 1) * HL0);
        } else {
          *aux = *(TE0 + i + j * L0 + (k - d) * HL0);
        }
      }
    }
  }
}

void regenerate_terre() {
  //LogPrintf("random seed: %d\n", graine);
  srand48(graine);

  Hd2 = (int) H0 / 2;
  Ld2 = (int) L0 / 2;
  Pd2 = (int) D0 / 2;
  Hd3 = (int) 2 * H0 / 3;
  Hd4 = (int) H0 / 4;

  // uncomment the functions that should be applied

  /***** Rotation d'angle alpha avec table tournante *****/
  if (alpha) {
    rotation();
  }

  /***** Rotation d'angle alpha avec bouclage *****/
  //if (alpha) rotation_cyclic();

  /***** Rotation complete d'angle n*90 degres *****/
  //rotation90(1);

  /***** Translation vers l'ouest *****/
  //translation();

  /***** Translation dans toute direction *****/
  //translation_alldir();

  /***** Translation est-ouest, avec bouclage (conservation des grains) *****/
  //translation_cyclic();

  /***** Translation nord-sud, avec bouclage (conservation des grains) *****/
  //translation_cyclic_ns();

  /***** Translation, homothetie vers l'ouest *****/
  //transcale1();

  /***** Translation, homothetie vers l'ouest et conservation de l'ancienne topo *****/
  //transcale2();

  /***** Dilatation verticale *****/
  //rescale_height();

  /***** Augmentation de la hauteur H *****/
  //change_height();

  /***** Augmentation de la longueur L *****/
  //change_width();

  /***** Augmentation de la profondeur D *****/
  //change_depth();

  /***** Ajout de couches horizontales en haut et en bas *****/
  //add_layers();

  /***** Modification du type des cellules *****/
  //change_cel();

  LogPrintf("size of cellular space: %d x %d x %d\n", H, L, D);
}


void dump_terre() {
  char out_filename[256];

  if (opt_prefix) {
    sprintf(out_filename, "Regen_%s", basename(csp_filename));
  } else {
    strcpy(out_filename, "Regen.csp");
  }

  csp_set_bounds(0, 0, 0);

  write_csp(out_format, out_filename);
}

void usage() {

  printf("usage: \n regenesis <CSP file> [OPTIONS]\n");
  printf(" regenesis <binary file> H L D [OPTIONS]\n");

  printf("  H \theight\n");
  printf("  L \tlenght\n");
  printf("  D \tdepth\n");

  printf("OPTIONS\n");
  printf("  -rot <ang>\trotation of angle <ang> (degrees)\n");
  printf("  -s <n> \trandom seed <n>\n");
  printf("  -prefix \tadd prefix \"Regen_\"\n");
  fflush(stdout);
  exit(-1);
}

void read_args(int32_t argc, char **argv) {
  int32_t n;

  n = 1;
  while ((n < argc) && (*argv[n] != '-')) {
    n++;
  }

  if ((n != 2) && (n != 5)) {
    usage();
  }

  csp_filename = argv[1];

  if (n == 5) {
    H0 = atoi(argv[2]);
    L0 = atoi(argv[3]);
    D0 = atoi(argv[4]);
  }

  for (; n < argc; n++) {
    if (!strcmp(argv[n], "-s") || !strcmp(argv[n], "-g")) {
      graine = atoi(argv[++n]);
    } else if (!strcmp(argv[n], "-rot")) {
      alpha = atof(argv[++n]);
      alpha = (alpha / 180) * PI;
    } else if (!strcmp(argv[n], "-prefix")) {
      opt_prefix = 1;
    } else {
      ErrPrintf("Error: option %s unknown\n", argv[n]);
      exit(-1);
    }
  }
}

int32_t main(int32_t argc, char **argv) {
  LogPrintf("regenesis %s\n", VER_NUM);

  read_args(argc, argv);

  read_terre();

  regenerate_terre();

  dump_terre();

  LogPrintf("regenesis : fin\n");

  return 1;
}

