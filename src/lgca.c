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


#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <assert.h>
#include <math.h>
#include <omp.h>
#include <stdint.h>

#ifdef __AVX__
#include <stdio.h>
#include <immintrin.h>
#endif

#include "defs.h"
#include "macros.h"
#include "param.h"
#include "cells.h"
#include "space.h"
#include "surface.h"
#include "transitions.h"
#include "lgca.h"
#include "trace.h"
#include "simul.h"

//#undef CYCLAGE_HOR	//pas de cyclage sur le flux (optionnel) -> sinon essoufflement du flux (frottements)

#ifdef LGCA

//#ifdef MVT_VFOR
#if defined(CYCLAGE_HOR) || defined(MVT_VFOR) || defined(MVT_GFOR)
//#define CYCLAGE_HOR //mode MVT_VFOR => cyclage horizontal sur le flux
#define CYCLAGE_MVT
#endif

extern int32_t pbc_mode;  //periodic boundary conditions
extern int32_t boundary;  //boundary conditions
extern int8_t *rot_map;        // periodic mapping of the rotating space
extern Pos2 *rot_map_pos;    // periodic mapping of the rotating space

#ifdef MVT_VFOR
float intensite_vfor = 0.5; //instensite du forcage vertical
#endif

#ifdef MVT_GFOR
float intensite_gfor = 0.001; //instensite du forcage global
#endif


extern uint8_t opt_info;
extern int32_t H, L, D, HL, HLD;       // les dimensions de la terre
extern int32_t LN, LS, LNS, HLN;    //couloir est-ouest (limite nord, limite sud, largeur nord-sud ...)
extern Cell  *TE;	           // la 'terre'
extern int32_t Ncel[];                       // nombre de cellules par type
extern uint8_t opt_vel;       //affichage des vitesses moyennes
extern Vec2 *norm2d;              //normale a la surface
//extern float *norm2dx;          // normale a la surface  (composante X)
//extern float *norm2dy;          // normale a la surface  (composante Y)
extern int32_t h_ceil;  //epaisseur du plafond
extern int8_t *psol;                // indicateur des plans solides verticaux est-ouest
extern const uint8_t Phase[MAX_CELL];	//phase (fluide ou solide) des types de cellules

MvtField *CelMvt;	//cellules en mouvements
MvtField *CelMvt2;	//cellules en mouvements
MvtField Collisions[SIZE_MVT_FIELD];  //tableau des collisions fuilde-fluide
#ifdef SOLSLOW
MvtField Collisions_SolSlow[16][SIZE_MVT_FIELD]; //tableau des collisions fluide-solide (cas d'une unique particule lente)
#endif
#ifdef SOLFAST
MvtField Collisions_SolFast[16][SIZE_MVT_FIELD]; //tableau des collisions fluide-solide (cas d'une unique particule rapide)
#endif
int32_t CelVelx[SIZE_MVT_FIELD];  //precomputed sums of horizontal velocities for a single node
int32_t CelVely[SIZE_MVT_FIELD];  //precomputed sums of vertical velocities for a single node

int32_t *Velx_time=NULL, *Vely_time=NULL;	//vitesses locales instantanees
int32_t *Velx=NULL, *Vely=NULL;	//vitesses locales moyennees (espace et temps)
int32_t *Veln=NULL;  //sliding window for the number of fluid nodes
float *Velx_interp=NULL, *Vely_interp=NULL;	//vitesses locales interpolees (espace) et moyennees (temps)
int32_t *Velx_vector=NULL, *Vely_vector=NULL; //sliding vectors of velocities
int32_t *Veln_vector=NULL; //sliding vector for the number of fluid nodes
int32_t CH, CLEO, CLNS, CHL, CHLD, CLN, CLS;     // dimensions du reseau de collisions (CLEO/CLNS = largeur du reseau de collisions en est-ouest/nord-sud)
int32_t VH, VL, VHL, VHLD;	// dimensions du tableau des vitesses
int32_t HLN;    //offset dans le reseau classique
int64_t nb_cel_fluide=0;	//nombre de cellules fluides
int64_t *NbMvt=NULL;		//nombre de particules en mouvement dans chaque couche du reseau
int64_t nb_mvt_in=0;		//nombre de particules entrant dans le reseau
int64_t nb_mvt_out=0;		//nombre de particules sortant du reseau
#ifdef MVT_REGUL
int64_t *NbMvt0=NULL;		//nombre de particules en mouvement au depart dans chaque couche du reseau
#endif
uint8_t NbMvtCell[SIZE_MVT_FIELD];  //nombre de particules en mouvement dans chaque etat d'une cellule
float densite=0;	//densite moyenne de particules en mouvement
float maxvel = 0.0;  //vitesse interpolee maximale apres stabilisation
float meanvel = 0.0;  //vitesse interpolee moyenne apres stabilisation
int32_t col_iter=0;   //nombre de cycles de collisions
int8_t *nom_fic_mvt = NULL; //nom du fichier VELOC (initial)
int32_t lgca_reset = 1;
int32_t lgca_speedup = 0; //initial speedup of the stabilization of lattice gas
int32_t lgca_mixing = 1; //mixing of the flux (on the left boundary)

void lecture_mvt();
void do_mvt_in(int32_t ix);
void do_regul_mvt();
void do_force_mvt();
void do_mixing();
void compute_vel_interp();
//void do_dump_mvt(int32_t inter);
void dump_densite();
void dump_mvt_in_out();
void dump_collisions();


void params_collisions()
{
  param_family("LGCA", "Lattice-gas parameters");

  /// MVT file
  parameter("Hpp_file", "initial lattice gas - optional", &nom_fic_mvt, PARAM_STRING, "LGCA");
  parameter_hidden("Mvt_file", "Mvt_file parameter is obsolete - Please use Hpp_file parameter instead !", &nom_fic_mvt, PARAM_STRING, "LGCA");

#ifdef MVT_GFOR
  parameter("Lgca_gfor", "global forcing of the flow - optional", &intensite_gfor, PARAM_FLOAT, "LGCA");
  parameter_hidden("Mvt_gfor", "Mvt_gfor parameter is obsolete - Please use Lgca_gfor parameter instead !", &intensite_gfor, PARAM_FLOAT, "LGCA");
#endif
  parameter("Lgca_reset", "reinitialize the flow after each rotation or translation (YES|NO) - optional", &lgca_reset, PARAM_BOOLEAN, "LGCA");
  parameter("Lgca_speedup", "initial speedup of the stabilization of lattice gas - optional", &lgca_speedup, PARAM_INT, "LGCA");
}

void init_collisions()
{
  int32_t i;

  /// size of the lattice
  CH=H*NB_MVT_VER;
  CLEO=L*NB_MVT_EO;
  CLNS=(LNS+DIST_MVT_NS-1)/DIST_MVT_NS;
  //CLNS_MAX = (rot_map)? D-2 : CLNS;
  CHL=CH*CLEO;

  /// largeur du couloir
  //for(i=0; (i<L) && plan_solide(i); i++);
  //CLN=i;
  //for(; (i<L) && !plan_solide(i); i++);
  //CLS=i-1;
  //CLNS=CLS-CLN+1;
  CLN = 0;
  CLS = CLNS-1;
  CHLD=CH*CLEO*CLNS;

  LogPrintf("init_collisions\n");
  LogPrintf("NB_MVT_VER = %d\n", NB_MVT_VER);
  LogPrintf("NB_MVT_EO = %d\n", NB_MVT_EO);
  LogPrintf("DIST_MVT_NS = %d\n", DIST_MVT_NS);
  LogPrintf("CLN = %d\n", CLN);
  LogPrintf("CLS = %d\n", CLS);
  LogPrintf("CH = %d\n", CH);
  LogPrintf("CLEO = %d\n", CLEO);
  LogPrintf("CLNS = %d\n", CLNS);
  LogPrintf("VSTEP_H = %d\n", VSTEP_H);
  LogPrintf("VSTEP_L = %d\n", VSTEP_L);
  LogPrintf("VSTEP_TIME = %d\n", VSTEP_TIME);
  VH = (CH + VSTEP_H - 1)/VSTEP_H;
  VL = (CLEO + VSTEP_L - 1)/VSTEP_L;
  LogPrintf("VH = %d\n", VH);
  LogPrintf("VL = %d\n", VL);
  VHL = VH*VL;
  VHLD = VHL*CLNS;
#if (VSTEP == 1)
  LogPrintf("VMEAN_REF = %d\n", VMEAN_REF);
#endif
#ifdef MVT_REGUL
  LogPrintf("regulation flux\n");
#endif
#ifdef MVT_HFOR
  LogPrintf("forcage flux par tapis roulant\n");
#endif
#ifdef MVT_VFOR
  LogPrintf("forcage flux sur le bord gauche, intensite = %f\n", intensite_vfor);
#endif
#ifdef MVT_GFOR
  LogPrintf("forcage global du flux, intensite = %f\n", intensite_gfor);
#endif
  lgca_mixing = !pbc_mode; //mixing when cellular space is not periodic
  if (lgca_mixing) LogPrintf("brouillage flux sur le bord gauche\n");
#ifdef SOLSLOW
  LogPrintf("rebonds realistes SOLSLOW\n");
#endif
#ifdef SOLFAST
  LogPrintf("rebonds realistes SOLFAST\n");
#endif
#ifdef SF_REAL
  LogPrintf("rebonds realistes SF_REAL\n");
//  LogPrintf("epaisseur du plafond = %d\n", h_ceil);
#endif
#ifdef PLAF_REAL
  LogPrintf("rebonds realistes PLAF_REAL\n");
//  LogPrintf("epaisseur du plafond = %d\n", h_ceil);
#endif
#ifdef VEL_SLIDE
  LogPrintf("sliding window for the averaging of the velocities\n");
#endif

  //allocation tableau CelMvt[]
  AllocMemoryPrint("CelMvt", CelMvt, MvtField, CHLD);
  ResetMemory(CelMvt, MvtField, CHLD);

  //allocation tableau CelMvt2[]
  AllocMemoryPrint("CelMvt2", CelMvt2, MvtField, CHLD);
  ResetMemory(CelMvt2, MvtField, CHLD);

  //allocation tableau NbMvt[]
  AllocMemory(NbMvt, long, CLNS);
  ResetMemory(NbMvt, long, CLNS);

#ifdef MVT_REGUL
  AllocMemory(NbMvt0, long, CLNS);
  ResetMemory(NbMvt0, long, CLNS);
#endif
  PrintTotalMemory();

  if (nom_fic_mvt){
    lecture_mvt();
  }
  else{
    init_mvt();
  }

  // initialisation tableau NbMvtCell[]
  for (i=0; i<SIZE_MVT_FIELD; i++){
    NbMvtCell[i] = 0;
    if (i & MVT_E) NbMvtCell[i]++;
    if (i & MVT_O) NbMvtCell[i]++;
    if (i & MVT_B) NbMvtCell[i]++;
    if (i & MVT_H) NbMvtCell[i]++;
    if (i & MVT_EB) NbMvtCell[i]++;
    if (i & MVT_EH) NbMvtCell[i]++;
    if (i & MVT_OB) NbMvtCell[i]++;
    if (i & MVT_OH) NbMvtCell[i]++;
  }

  /// number of gas particules
#ifdef MVT_REGUL
  for(k=CLN; k<=CLS; k++){
    if (NbMvt[k]) //TODO
      NbMvt0[k] = NbMvt[k];
    else //approximation
#ifdef MODEL_DUN
      NbMvt0[k] = Ncel[EAUC]*DENSITE/CLNS;
#else
      NbMvt0[k] = Ncel[MOINS]*DENSITE/CLNS;
#endif
    LogPrintf("NbMvt0[%d] = %ld\n", k, NbMvt0[k]);
  }
#endif

  LogPrintf("nb_mvt_in = %ld\n", nb_mvt_in);

  // initialisation tableau Collisions[]
  for (i=0; i<SIZE_MVT_FIELD; i++){
    if (i & MVT_OUT){
      Collisions[i] = MVT_OUT;
    }
    else if (i & MVT_SOLID){
      Collisions[i] = MVT_SOLID;
      if (i & MVT_E) Collisions[i] |= MVT_O;
      //if (i & MVT_E) Collisions[i] |= MVT_H; //test
      if (i & MVT_O) Collisions[i] |= MVT_E;
      if (i & MVT_B) Collisions[i] |= MVT_H;
      if (i & MVT_H) Collisions[i] |= MVT_B;

      if (i & MVT_EB) Collisions[i] |= MVT_OH; //MVT_OH;
      if (i & MVT_EH) Collisions[i] |= MVT_OB; //MVT_OB;
      if (i & MVT_OB) Collisions[i] |= MVT_EH; //MVT_EH;
      if (i & MVT_OH) Collisions[i] |= MVT_EB; //MVT_EB;

      /*if (i & MVT_EB) Collisions[i] |= MVT_EH; //MVT_OH;
      if (i & MVT_EH) Collisions[i] |= MVT_EB; //MVT_OB;
      if (i & MVT_OB) Collisions[i] |= MVT_OH; //MVT_EH;
      if (i & MVT_OH) Collisions[i] |= MVT_OB; //MVT_EB;*/

      /*if (i & MVT_EB) Collisions[i] |= MVT_OB; //MVT_OH;
      if (i & MVT_EH) Collisions[i] |= MVT_EB; //MVT_OB;
      if (i & MVT_OB) Collisions[i] |= MVT_OH; //MVT_EH;
      if (i & MVT_OH) Collisions[i] |= MVT_EH; //MVT_EB;*/
    }
    else{
      Collisions[i] = 0;

      if ((i & MVT_SLOW) == MVT_E){
        if ((i & MVT_FAST) == MVT_OB) Collisions[i] |= (MVT_O | MVT_EB);
	else if ((i & MVT_FAST) == MVT_OH) Collisions[i] |= (MVT_O | MVT_EH);
	else Collisions[i] |= MVT_E;
      }
      else if ((i & MVT_SLOW) == MVT_O){
        if ((i & MVT_FAST) == MVT_EB) Collisions[i] |= (MVT_E | MVT_OB);
	else if ((i & MVT_FAST) == MVT_EH) Collisions[i] |= (MVT_E | MVT_OH);
	else Collisions[i] |= MVT_O;
      }
      else if ((i & MVT_SLOW) == MVT_B){
        if ((i & MVT_FAST) == MVT_EH) Collisions[i] |= (MVT_H | MVT_EB);
	else if ((i & MVT_FAST) == MVT_OH) Collisions[i] |= (MVT_H | MVT_OB);
	else Collisions[i] |= MVT_B;
      }
      else if ((i & MVT_SLOW) == MVT_H){
        if ((i & MVT_FAST) == MVT_EB) Collisions[i] |= (MVT_B | MVT_EH);
	else if ((i & MVT_FAST) == MVT_OB) Collisions[i] |= (MVT_B | MVT_OH);
	else Collisions[i] |= MVT_H;
      }
      else if ((i & MVT_SLOW) == (MVT_E | MVT_O)) Collisions[i] |= (MVT_H | MVT_B);
      else if ((i & MVT_SLOW) == (MVT_H | MVT_B)) Collisions[i] |= (MVT_E | MVT_O);
      else Collisions[i] |= (i & MVT_SLOW);

      if ((i & MVT_FAST) == (MVT_EB | MVT_OH)) Collisions[i] |= (MVT_EH | MVT_OB);
      else if ((i & MVT_FAST) == (MVT_EH | MVT_OB)) Collisions[i] |= (MVT_EB | MVT_OH);
      else if (!(Collisions[i] & MVT_FAST)) Collisions[i] |= (i & MVT_FAST);
    }
  }

#ifdef SOLSLOW
  // initialisation (d'une partie) du tableau Collisions_SolSlow[]
  int32_t mvt_sol_E, mvt_sol_O, mvt_sol_H, mvt_sol_B, mvt_sol;
  LogPrintf("collisions realistes : solslow\n");
  memset(Collisions_SolSlow,0,sizeof(Collisions_SolSlow));
  for (mvt_sol=0; mvt_sol<16; mvt_sol++){
    mvt_sol_E = mvt_sol & (MVT_E >> 2);
    mvt_sol_O = mvt_sol & (MVT_O >> 2);
    mvt_sol_B = mvt_sol & (MVT_B >> 2);
    mvt_sol_H = mvt_sol & (MVT_H >> 2);

    Collisions_SolSlow[mvt_sol][MVT_E] = (!mvt_sol_H) ? MVT_H : (!mvt_sol_B) ? MVT_B : MVT_O;
    Collisions_SolSlow[mvt_sol][MVT_O] = (!mvt_sol_H) ? MVT_H : (!mvt_sol_B) ? MVT_B : MVT_E;
    Collisions_SolSlow[mvt_sol][MVT_B] = (!mvt_sol_E) ? MVT_E : (!mvt_sol_O) ? MVT_O : MVT_H;
    //Collisions_SolSlow[mvt_sol][MVT_B] = MVT_H;
    Collisions_SolSlow[mvt_sol][MVT_H] = (!mvt_sol_E) ? MVT_E : (!mvt_sol_O) ? MVT_O : MVT_B;
  }
#endif

#ifdef SOLFAST
  // initialisation (d'une partie) du tableau Collisions_SolFast[]
  int32_t mvt_sol_EB, mvt_sol_OB, mvt_sol_EH, mvt_sol_OH;//, mvt_sol;
  LogPrintf("collisions realistes : solfast\n");
  memset(Collisions_SolFast,0,sizeof(Collisions_SolFast));
  for (mvt_sol=0; mvt_sol<16; mvt_sol++){
    mvt_sol_EB = mvt_sol & (MVT_EB >> 6);
    mvt_sol_EH = mvt_sol & (MVT_EH >> 6);
    mvt_sol_OB = mvt_sol & (MVT_OB >> 6);
    mvt_sol_OH = mvt_sol & (MVT_OH >> 6);

    Collisions_SolFast[mvt_sol][MVT_EB] = (!mvt_sol_EH) ? MVT_EH : (!mvt_sol_OB) ? MVT_OB : MVT_OH;
    Collisions_SolFast[mvt_sol][MVT_EH] = (!mvt_sol_OH) ? MVT_OH : (!mvt_sol_EB) ? MVT_EB : MVT_OB;
    Collisions_SolFast[mvt_sol][MVT_OB] = (!mvt_sol_OH) ? MVT_OH : (!mvt_sol_EB) ? MVT_EB : MVT_EH;
    Collisions_SolFast[mvt_sol][MVT_OH] = (!mvt_sol_EH) ? MVT_EH : (!mvt_sol_OB) ? MVT_OB : MVT_EB;
    //Collisions_SolFast[mvt_sol][MVT_EB] = MVT_EH;
    //Collisions_SolFast[mvt_sol][MVT_EH] = MVT_EB;
    //Collisions_SolFast[mvt_sol][MVT_OB] = MVT_OH;
    //Collisions_SolFast[mvt_sol][MVT_OH] = MVT_OB;
  }
#endif

  if (opt_info) dump_collisions();

  // initialization of precomputed sums of velocities
  int32_t *pvx = CelVelx;
  int32_t *pvy = CelVely;
  for(i=0; i<SIZE_MVT_FIELD; i++, pvx++, pvy++){
#ifdef VEL_SLIDE
    if (i & MVT_SOLID) continue;
#endif
    if (i & MVT_SLOW){
      if (i & MVT_E) (*pvx)++;
      if (i & MVT_O) (*pvx)--;
      if (i & MVT_B) (*pvy)++;
      if (i & MVT_H) (*pvy)--;
    }
    if (i & MVT_FAST){
      if (i & MVT_EB) {(*pvx)++; (*pvy)++;}
      if (i & MVT_EH) {(*pvx)++; (*pvy)--;}
      if (i & MVT_OB) {(*pvx)--; (*pvy)++;}
      if (i & MVT_OH) {(*pvx)--; (*pvy)--;}
    }
  }

  //dump_densite();

#ifdef DUMP_SIGNATURE
  //empreinte gazeuse
  if (opt_info) dump_signature_mvt();
#endif

#ifdef OPENMP
#pragma omp parallel
  {
    int32_t tid = omp_get_thread_num();
    if (tid==0){
      int32_t nthreads = omp_get_num_threads();
      LogPrintf("number of OMP threads : %d\n", nthreads);
    }
  }
#endif
}


// initialize one site of lattice gas
void init_mvt_cell(MvtField *mvt)
{
  static MvtField mvt_slow[] = {MVT_E, MVT_O, MVT_B, MVT_H};
  static MvtField mvt_fast[] = {MVT_EB, MVT_EH, MVT_OB, MVT_OH};

  int32_t alea;
  //if (drand48() < 0.5)
  /*if (j>0.6*H)*/{
    //do_mvt_in(k*HL+j*L+i);
    /* *mvt = MVT_E;
    nb_mvt++;
    float aleat = drand48();
    if (aleat < 0.25) {*mvt |= MVT_EB; nb_mvt++;}
    else if (aleat < 0.5) {*mvt |= MVT_EH; nb_mvt++;}*/

    // 1 particule
    /*int32_t alea8=drand48()*8;
    *aux_mvt = (alea8<1) ? MVT_E : 0;
    // *aux_mvt = (alea8<4) ? mvt_slow[alea8] : (alea8<8) ? mvt_fast[alea8-4] : 0;
            if (*aux_mvt) nb_mvt_in++;*/

    // 1 particule lente
    alea=drand48()*4;
    if (alea<4){
      *mvt = mvt_slow[alea];
      nb_mvt_in++;
    }
    //*aux_mvt = MVT_E;
    //*aux_mvt = (alea4) ? MVT_E : 0;

    // 1 particule rapide
    //alea4=drand48()*4;
    //*aux_mvt |= mvt_fast[alea4];
    //nb_mvt_in +=2;

    // 1/2 particule rapide
    alea=drand48()*8;
    if (alea<4){ // 1/2 particule rapide
      *mvt |= mvt_fast[alea];
      nb_mvt_in++;
    }

    // 2 particules lentes
    /*int32_t alea2=drand48()*2;
    *aux_mvt = mvt_slow[alea2]; // 1 particule lente horizontale
    alea2=drand48()*2;
    *aux_mvt |= mvt_slow[2+alea2]; // 1 particule lente verticale
    nb_mvt_in +=2;*/

    // 2 particules rapides
    /*alea2=drand48()*2;
    *aux_mvt |= mvt_fast[alea2]; // 1 particule rapide vers le bas
    alea2=drand48()*2;
    *aux_mvt |= mvt_fast[2+alea2]; // 1 particule rapide vers le haut
    nb_mvt_in +=2;*/

    // 6 particules
    /*int32_t alea4=drand48()*4;
    *aux_mvt = MVT_SLOW ^ mvt_slow[alea4]; // 3 particules lentes
    alea4=drand48()*4;
    *aux_mvt |= MVT_FAST ^ mvt_fast[alea4]; // 3 particules rapides
    nb_mvt_in +=6;*/
  }
}

// initialisation du tableau MvtCel des particules de gaz sur reseau
void init_mvt()
{
  int32_t i, j, k, ic;
  Cell *aux;
  MvtField *aux_mvt;

  LogPrintf("reset lattice gas\n");

  if (rot_map){
    /// set the north-south boundaries (we have to keep an extra vertical plan on north and south sides for the interpolation)
    CLN = 0;
    for(k=0; psol[LN + k*DIST_MVT_NS]; k++) CLN = k;
    for(; !psol[LN + k*DIST_MVT_NS]; k++);
    CLS = k;
    assert(CLS < CLNS);
  }
  //LogPrintf("init_mvt: CLN=%d, CLS=%d\n", CLN, CLS);

  // initialisation tableau CelMvt[]
  nb_cel_fluide = 0;
#if defined(OPENMP) && !defined(DETERMINISTIC)
#pragma omp parallel for private(i,j,k,ic,aux,aux_mvt) reduction(+:nb_cel_fluide)
#endif
  for(k=CLN; k<=CLS; k++){ //profondeur
    aux = TE + HLN + k*HL*DIST_MVT_NS;
    aux_mvt = CelMvt + k*HL*NB_MVT_EO;
    for(j=0; j<H; j++){ //hauteur
      for(i=0; i<L; i++, aux++){ //largeur
        for(ic=0; ic<NB_MVT_EO; ic++, aux_mvt++){
          if (Phase[aux->celltype] == FLUID){
            *aux_mvt=0;
            nb_cel_fluide++;
            init_mvt_cell(aux_mvt);
          }
          else if (Phase[aux->celltype] == SOLID)
            *aux_mvt = MVT_SOLID;

          //if ((j==0) || (j==H-1)) //bords horizontaux
          //if ((j==1) || (j==H-2)) //bords horizontaux
          //  *aux_mvt = MVT_SOLID;
//#if defined(MVT_REGUL) && !defined(CYCLAGE_HOR)
#if !defined(CYCLAGE_MVT)
#ifdef MVT_REGUL
          else if ((i==0) || (i==L-1)) //bords est-ouest
            *aux_mvt = MVT_OUT;
#else
          else if (i==L-1) // bord est
            *aux_mvt = MVT_OUT;
#endif
#endif
        }
      }
    }
  }

  if (rot_map) out_of_space_mvt(1);

  //speed up the stabilization of the flow
  if (lgca_speedup){
    //int32_t nb_cyc_for=VSTEP_TIME*10;
    LogPrintf("speedup of the lattice gas: %d\n", lgca_speedup);
    for(i=0; i<lgca_speedup; i++) do_force_mvt();
  }

  if (Velx_interp){
    ResetMemory(Velx_interp, float, HLD);
    ResetMemory(Vely_interp, float, HLD);
  }
}

// set the lattice gas out of the rotating space, using periodic boundary conditions
void out_of_space_mvt(int32_t reset_mvt)
{
  int32_t i, j, k, ic, ck, typ;
  MvtField *aux_mvt;
  Pos2 cp;

  if (reset_mvt) LogPrintf("initialization of the lattice gas outside the rotating space\n");

#if defined(OPENMP) && !defined(DETERMINISTIC)
#pragma omp parallel for private(i,j,k,ck,ic,aux_mvt,cp,typ) reduction(+:nb_cel_fluide)
#endif
  for(ck=CLN; ck<=CLS; ck++){ //profondeur
    k = LN + ck*DIST_MVT_NS;
    for(j=0; j<H; j++){ //hauteur
      aux_mvt = CelMvt + j*CLEO*NB_MVT_VER + ck*CHL;
      for(i=0; i<L; i++){ //largeur
        if (OutOfSpace(i, k)){
          if (pbc_mode){
            cp = RotMapPos(i, k);
            typ = CellType(cp.x, j, cp.y);
            for(ic=0; ic<NB_MVT_EO; ic++, aux_mvt++){
              if (Phase[typ] == SOLID)
                SetMask(*aux_mvt, MVT_SOLID);
              else{
                UnsetMask(*aux_mvt, MVT_SOLID);
                if (reset_mvt){
                  init_mvt_cell(aux_mvt);
                  nb_cel_fluide++;
                }
              }
            }
          }
          /*else if (boundary == BC_CLOSE){
            for(ic=0; ic<NB_MVT_EO; ic++, aux_mvt++){
              SetMask(*aux_mvt, MVT_SOLID);
            }
          }*/
          else{
            // Open boundary conditions: by default we put solid ground and ceiling
            for(ic=0; ic<NB_MVT_EO; ic++, aux_mvt++){
              if ((j<=1) || (j>=H-2)){
                SetMask(*aux_mvt, MVT_SOLID);
              }
              else{
                UnsetMask(*aux_mvt, MVT_SOLID);
                if (reset_mvt){
                  init_mvt_cell(aux_mvt);
                  nb_cel_fluide++;
                }
              }
            }
          }
        }
        else aux_mvt += NB_MVT_EO;
      }
    }
  }
}

// modification globale de la terre
void collisions_mod_terre()
{
  int32_t i, j, k, ic;
  Cell *aux;
  MvtField *aux_mvt;

  nb_cel_fluide = 0;
  aux = TE+HLN;
  aux_mvt = CelMvt;
  for(k=0; k<CLNS; k++, aux+=HL*(DIST_MVT_NS-1)){ //profondeur
    for(j=0; j<H; j++){ //hauteur
      for(i=0; i<L; i++, aux++){ //largeur
        for(ic=0; ic<NB_MVT_EO; ic++, aux_mvt++){
          if (Phase[aux->celltype] == SOLID){
            SetMask(*aux_mvt, MVT_SOLID);
          }
          else{
            UnsetMask(*aux_mvt, MVT_SOLID);
            nb_cel_fluide++;
          }
        }
      }
    }
  }
}

void lecture_mvt()
{
  FILE *fp;
  int32_t i, ix;

  LogPrintf("lecture fichier HPP : %s\n",nom_fic_mvt);

  fp = fopen(nom_fic_mvt,"r");
  if ( ! fp ){
    ErrPrintf("ERROR: cannot open HPP file %s\n", nom_fic_mvt);
    exit(-4);
  }

  fread(CelMvt, sizeof(MvtField), CHLD, fp);

  fclose(fp);

  // initialisation nb_cel_fluide
  for (ix=0; ix<CHLD; ix++){
    if (!(CelMvt[ix] & MVT_SOLID)) nb_cel_fluide++;
  }
  /*nb_cel_fluide = 0;*/

  // moyennage de la vitesse initiale
  /*for (i=0; i<VSTEP_TIME; i++)*/ compute_vel(0);
}


int32_t mvt(int32_t index)
{
  //return (CelMvt[index] & (MVT_E | MVT_O | MVT_B | MVT_H));
  return (CelMvt[index] & MVT_ALLDIR);
  //return (NbMvtCell[CelMvt[index]] > densite);
}


int32_t mvt_solid(int32_t ii)
{
  return (CelMvt[ii] & MVT_SOLID);
}

int32_t mvt_bas(int32_t ii)
{
  return (CelMvt[ii] & (MVT_B | MVT_EB | MVT_OB));
}

int32_t mvt_haut(int32_t ii)
{
  return (CelMvt[ii] & (MVT_H | MVT_EH | MVT_OH));
}

int32_t mvt_est(int32_t ii)
{
  //return (CelMvt[ii] & (MVT_EB | MVT_EH));
  return (CelMvt[ii] & (MVT_E | MVT_EB | MVT_EH));
}

int32_t mvt_ouest(int32_t ii)
{
  //return (CelMvt[ii] & MVT_O);
  //return (CelMvt[ii] & (MVT_OB | MVT_OH));
  return (CelMvt[ii] & (MVT_O | MVT_OB | MVT_OH));
}

int32_t check_mvt(int32_t index, void *data)
{
  int32_t ii;
  Calcule_cix(index, ii);
  return (CelMvt[ii] & MVT_ALLDIR);
}

int32_t check_mvt_solid(int32_t index, void *data)
{
  int32_t ii;
  Calcule_cix(index, ii);
  return (CelMvt[ii] & MVT_SOLID);
}

int32_t check_mvt_bas(int32_t index, void *data)
{
  int32_t ii;
  Calcule_cix(index, ii);
  /*
  int32_t val = 0;
  if (CelMvt[index] & MVT_B) val++;
  if (CelMvt[index] & MVT_EB) val++;
  if (CelMvt[index] & MVT_OB) val++;
  if (CelMvt[index] & MVT_H) val--;
  if (CelMvt[index] & MVT_EH) val--;
  if (CelMvt[index] & MVT_OH) val--;

  return (val > 0);
  */
  return (CelMvt[ii] & (MVT_B | MVT_EB | MVT_OB));
}

int32_t check_mvt_haut(int32_t index, void *data)
{
  int32_t ii;
  Calcule_cix(index, ii);
  /*
  int32_t val = 0;
  if (CelMvt[ii] & MVT_H) val++;
  if (CelMvt[ii] & MVT_EH) val++;
  if (CelMvt[ii] & MVT_OH) val++;
  if (CelMvt[ii] & MVT_B) val--;
  if (CelMvt[ii] & MVT_EB) val--;
  if (CelMvt[ii] & MVT_OB) val--;

  return (val > 0);
  */
  return (CelMvt[ii] & (MVT_H | MVT_EH | MVT_OH));
}

int32_t check_mvt_est(int32_t index, void *data)
{
  int32_t ii;
  Calcule_cix(index, ii);
  /*
  int32_t val = 0;
  if (CelMvt[index] & MVT_E) val++;
  if (CelMvt[index] & MVT_EB) val++;
  if (CelMvt[index] & MVT_EH) val++;
  if (CelMvt[index] & MVT_O) val--;
  if (CelMvt[index] & MVT_OB) val--;
  if (CelMvt[index] & MVT_OH) val--;

  return (val > 0);
  */
  return (CelMvt[ii] & (MVT_E | MVT_EB | MVT_EH));
  //return (CelMvt[ii] & MVT_E);
}

int32_t check_mvt_ouest(int32_t index, void *data)
{
  int32_t ii;
  Calcule_cix(index, ii);
  /*int32_t val = 0;
  if (CelMvt[ii] & MVT_O) val++;
  if (CelMvt[ii] & MVT_OB) val++;
  if (CelMvt[ii] & MVT_OH) val++;
  if (CelMvt[ii] & MVT_E) val--;
  if (CelMvt[ii] & MVT_EB) val--;
  if (CelMvt[ii] & MVT_EH) val--;

  return (val > 0);*/
  return (CelMvt[ii] & (MVT_O | MVT_OB | MVT_OH));
}

int32_t check_mvt_est_et_bas(int32_t index, void *data)
{
  return check_mvt_est(index, NULL) && check_mvt_bas(index, NULL);
}

int32_t check_mvt_EB(int32_t index, void *data)
{
  int32_t ii;
  Calcule_cix(index, ii);
  //return (CelMvt[ii] & MVT_B | MVT_E | MVT_EB); //bug
  //return (CelMvt[ii] & (MVT_B | MVT_E | MVT_EB));
  return (CelMvt[ii-10*CLEO] & (MVT_B | MVT_E | MVT_EB)); //magouille
}

int32_t check_no_mvt_est_et_bas(int32_t index, void *data)
{
  return !check_mvt_est_et_bas(index, NULL);
}

int32_t check_mvt_mean_bas(int32_t index, void *data)
{
#if (VSTEP > 1)
  int32_t cx, cy, cz;
  Calcule_cxyz(index, cx, cy, cz);
  return (Vely[(cx/VSTEP_L)+(cy/VSTEP_H)*VL+cz*VHL] > 0);
#else
  return (Vely[index] >0);
#endif
}

int32_t check_mvt_mean_est(int32_t index, void *data)
{
#if (VSTEP > 1)
  int32_t cx, cy, cz;
  Calcule_cxyz(index, cx, cy, cz);
  return (Velx[(cx/VSTEP_L)+(cy/VSTEP_H)*VL+cz*VHL] > 0);
#else
  return (Velx[index] >0);
#endif
}

int32_t check_no_mvt_mean_bas(int32_t index, void *data)
{
  return !check_mvt_mean_bas(index, NULL);
}

void collisions_modcell(int32_t type, int32_t index)
{
  //int32_t ii=(index-HLN)*NB_MVT_EO;
  int32_t x, y, z, ii, ic, cz;
  Calcule_xyz(index, x, y, z);
  cz = (z-LN)/DIST_MVT_NS;
  if (z != cz*DIST_MVT_NS + LN){
    //ErrPrintf("collisions_modcell : mauvais plan ! (%d)\n", z);
    return;
  }
  ii = x*NB_MVT_EO + y*CLEO + cz*CHL;
  //Calcul_cix(index, ii);

  for(ic=0; ic<NB_MVT_EO; ic++, ii++){
    if (Phase[type] == SOLID){
      if (!(CelMvt[ii] & MVT_SOLID)) nb_cel_fluide--;
      SetMask(CelMvt[ii], MVT_SOLID);
      //if (mvt(index)) {LogPrintf("vert piege !!\n");}
    }
    else{
      if (CelMvt[ii] & MVT_SOLID) nb_cel_fluide++;
      UnsetMask(CelMvt[ii], MVT_SOLID);
    }
  }
}


void do_collisions()
{
  int32_t i, j, k, ii, ii2;

#ifdef SF_REAL
  //const float PI = 3.14159;
  // angle limite pour les rebonds = PI/8
  // si l'angle avec la normale < angle limite, alors direction inversee, sinon changement de direction a angle droit
  //const float TP8 = tan(PI/8.0); //~0.41
  const float CP8 = cos(PI/8.0);
  const float SP8 = sin(PI/8.0);
  int8_t sfmod = 0;
#endif

#ifdef OPENMP
#pragma omp parallel for private(i,ii,j,k)
#endif
//shared(CLNS, CLEO, CH, Collisions, CelMvt, CelMvt2, NbMvtCell, h_ceil)
  //LogPrintf("do_collisions\n");
  for(k=CLN; k<=CLS; k++){ // profondeur
    /* Obtain thread number */
    //int32_t tid = omp_get_thread_num();

    ii = CLEO*CH*k + CLEO;
    //LogPrintf("tid=%d k=%d ii=%d\n", tid, k, ii);
    //ii+=CLEO;
    for(j=1; j<CH-1; j++){ // hauteur
      ii += NB_MVT_EO;
      for(i=NB_MVT_EO; i<CLEO-NB_MVT_EO; i++, ii++){ //largeur
        CelMvt2[ii] = Collisions[CelMvt[ii]];
        //CelMvt2[ii] = CelMvt[ii];

#ifdef PLAF_REAL
        if ((j <= h_ceil) && (NbMvtCell[CelMvt[ii] & MVT_FAST]==1)){
          //rebond realiste au plafond
          MvtField mvt2 = 0;
          switch (CelMvt[ii] & MVT_FAST)
          {
            case MVT_EH:
              mvt2 = MVT_EB;
              break;
            case MVT_OH:
              mvt2 = MVT_OB;
              break;
          }
          if (mvt2){
            UnsetMask(CelMvt2[ii], MVT_FAST);
            SetMask(CelMvt2[ii], mvt2);
          }
        }
#endif //PLAF_REAL

  //int32_t mvt_sol_slow = 0;//debug
  //int32_t mvt_sol_fast = 0;//debug
  //int32_t mvt_sfr_slow = 0;//debug
  //int32_t mvt_sfr_fast = 0;//debug

#ifdef SOLSLOW
        if (mvt_solid(ii) && (NbMvtCell[CelMvt[ii] & MVT_SLOW]==1)){ //Collisions ameliorees (fuilde-solide)
          //EST
          ii2 = ii+1;
#ifdef CYCLAGE_HOR
          if (pbc_mode && (i==CLEO-NB_MVT_EO-1)) ii2 -= CLEO-2*NB_MVT_EO;
#endif
          //OUEST
          int32_t mvt_sol_E = mvt_solid(ii2) ? (MVT_E >> 2) : 0;
          ii2 = ii-1;
#ifdef CYCLAGE_HOR
          if (pbc_mode && (i==NB_MVT_EO)) ii2 += CLEO-2*NB_MVT_EO;
#endif
          int32_t mvt_sol_O = mvt_solid(ii-1) ? (MVT_O >> 2) : 0;
          //BAS
          int32_t mvt_sol_B = mvt_solid(ii+CLEO) ? (MVT_B >> 2) : 0;
          //HAUT
          int32_t mvt_sol_H = mvt_solid(ii-CLEO) ? (MVT_H >> 2) : 0;
          int32_t mvt_sol = mvt_sol_E | mvt_sol_O | mvt_sol_B | mvt_sol_H;
          UnsetMask(CelMvt2[ii], MVT_SLOW);
          SetMask(CelMvt2[ii], Collisions_SolSlow[mvt_sol][CelMvt[ii] & MVT_SLOW]);
          //mvt_sol_slow = Collisions_SolSlow[mvt_sol][CelMvt[ii] & MVT_SLOW];  //debug
        }
#endif //SOLSLOW

#ifdef SOLFAST
        //if ((CelMvt[ii] & MVT_SOLID) && (NbMvtCell[CelMvt[ii] & MVT_FAST]==1)){ //Collisions ameliorees (fuilde-solide)
        if (mvt_solid(ii) && (NbMvtCell[CelMvt[ii] & MVT_FAST]==1)){ //Collisions ameliorees (fuilde-solide)
        //if (mvt_solid(ii) && (NbMvtCell[CelMvt[ii] & MVT_FAST]>=1)){ //Collisions ameliorees (fuilde-solide)
          /*int32_t mvt_sol_E = CelMvt[ii+1] & MVT_SOLID;
          int32_t mvt_sol_O = CelMvt[ii-1] & MVT_SOLID;
          int32_t mvt_sol_B = CelMvt[ii+CLEO] & MVT_SOLID;
          int32_t mvt_sol_H = CelMvt[ii-CLEO] & MVT_SOLID;*/

          //EST-BAS
          ii2 = ii+CLEO+1;
#ifdef CYCLAGE_HOR
          if (pbc_mode && (i==CLEO-NB_MVT_EO-1)) ii2 -= CLEO-2*NB_MVT_EO;
#endif
          int32_t mvt_sol_EB = mvt_solid(ii2) ? (MVT_EB >> 6) : 0;  //CelMvt[ii+CLEO+1] & MVT_SOLID;

          //EST-HAUT
          ii2 = ii-CLEO+1;
#ifdef CYCLAGE_HOR
          if (pbc_mode && (i==CLEO-NB_MVT_EO-1)) ii2 -= CLEO-2*NB_MVT_EO;
#endif
          int32_t mvt_sol_EH = mvt_solid(ii2) ? (MVT_EH >> 6) : 0;  //CelMvt[ii-CLEO+1] & MVT_SOLID;

          //OUEST-BAS
          ii2 = ii+CLEO-1;
#ifdef CYCLAGE_HOR
          if (pbc_mode && (i==NB_MVT_EO)) ii2 += CLEO-2*NB_MVT_EO;
#endif
          int32_t mvt_sol_OB = mvt_solid(ii2) ? (MVT_OB >> 6) : 0;  //CelMvt[ii+CLEO-1] & MVT_SOLID;

          //OUEST-HAUT
          ii2 = ii-CLEO-1;
#ifdef CYCLAGE_HOR
          if (pbc_mode && (i==NB_MVT_EO)) ii2 += CLEO-2*NB_MVT_EO;
#endif
          int32_t mvt_sol_OH = mvt_solid(ii2) ? (MVT_OH >> 6) : 0;  //CelMvt[ii-CLEO-1] & MVT_SOLID;
          //int32_t mvt_sol = (mvt_sol_E << 3) | (mvt_sol_O << 2) | (mvt_sol_B << 1) | mvt_sol_H;
          //int32_t mvt_sol = (mvt_sol_EB << 3) | (mvt_sol_OB << 2) | (mvt_sol_EH << 1) | mvt_sol_OH;
          int32_t mvt_sol = mvt_sol_EB | mvt_sol_EH | mvt_sol_OB | mvt_sol_OH;
          UnsetMask(CelMvt2[ii], MVT_FAST);
          SetMask(CelMvt2[ii], Collisions_SolFast[mvt_sol][CelMvt[ii] & MVT_FAST]);
          //mvt_sol_fast = Collisions_SolFast[mvt_sol][CelMvt[ii] & MVT_FAST];  //debug
        }
#endif //SOLFAST

#ifdef SF_REAL  //Collisions realistes (solide-fluide)
        const float slow_val = SP8; //0.01; //SP8; //~0.38
        const float fast_val = CP8-SP8; //0.01; //CP8-SP8; //~0.54
        // normale precalculee
        int32_t ni = (i/NB_MVT_EO) - 1 + (k*DIST_MVT_NS)*(L-2);
        //LogPrintf("ni = %d \n", ni);
        float nx = norm2d[ni].x;
        float ny = norm2d[ni].y;
        //cas d'une particule lente
        //if (0){
        if (mvt_solid(ii) && (NbMvtCell[CelMvt[ii] & MVT_SLOW]==1)){
          MvtField mvt2 = 0;
          sfmod = 0;
          // calcul de la direction apres collision
          switch (CelMvt[ii] & MVT_SLOW)
          {
            case MVT_E:
              //LogPrintf("i=%d  k=%d   nx = %f   ny = %f\n", i, k, nx, ny);
              //if (fabs(ny) > slow_coef*fabs(nx)){
              if (fabs(ny) > slow_val){
              //if (ny){
                mvt2 = (ny > 0) ? MVT_H : MVT_B;
                sfmod = 1;
              }
              break;

            case MVT_O:
              if (fabs(ny) > slow_val){
              //if (ny){
                mvt2 = (ny > 0) ? MVT_H : MVT_B;
                sfmod = 1;
              }
              break;

            case MVT_B:
              /*if (fabs(nx) > slow_val){
              //if (nx){
                mvt2 = (nx > 0) ? MVT_E : MVT_O;
                sfmod = 1;
          }*/
              break;

            case MVT_H:
              //on conserve le rebond classique (MVT_B), on suppose qu'il s'agit du plafond ...
              break;
          }
          if (sfmod){
            // verification de l'absence de particule solide dans la direction calculee
            switch (mvt2)
            {
              case MVT_E:
                ii2 = ii+1;
#ifdef CYCLAGE_HOR
                if (pbc_mode && (i==CLEO-NB_MVT_EO-1)) ii2 -= CLEO-2*NB_MVT_EO;
#endif
                break;
              case MVT_O:
                ii2 = ii-1;
#ifdef CYCLAGE_HOR
                if (pbc_mode && (i==NB_MVT_EO)) ii2 += CLEO-2*NB_MVT_EO;
#endif
                break;
              case MVT_B:
                ii2 = ii+CLEO;
                break;
              case MVT_H:
                ii2 = ii-CLEO;
                break;
            }
            if (mvt_solid(ii2)) sfmod = 0;
          }
          if (sfmod){
            // modification des collisions par defaut
            //LogPrintf("collision modifiee : %d\n", mvt2);
            UnsetMask(CelMvt2[ii], MVT_SLOW);
            SetMask(CelMvt2[ii], mvt2);
            //mvt_sfr_slow = mvt2;  //debug
          }
        }

        //cas d'une particule rapide
        if (mvt_solid(ii) && (NbMvtCell[CelMvt[ii] & MVT_FAST]==1)){
          MvtField mvt2 = 0;
          sfmod = 0;
          // calcul de la direction apres collision
          switch (CelMvt[ii] & MVT_FAST)
          {
            case MVT_EB:
              if (fabs(nx+ny) > fast_val){
              //if (nx+ny){
              //if (1){
                mvt2 = (nx+ny > 0) ? MVT_EH : MVT_OB;
                sfmod = 1;
              }
              break;

            case MVT_EH:
              //if (!(ny-nx))
              //  continue; //on conserve le rebond classique (MVT_OB)
              //else
                //mvt2 = (ny-nx > 0) ? MVT_EH : MVT_OB;
              mvt2 = (j <= h_ceil) ? MVT_EB : MVT_OH; //selon l'altitude on considere que l'on rebondit sur le plafond ou sur une surface quasi-verticale
              sfmod = 1;
              break;

            case MVT_OB:
              if (fabs(ny-nx) > fast_val){
              //if (ny-nx){
                mvt2 = (ny-nx > 0) ? MVT_OH : MVT_EB;
                sfmod = 1;
              }
              break;

            case MVT_OH:
              mvt2 = (j <= h_ceil) ? MVT_OB : MVT_EH; //selon l'altitude on considere que l'on rebondit sur le plafond ou sur une surface quasi-verticale
              sfmod = 1;
          }
          if (sfmod){
            // verification de l'absence de particule solide dans la direction calculee
            switch (mvt2)
            {
              case MVT_EB:
                ii2=ii+1+CLEO;
#ifdef CYCLAGE_HOR
                if (pbc_mode && (i==CLEO-NB_MVT_EO-1)) ii2 -= CLEO-2*NB_MVT_EO;
#endif
                break;
              case MVT_EH:
                ii2=ii+1-CLEO;
#ifdef CYCLAGE_HOR
                if (pbc_mode && (i==CLEO-NB_MVT_EO-1)) ii2 -= CLEO-2*NB_MVT_EO;
#endif
                break;
              case MVT_OB:
                ii2=ii-1+CLEO;
#ifdef CYCLAGE_HOR
                if (pbc_mode && (i==NB_MVT_EO)) ii2 += CLEO-2*NB_MVT_EO;
#endif
                break;
              case MVT_OH:
                ii2=ii-1-CLEO;
#ifdef CYCLAGE_HOR
                if (pbc_mode && (i==NB_MVT_EO)) ii2 += CLEO-2*NB_MVT_EO;
#endif
                break;
            }
            if (mvt_solid(ii2)) sfmod = 0;
          }
          if (sfmod){
            // modification des collisions par defaut
            //LogPrintf("collision modifiee : %d (k=%d   j=%d   i=%d)\n", mvt2, k, j, i);
            UnsetMask(CelMvt2[ii], MVT_FAST);
            SetMask(CelMvt2[ii], mvt2);
            //mvt_sfr_fast = mvt2;  //debug
          }
        }
#endif //SF_REAL

        //comparaison SOL vs SF_REAL
        //if (mvt_sol_slow != mvt_sfr_slow) {LogPrintf("col_iter = %d   ii = %d   CelMvt[ii] = %d   CelMvt2[ii] = %d   mvt_sol_slow = %d   mvt_sfr_slow = %d   nx = %f   ny = %f   ni = %d\n", col_iter, ii, CelMvt[ii], CelMvt2[ii], mvt_sol_slow, mvt_sfr_slow, nx, ny, ni);}
        //if (mvt_sol_fast != mvt_sfr_fast) {LogPrintf("col_iter = %d   ii = %d   CelMvt[ii] = %d   CelMvt2[ii] = %d   mvt_sol_fast = %d   mvt_sfr_fast = %d   nx = %f   ny = %f   ni = %d\n", col_iter, ii, CelMvt[ii], CelMvt2[ii], mvt_sol_fast, mvt_sfr_fast, nx, ny, ni);}

#if !defined(CYCLAGE_MVT)
        //if ((CelMvt[ix] & MVT_OUT) && mvt(ix)) nb_mvt--;
        if ((CelMvt[ii] & MVT_OUT) && mvt(ii)){
          NbMvt[k] -= NbMvtCell[CelMvt[ii]];
          nb_mvt_out += NbMvtCell[CelMvt[ii]];
        }
    //if ((CelMvt[ix] & MVT_OUT) && (CelMvt[ix] & (MVT_E | MVT_O | MVT_B | MVT_H))) nb_mvt--;
#endif
      }
      ii += NB_MVT_EO;
    }
    ii+=CLEO;
  }
  col_iter++;
}

void propagate_mvtcel(int32_t i, int32_t ii)
{
  int32_t ii0;

  CelMvt[ii]=0;

  //SOLID+OUT
  CelMvt[ii] |= (CelMvt2[ii] & (MVT_SOLID | MVT_OUT));

  //EST
  ii0 = ii-1;
#ifdef CYCLAGE_MVT
  if (i==NB_MVT_EO) ii0 += CLEO-2*NB_MVT_EO;
#endif
  CelMvt[ii] |= (CelMvt2[ii0] & MVT_E);

  //OUEST
  ii0 = ii+1;
#ifdef CYCLAGE_MVT
  if (i==CLEO-NB_MVT_EO-1) ii0 -= CLEO-2*NB_MVT_EO;
#endif
  CelMvt[ii] |= (CelMvt2[ii0] & MVT_O);

  //BAS
  ii0 = ii-CLEO;
  CelMvt[ii] |= (CelMvt2[ii0] & MVT_B);

  //HAUT
  ii0 = ii+CLEO;
  CelMvt[ii] |= (CelMvt2[ii0] & MVT_H);

  //EST+BAS
  ii0 = ii-CLEO-1;
#ifdef CYCLAGE_MVT
  if (i==NB_MVT_EO) ii0 += CLEO-2*NB_MVT_EO;
#endif
  CelMvt[ii] |= (CelMvt2[ii0] & MVT_EB);

  //EST+HAUT
  ii0 = ii+CLEO-1;
#ifdef CYCLAGE_MVT
  if (i==NB_MVT_EO) ii0 += CLEO-2*NB_MVT_EO;
#endif
  CelMvt[ii] |= (CelMvt2[ii0] & MVT_EH);

  //OUEST+BAS
  ii0 = ii-CLEO+1;
#ifdef CYCLAGE_MVT
  if (i==CLEO-NB_MVT_EO-1) ii0 -= CLEO-2*NB_MVT_EO;
#endif
  CelMvt[ii] |= (CelMvt2[ii0] & MVT_OB);

  //OUEST+HAUT
  ii0 = ii+CLEO+1;
#ifdef CYCLAGE_MVT
  if (i==CLEO-NB_MVT_EO-1) ii0 -= CLEO-2*NB_MVT_EO;
#endif
  CelMvt[ii] |= (CelMvt2[ii0] & MVT_OH);

}

#ifdef USE_AVX
void propagate_mvtcel_vec(int32_t i, int32_t ii)
{
  static int32_t step = sizeof(__m256d)/sizeof(MvtField);
  static int32_t start = 1;
  static __m256i mm_maskv_so, mm_maskv_e, mm_maskv_w, mm_maskv_u, mm_maskv_d, mm_maskv_ue, mm_maskv_uw, mm_maskv_de, mm_maskv_dw;
  __m256i *mm_valp, mm_a, mm_b, mm_res, *mm_curp;

  if (start){
#ifdef OPENMP
    int32_t tid=omp_get_thread_num();
    //LogPrintf("thread_num=%d\n", tid);
    if (tid==0)
#endif
    {
      start = 0;
      LogPrintf("AVX optimization: %d MVT cells per vector\n", step);
    }
    mm_maskv_so = _mm256_set1_epi16(MVT_SOLID | MVT_OUT);
    mm_maskv_e = _mm256_set1_epi16(MVT_E);
    mm_maskv_w = _mm256_set1_epi16(MVT_O);
    mm_maskv_u = _mm256_set1_epi16(MVT_H);
    mm_maskv_d = _mm256_set1_epi16(MVT_B);
    mm_maskv_ue = _mm256_set1_epi16(MVT_EH);
    mm_maskv_uw = _mm256_set1_epi16(MVT_OH);
    mm_maskv_de = _mm256_set1_epi16(MVT_EB);
    mm_maskv_dw = _mm256_set1_epi16(MVT_OB);
  }

  if (!CheckAlign32(CelMvt2+ii)){
    ErrPrintf("ERROR: address of (CelMvt2+ii) is misaligned for AVX (256b)\n");
    exit(-1);
  };


  /// SOLID+OUT
  mm_valp = (__m256i*) (CelMvt2+ii);
  mm_a = _mm256_load_si256(mm_valp);
  mm_res = (__m256i) _mm256_and_pd((__m256d)mm_a, (__m256d)mm_maskv_so);

  /// EAST
  mm_valp = (__m256i*) (CelMvt2+ii-1);
  mm_a = _mm256_loadu_si256(mm_valp);
  mm_b = (__m256i) _mm256_and_pd((__m256d)mm_a, (__m256d)mm_maskv_e);
  mm_res = (__m256i) _mm256_or_pd((__m256d)mm_b, (__m256d)mm_res);

  /// WEST
  mm_valp = (__m256i*) (CelMvt2+ii+1);
  mm_a = _mm256_loadu_si256(mm_valp);
  mm_b = (__m256i) _mm256_and_pd((__m256d)mm_a, (__m256d)mm_maskv_w);
  mm_res = (__m256i) _mm256_or_pd((__m256d)mm_b, (__m256d)mm_res);

  /// UP
  mm_valp = (__m256i*) (CelMvt2+ii+CLEO);
  mm_a = _mm256_loadu_si256(mm_valp);
  mm_b = (__m256i) _mm256_and_pd((__m256d)mm_a, (__m256d)mm_maskv_u);
  mm_res = (__m256i) _mm256_or_pd((__m256d)mm_b, (__m256d)mm_res);

  /// DOWN
  mm_valp = (__m256i*) (CelMvt2+ii-CLEO);
  mm_a = _mm256_loadu_si256(mm_valp);
  mm_b = (__m256i) _mm256_and_pd((__m256d)mm_a, (__m256d)mm_maskv_d);
  mm_res = (__m256i) _mm256_or_pd((__m256d)mm_b, (__m256d)mm_res);

  /// UP-EAST
  mm_valp = (__m256i*) (CelMvt2+ii+CLEO-1);
  mm_a = _mm256_loadu_si256(mm_valp);
  mm_b = (__m256i) _mm256_and_pd((__m256d)mm_a, (__m256d)mm_maskv_ue);
  mm_res = (__m256i) _mm256_or_pd((__m256d)mm_b, (__m256d)mm_res);

  /// UP-WEST
  mm_valp = (__m256i*) (CelMvt2+ii+CLEO+1);
  mm_a = _mm256_loadu_si256(mm_valp);
  mm_b = (__m256i) _mm256_and_pd((__m256d)mm_a, (__m256d)mm_maskv_uw);
  mm_res = (__m256i) _mm256_or_pd((__m256d)mm_b, (__m256d)mm_res);

  /// DOWN-EAST
  mm_valp = (__m256i*) (CelMvt2+ii-CLEO-1);
  mm_a = _mm256_loadu_si256(mm_valp);
  mm_b = (__m256i) _mm256_and_pd((__m256d)mm_a, (__m256d)mm_maskv_de);
  mm_res = (__m256i) _mm256_or_pd((__m256d)mm_b, (__m256d)mm_res);

  /// DOWN-WEST
  mm_valp = (__m256i*) (CelMvt2+ii-CLEO+1);
  mm_a = _mm256_loadu_si256(mm_valp);
  mm_b = (__m256i) _mm256_and_pd((__m256d)mm_a, (__m256d)mm_maskv_dw);
  mm_res = (__m256i) _mm256_or_pd((__m256d)mm_b, (__m256d)mm_res);

  /// UPDATE CURRENT VECTOR
  mm_curp = (__m256i*) (CelMvt+ii);
  _mm256_store_si256(mm_curp, mm_res);
}
#endif // USE_AVX

void do_propagations()
{
  int32_t i, j, k, ii;
  int32_t n=0;
#ifdef USE_AVX
  static int32_t step = sizeof(__m256d)/sizeof(MvtField);
#endif

#ifdef OPENMP
#pragma omp parallel for private(i,j,k,ii)
#endif
  for(k=CLN; k<=CLS; k++){ // profondeur
    for(j=1; j<CH-1; j++){ // hauteur
      i = NB_MVT_EO;
      ii = CLEO*CH*k + CLEO*j + i;
      while(i<CLEO-NB_MVT_EO){ //largeur
#ifdef USE_AVX
        if (CheckAlign32(CelMvt+ii) && (i>1) && (i+step < CLEO-NB_MVT_EO)){
          propagate_mvtcel_vec(i ,ii);
          i += step;
          ii += step;
        }
        else{
          propagate_mvtcel(i++, ii++);
        }
#else
        propagate_mvtcel(i++, ii++);
#endif // USE_AVX
      }
    }
  }

  //Retroaction transport/matiere
  /*
  ix = 0;
  for(j=0; j<H; j++){ //hauteur
    for(i=0; i<L; i++, ix++){ //largeur
      if (CelMvt[ix] & ((-1) ^ MVT_SOLID)) { if (TE[ix].celltype == MOINS) {TE[ix].celltype = ZERO; n++;} }
      else if (TE[ix].celltype == ZERO) {TE[ix].celltype = MOINS; n++;}
    }
  }
  */

  //LogPrintf("propagation de %d cellules\n", n);

#ifdef MVT_REGUL
  do_regul_mvt();
#endif
#if defined(MVT_HFOR) || defined(MVT_VFOR) || defined(MVT_GFOR) //forcage flux
  do_force_mvt();
#endif

  if (lgca_mixing) do_mixing();

#ifdef DUMP_SIGNATURE
  //empreinte gazeuse
  if (opt_info && !(col_iter & 0xf))dump_signature_mvt(); //1 cycle sur 16
#endif

  if (opt_info) dump_mvt_in_out();
}

void do_mvt_in(int32_t ix)
{
  float aleat = drand48();
  int32_t k = (int)(ix/CHL);
  //MvtField mvt0 = CelMvt[ix];
  //injection d'une particule lente
  if (aleat < 0.5*DENSITE)
  if (!(CelMvt[ix] & (MVT_SOLID | MVT_OUT | MVT_E | MVT_O))){
    CelMvt[ix] |= MVT_E;
    NbMvt[k]++;
    nb_mvt_in++;
  }

  //injection d'une particule rapide
  if (aleat < 0.25*DENSITE){
    if (!(CelMvt[ix] & (MVT_SOLID | MVT_OUT | MVT_EB | MVT_OH))){
      CelMvt[ix] |= MVT_EB;
      NbMvt[k]++;
      nb_mvt_in++;
    }
  }
  else if (aleat < 0.5*DENSITE){
    if (!(CelMvt[ix] & (MVT_SOLID | MVT_OUT | MVT_EH | MVT_OB))){
      CelMvt[ix] |= MVT_EH;
      NbMvt[k]++;
      nb_mvt_in++;
    }
  }
}

#ifdef MVT_REGUL
void do_regul_mvt()
{
  int32_t k, ix, y, cpt;

  //LogPrintf("do_regul_mvt : %d -> %d\n", NbMvt, NbMvt0);

  for(k=CLN; k<=CLS; k++){
    //on injecte des particules sur le bord gauche
    cpt = 0;
    while ((NbMvt[k] < NbMvt0[k]) && (cpt<CH*2/*CH/2*/)){
      y = CH*drand48();
      ix = NB_MVT_EO+y*CLEO+k*CHL;
      do_mvt_in(ix);
      cpt++;
    }

    if (NbMvt[k] < NbMvt0[k]*0.9) LogPrintf("NbMvt[%d] = %ld\n", k, NbMvt[k])
  }

  //LogPrintf("fin regul\n");
}
#endif

void do_force_mvt()
{
  int32_t x, y, z;
  int32_t ix;
  // injection flux maximum sur le bord gauche
  /*for (y=0; y<CH; y++){
    ix = 1+y*CLEO;
    do_mvt_in(ix);
  }*/
#ifdef MVT_HFOR
  // forcage flux ouest-est sur la couche superieure
  //for(y=1; y<CH/2; y++)

  for(z=CLN; z<=CLS; z++){
    for(x=0; x<CLEO; x++){
      ix = x+CLEO*NB_MVT_VER*(h_ceil+1)+z*CHL;
      if ((CelMvt[ix] & (MVT_E | MVT_O)) == MVT_O){
        UnsetMask(CelMvt[ix], MVT_O);
        SetMask(CelMvt[ix], MVT_E);
      }
      if ((CelMvt[ix] & (MVT_EH | MVT_OH)) == MVT_OH){
        UnsetMask(CelMvt[ix], MVT_OH);
        SetMask(CelMvt[ix], MVT_EH);
      }
      if ((CelMvt[ix] & (MVT_EB | MVT_OB)) == MVT_OB){
        UnsetMask(CelMvt[ix], MVT_OB);
        SetMask(CelMvt[ix], MVT_EB);
      }
    }
  }


  // forcage flux est-ouest sur la couche superieure
  /*  for(x=0; x<L; x++){
      ix = x+L;
      if ((CelMvt[ix] & (MVT_E | MVT_O)) == MVT_E){
        CelMvt[ix] &= (-1) ^ MVT_E;
        CelMvt[ix] |= MVT_O;
      }
      if ((CelMvt[ix] & (MVT_EH | MVT_OH)) == MVT_EH){
        CelMvt[ix] &= (-1) ^ MVT_EH;
        CelMvt[ix] |= MVT_OH;
      }
      if ((CelMvt[ix] & (MVT_EB | MVT_OB)) == MVT_EB){
        CelMvt[ix] &= (-1) ^ MVT_EB;
        CelMvt[ix] |= MVT_OB;
      }
  }*/

  // forcage flux ouest-est sur la couche inferieure
  /*for(x=0; x<L; x++){
    ix = x+HL-2*L;
    if ((CelMvt[ix] & (MVT_E | MVT_O)) == MVT_O){
      CelMvt[ix] &= (-1) ^ MVT_O;
      CelMvt[ix] |= MVT_E;
    }
    if ((CelMvt[ix] & (MVT_EH | MVT_OH)) == MVT_OH){
      CelMvt[ix] &= (-1) ^ MVT_OH;
      CelMvt[ix] |= MVT_EH;
    }
    if ((CelMvt[ix] & (MVT_EB | MVT_OB)) == MVT_OB){
      CelMvt[ix] &= (-1) ^ MVT_OB;
      CelMvt[ix] |= MVT_EB;
    }
  }*/
#endif

#ifdef MVT_VFOR
  // forcage partiel du flux sur le bord gauche
  int32_t cpt_vfor=0;

  for(z=CLN; z<=CLS; z++){
    for(x=0; x<1; x++){
      ix = z*CHL+5*CLEO+NB_MVT_EO+x;
      //for (y=5; y<CH-5; y++, ix+=CLEO){
      for (y=0; y<CH; y++, ix+=CLEO){
        //if (drand48()<1.0){
        //forcage
        if (drand48() < intensite_vfor){
        //if (1){
          if ((CelMvt[ix] & (MVT_E | MVT_O)) == MVT_O){
            UnsetMask(CelMvt[ix], MVT_O);
            SetMask(CelMvt[ix], MVT_E);
            cpt_vfor++;
          }

          if ((CelMvt[ix] & (MVT_EH | MVT_OH)) == MVT_OH){
            UnsetMask(CelMvt[ix], MVT_OH);
            SetMask(CelMvt[ix], MVT_EH);
            cpt_vfor++;
          }

          if ((CelMvt[ix] & (MVT_EB | MVT_OB)) == MVT_OB){
            UnsetMask(CelMvt[ix], MVT_OB);
            SetMask(CelMvt[ix], MVT_EB);
            cpt_vfor++;
          }
        }
        //if (!mvt_solid(ix)) CelMvt[ix] = MVT_E;
        //if (!mvt_solid(ix+1)) CelMvt[ix+1] = MVT_E;
        //if (!mvt_solid(ix+2)) CelMvt[ix+2] = MVT_E;
      }
    }
  }
#endif //MVT_VFOR

#ifdef MVT_GFOR
  // forcage global du flux
  int32_t cpt_for = 0;
  int32_t nb_mvt_for = nb_cel_fluide * intensite_gfor;

  //LogPrintf ("forcage sur %d cellules MVT\n", nb_mvt_for);

  while (nb_mvt_for > 0){
    //ix = drand48()*CHLD;
    ix = CLN*CHL + drand48()*((CLS-CLN+1)*CHL);
    if (CelMvt[ix] & MVT_SOLID) continue;
    if ((CelMvt[ix] & (MVT_E | MVT_O)) == MVT_O){
      UnsetMask(CelMvt[ix], MVT_O);
      SetMask(CelMvt[ix], MVT_E);
      cpt_for++;
    }

    if ((CelMvt[ix] & (MVT_EH | MVT_OH)) == MVT_OH){
      UnsetMask(CelMvt[ix], MVT_OH);
      SetMask(CelMvt[ix], MVT_EH);
      cpt_for++;
    }

    if ((CelMvt[ix] & (MVT_EB | MVT_OB)) == MVT_OB){
      UnsetMask(CelMvt[ix], MVT_OB);
      SetMask(CelMvt[ix], MVT_EB);
      cpt_for++;
    }
    nb_mvt_for--;
  }
  //LogPrintf ("forcage de %d particules MVT\n", cpt_for);
#endif
}

//#ifdef BROUILLAGE_MVT
void do_mixing()
{
  int32_t i, y, z, y1, y2;
  int32_t cix, cix1, cix2, mvt;

  for(z=CLN; z<=CLS; z++){
    for (i=0; i<(CH/2); i++){

      //brouillage directionnel (bord gauche)
      y = 5+(CH-10)*drand48();
      cix = z*CHL+y*CLEO+NB_MVT_EO;
      mvt = CelMvt[cix] & (MVT_EB | MVT_EH);
      if (mvt == MVT_EB){
        UnsetMask(CelMvt[cix], MVT_EB);
        SetMask(CelMvt[cix], MVT_EH);
      }
      else if (mvt == MVT_EH){
        UnsetMask(CelMvt[cix], MVT_EH);
        SetMask(CelMvt[cix], MVT_EB);
      }

      //brouillage positionnel (bord gauche)
      y1 = 5+(CH-10)*drand48();
      y2 = 5+(CH-10)*drand48();
      cix1 = z*CHL+y1*CLEO+NB_MVT_EO;
      cix2 = z*CHL+y2*CLEO+NB_MVT_EO;
      if (!mvt_solid(cix1) && !mvt_solid(cix2)){
        //permutation
        mvt = CelMvt[cix1];
        CelMvt[cix1] = CelMvt[cix2];
        CelMvt[cix2] = mvt;
      }
    }
  }
}
//#endif


void compute_vel(int8_t flag_interp)
{
  static int32_t cpt=0;
  //static int32_t cpt_dump=0;
  int32_t i, j, k, ix, cpt1;
  int32_t *pvx, *pvy;
  int32_t vel_size;
  float fluid_ratio;

#ifndef VEL_SLIDE
  vel_size = VHLD;
#else
  vel_size = CHLD;
#endif
  //LogPrintf("vel_size=%d\n", vel_size);

  if (!Velx){
    AllocMemoryPrint("Velx", Velx, int, vel_size);
    ResetMemory(Velx, int, vel_size);
  }
  if (!Vely){
    AllocMemoryPrint("Vely", Vely, int, vel_size);
    ResetMemory(Vely, int, vel_size);
  }
  if (!Velx_time){
    if (VSTEP_TIME > 1){
      AllocMemoryPrint("Velx_time", Velx_time, int, vel_size*VSTEP_TIME);
      ResetMemory(Velx_time, int, vel_size*VSTEP_TIME);
    }
    else{
      Velx_time = Velx;
    }
  }
  if (!Vely_time){
    if (VSTEP_TIME > 1){
      AllocMemoryPrint("Vely_time", Vely_time, int, vel_size*VSTEP_TIME);
      ResetMemory(Vely_time, int, vel_size*VSTEP_TIME);
      PrintTotalMemory();
    }
    else{
      Vely_time = Vely;
    }
  }

#ifdef VEL_SLIDE
  /// sliding window
  if (!Veln){
    AllocMemoryPrint("Veln", Veln, int, vel_size);
    ResetMemory(Veln, int, vel_size);
  }
  /// sliding vectors
  if (!Velx_vector){
    AllocMemoryPrint("Velx_vector", Velx_vector, int, CHLD);
    ResetMemory(Velx_vector, int, CHLD);
  }
  if (!Vely_vector){
    AllocMemoryPrint("Vely_vector", Vely_vector, int, CHLD);
    ResetMemory(Vely_vector, int, CHLD);
  }
  if (!Veln_vector){
    AllocMemoryPrint("Veln_vector", Veln_vector, int, CHLD);
    ResetMemory(Veln_vector, int, CHLD);
  }
#endif

#ifndef VEL_SLIDE
  /// averaging over neighboring space by using fixed windows
  //LogPrintf("compute_vel: moyennage dans l'espace\n");
  int32_t iv, jv, imax, jmax;
  int32_t nb_mvt_fluide;

#ifdef OPENMP
#pragma omp parallel for private(i,j,k,iv,jv,ix,pvx,pvy,imax,jmax,nb_mvt_fluide,fluid_ratio)
#endif
  for(k=CLN; k<=CLS; k++){ //profondeur
    pvx = Velx_time + VHLD*cpt + VHL*k;
    pvy = Vely_time + VHLD*cpt + VHL*k;
    for(jv=0; jv<VH; jv++){ //hauteur
      for(iv=0; iv<VL; iv++, pvx++, pvy++){ //largeur
        // moyenne des vitesses dans un rectangle de taille VSTEP_H x VSTEP_L
        ix = k*CHL + jv*CLEO*VSTEP_H + iv*VSTEP_L;
        if (jv < VH-1) jmax = VSTEP_H; else jmax = CH-jv*VSTEP_H;
        if (iv < VL-1) imax = VSTEP_L; else imax = CLEO-iv*VSTEP_L;
        *pvx = *pvy = 0;
        nb_mvt_fluide = 0;
        for(j=0; j<jmax; j++, ix+=CLEO-imax){ //hauteur
          for(i=0; i<imax; i++, ix++){ //largeur
            *pvx += CelVelx[CelMvt[ix]];
            *pvy += CelVely[CelMvt[ix]];

            if (!(CelMvt[ix] & MVT_SOLID)){
              nb_mvt_fluide++;
            }
          }
        }
        //moyennage sur les cellules fluides uniquement
        fluid_ratio = (float)nb_mvt_fluide/(imax*jmax);
        if (fluid_ratio < 0.2) fluid_ratio = 0.2;
        *pvx = (*pvx)/fluid_ratio;
        *pvy = (*pvy)/fluid_ratio;
      }
    }
  }
#else
  /// averaging over neighboring space by using sliding windows

  int32_t *pwx, *pwy, *pwn, *pvn;
  int32_t *pwx0, *pwy0, *pwx1, *pwy1;
  int32_t iw, jw;
  int32_t i0 = VSTEP_L/2;
  int32_t j0 = VSTEP_H/2;
  int32_t i1 = VSTEP_L - i0;

  pvx = pvy = pvn = NULL;

#ifdef OPENMP
#pragma omp parallel for private(i,j,k,iw,jw,ix,pvx,pvy,pvn,pwx,pwy,pwn)
#endif
  for(k=CLN; k<=CLS; k++){ //depth
    pvx = Velx_vector + CHL*k;
    pvy = Vely_vector + CHL*k;
    pvn = Veln_vector + CHL*k;
    for(j=0; j<CH-VSTEP_H; j++){ //vertical slide
      //sliding window of size (VSTEP_H x VSTEP_L)
      pwx = Velx_time + CHLD*cpt + CHL*k + CLEO*(j+j0) + i0;
      pwy = Vely_time + CHLD*cpt + CHL*k + CLEO*(j+j0) + i0;
      pwn = Veln + CHL*k + CLEO*(j+j0) + i0;
      for(i=0; i<CLEO-VSTEP_L; i++, pvx++, pvy++, pvn++, pwx++, pwy++, pwn++){ //horizontal slide
        //averaging on the sliding window
        if (j==0){
          if (i==0){
            //upper left window of the plane
            ix = k*CHL;
            *pwx = *pwy = *pwn = 0;
            for(jw=0; jw<VSTEP_H; jw++, ix+=CLEO-VSTEP_L, pvx+=CLEO, pvy+=CLEO, pvn+=CLEO){ //height
              //sliding vectors of size VSTEP_L
              *pvx = *pvy = *pvn = 0;
              for(iw=0; iw<VSTEP_L; iw++, ix++){ //width
                if (CelMvt[ix] & MVT_SOLID) continue;
                *pvx += CelVelx[CelMvt[ix]];
                *pvy += CelVely[CelMvt[ix]];
                (*pvn)++;
              }
              *pwx += *pvx;
              *pwy += *pvy;
              *pwn += *pvn;
            }
          }
          else{
            //upper window
            ix = k*CHL + (VSTEP_L-1+i);
            *pwx = *pwy = *pwn = 0;
            for(jw=0; jw<VSTEP_H; jw++, ix+=CLEO){ //height
              pvx = Velx_vector + CHL*k + CLEO*jw + i;
              pvy = Vely_vector + CHL*k + CLEO*jw + i;
              pvn = Veln_vector + CHL*k + CLEO*jw + i;

              //slide vectors to the right
              *pvx = *(pvx-1) - CelVelx[CelMvt[ix-VSTEP_L]] + CelVelx[CelMvt[ix]];
              *pvy = *(pvy-1) - CelVely[CelMvt[ix-VSTEP_L]] + CelVely[CelMvt[ix]];
              *pvn = *(pvn-1) + (CelMvt[ix-VSTEP_L] & MVT_SOLID) - (CelMvt[ix] & MVT_SOLID);
              *pwx += *pvx;
              *pwy += *pvy;
              *pwn += *pvn;
            }
          }
        }
        else{
          if (i==0){
            //left window
            ix = k*CHL + (VSTEP_H-1+j)*CLEO;

            //new sliding vectors
            pvx = Velx_vector + CHL*k + CLEO*(VSTEP_H-1+j);
            pvy = Vely_vector + CHL*k + CLEO*(VSTEP_H-1+j);
            pvn = Veln_vector + CHL*k + CLEO*(VSTEP_H-1+j);
            *pvx = *pvy = *pvn = 0;

            for(iw=0; iw<VSTEP_L; iw++, ix++){ //width
              if (CelMvt[ix] & MVT_SOLID) continue;
              *pvx += CelVelx[CelMvt[ix]];
              *pvy += CelVely[CelMvt[ix]];
              (*pvn)++;
            }
          }
          else{
            //sliding window
            ix = k*CHL + (VSTEP_H-1+j)*CLEO + (VSTEP_L-1+i);

            //slide vectors to the right
            *pvx = *(pvx-1) - CelVelx[CelMvt[ix-VSTEP_L]] + CelVelx[CelMvt[ix]];
            *pvy = *(pvy-1) - CelVely[CelMvt[ix-VSTEP_L]] + CelVely[CelMvt[ix]];
            *pvn = *(pvn-1) + (CelMvt[ix-VSTEP_L] & MVT_SOLID) - (CelMvt[ix] & MVT_SOLID);
          }

          //slide down the window
          *pwx = *(pwx-CLEO) - *(pvx-CLEO*VSTEP_H) + *pvx;
          *pwy = *(pwy-CLEO) - *(pvy-CLEO*VSTEP_H) + *pvy;
          *pwn = *(pwn-CLEO) - *(pvn-CLEO*VSTEP_H) + *pvn;

          //full computation
          /*ix = k*CHL + j*L + i;
          *pwx = *pwy = *pwn = 0;
          for(jw=0; jw<VSTEP_H; jw++, ix+=CLEO-VSTEP_L){ //height
            //sliding vectors of size VSTEP_L
            pvx = Velx_vector + CHL*k + CLEO*(j+jw) + i;
            pvy = Vely_vector + CHL*k + CLEO*(j+jw) + i;
            pvn = Veln_vector + CHL*k + CLEO*(j+jw) + i;
            *pvx = *pvy = *pvn = 0;
            for(iw=0; iw<VSTEP_L; iw++, ix++){ //width
              if (CelMvt[ix] & MVT_SOLID) continue;
              *pvx += CelVelx[CelMvt[ix]];
              *pvy += CelVely[CelMvt[ix]];
              (*pvn)++;
            }
            *pwx += *pvx;
            *pwy += *pvy;
            *pwn += *pvn;
          }*/
        }
        /*if ((k==0) && (i<=2)){
          LogPrintf("i=%d  j=%d  pvn=%d  pvx=%d  pwn=%d  pwx=%d\n",i,j,*pvn,*pvx,*pwn,*pwx);
        }*/
        if (*pwn>VSTEP_L*VSTEP_H){
          ErrPrintf("ERROR: compute_vel - i=%d  j=%d  k=%d  pwn=%d\n",i,j,k,*pwn);
          exit(-3);
        }
      }
    }
  }

  /// normalize with the number of fluid cells in sliding windows
#ifdef OPENMP
#pragma omp parallel for private(i,j,k,pwx,pwy,pwn,fluid_ratio)
#endif
  for(k=CLN; k<=CLS; k++){
    for(j=0; j<CH-VSTEP_H; j++){
      pwx = Velx_time + CHLD*cpt + CHL*k + CLEO*(j+j0) + i0;
      pwy = Vely_time + CHLD*cpt + CHL*k + CLEO*(j+j0) + i0;
      pwn = Veln + CHL*k + CLEO*(j+j0) + i0;
      for(i=0; i<CLEO-VSTEP_L; i++, pwx++, pwy++, pwn++){
        fluid_ratio = (float)(*pwn)/(VSTEP_L*VSTEP_H);
        //if ((k==0) && (i==0)) LogPrintf("fluid_ratio=%f\n",fluid_ratio);
        if (fluid_ratio < 0.2) fluid_ratio = 0.2;
        if (fluid_ratio < 1.0){
          *pwx = (*pwx)/fluid_ratio;
          *pwy = (*pwy)/fluid_ratio;
        }
      }
    }
  }

  /// extend the values up to the boundaries
  /// east-west
  for(k=CLN; k<=CLS; k++){
    for(j=0; j<CH-VSTEP_H; j++){
      //left boundary
      pwx = Velx_time + CHLD*cpt + CHL*k + CLEO*(j+j0);
      pwy = Vely_time + CHLD*cpt + CHL*k + CLEO*(j+j0);
      pwx0 = pwx + i0;
      pwy0 = pwy + i0;
      for (i=0; i<i0; i++){
        *(pwx+i) = *pwx0;
        *(pwy+i) = *pwy0;
      }
      //rigth boundary
      pwx += CLEO-1;
      pwy += CLEO-1;
      pwx1 = pwx - i1;
      pwy1 = pwy - i1;
      for (i=0; i<i1; i++){
        *(pwx-i) = *pwx1;
        *(pwy-i) = *pwy1;
      }
    }
  }
  /// top
  for(k=CLN; k<=CLS; k++){
    for(i=0; i<CLEO; i++){
      pwx = Velx_time + CHLD*cpt + CHL*k + i;
      pwy = Vely_time + CHLD*cpt + CHL*k + i;
      pwx0 = pwx + CLEO*j0;
      pwy0 = pwy + CLEO*j0;
      for (j=0; j<j0; j++){
        *(pwx+j*CLEO) = *pwx0;
        *(pwy+j*CLEO) = *pwy0;
      }
    }
  }
#endif


  /// averaging over time by using a circular buffer
  //LogPrintf("compute_vel: moyennage dans le temps\n");
  cpt1 = cpt+1;
  if (cpt1==VSTEP_TIME) cpt1 = 0;
  if (VSTEP_TIME > 1){
    for(i=0; i<vel_size; i++){
      Velx[i] += Velx_time[cpt*vel_size+i] - Velx_time[cpt1*vel_size+i];
      Vely[i] += Vely_time[cpt*vel_size+i] - Vely_time[cpt1*vel_size+i];
    }
    //if (Velx[0]) {LogPrintf("Velx[0]=%d\tcpt=%d\n", Velx[0], cpt); }
  }

  /// computing mean and max velocity values over all space
  /*if (lgca_ready)*/{
    //meanvel = 0;
    float vel_max = 0;
    int32_t vel_size_2d = vel_size/CLNS;
    int32_t vx, vy;
    float vel, vel_sum;
    vel_sum = 0;
    maxvel = 0;
#ifdef OPENMP
#pragma omp parallel for private(i,k,vx,vy,vel,vel_max) reduction(+:vel_sum)
#endif
    for (k=CLN; k<=CLS; k++){
      for(i=0; i<vel_size_2d; i++){
        vx = Velx[i+k*vel_size_2d];
        vy = Vely[i+k*vel_size_2d];
        vel = sqrt(vx*vx + vy*vy);
        vel_sum += vel;
        if (vel_max<vel) vel_max = vel;
      }
#pragma omp critical
      {
        if (maxvel<vel_max) maxvel = vel_max;
      }
    }
    meanvel = vel_sum/(vel_size_2d*(CLS-CLN+1)); //VHLD;
  }

  if (flag_interp) compute_vel_interp();

  cpt = cpt1;
}

/// interpolation of the velocity
void compute_vel_interp()
{
  int32_t i, j, k;
  int32_t ix, ix0, ix1;

  if (!Velx_interp){
    AllocMemoryPrint("Velx_interp", Velx_interp, float, HLD);
    ResetMemory(Velx_interp, float, HLD);
  }
  if (!Vely_interp){
    AllocMemoryPrint("Vely_interp", Vely_interp, float, HLD);
    ResetMemory(Vely_interp, float, HLD);
    PrintTotalMemory();
  }

#ifndef VEL_SLIDE
  int32_t iv, jv;
  int32_t *pvx, *pvy;
  int32_t Velx00, Velx01, Velx10, Velx11;
  int32_t Vely00, Vely01, Vely10, Vely11;
  int32_t a00, a01, a10, a11;
  int32_t imax = VSTEP/NB_MVT_EO;

  /// interpolation in each vertical plane
  //LogPrintf("interpolation dans chaque plan de collisions\n");
#ifdef OPENMP
#pragma omp parallel for private(i,j,k,iv,jv,ix,pvx,pvy,Velx00,Velx01,Velx10,Velx11,Vely00,Vely01,Vely10,Vely11,a00,a01,a10,a11)
#endif
  for(k=CLN; k<=CLS; k++){ //profondeur
    for(jv=0; jv<VH-1; jv++){ //hauteur
      pvx = Velx+k*VHL+jv*VL;
      pvy = Vely+k*VHL+jv*VL;
      Velx00 = Velx01 = Velx10 = Velx11 = 0;
      Vely00 = Vely01 = Vely10 = Vely11 = 0;
      for(iv=0; iv<VL-2; iv++, pvx++, pvy++){ //largeur
        // interpolation des vitesses dans un rectangle de taille (imax)xVSTEP
        ix = (LN+k*DIST_MVT_NS)*HL + jv*L*VSTEP + iv*imax;// + (VSTEP>>1)*L + (VSTEP>>1);
        Velx00 = *pvx;
        Velx10 = *(pvx+VL);
        Vely00 = *pvy;
        Vely10 = *(pvy+VL);
        if (imax > 1){
          Velx01 = *(pvx+1);
          Velx11 = *(pvx+VL+1);
          Vely01 = *(pvy+1);
          Vely11 = *(pvy+VL+1);
        }
        for(j=0; j<VSTEP; j++, ix+=L-imax){ //hauteur
          if (imax > 1){ //interpolation entre les 4 sommets du rectangle englobant
            for(i=0; i<imax; i++, ix++){ //largeur
              a00 = (VSTEP-j)*(imax-i);
              a01 = (VSTEP-j)*i;
              a10 = j*(imax-i);
              a11 = j*i;
              Velx_interp[ix] = (float)(Velx00*a00 + Velx01*a01 + Velx10*a10 + Velx11*a11)/(VSTEP*imax);
              Vely_interp[ix] = (float)(Vely00*a00 + Vely01*a01 + Vely10*a10 + Vely11*a11)/(VSTEP*imax);
            }
          }
          else{ //interpolation verticale entre 2 sommets
            Velx_interp[ix] = (float)(Velx00*(VSTEP-j) + Velx10*j)/VSTEP;
            Vely_interp[ix] = (float)(Vely00*(VSTEP-j) + Vely10*j)/VSTEP;
            ix+=imax;
          }
        }
      }

      // interpolation a cote du bord droit
      // on repete la derniere valeur interpolee
      // c'est plus facile et on reduit les effets de bord ...
      ix = (LN+k*DIST_MVT_NS)*HL + jv*L*VSTEP + iv*imax;// + (VSTEP>>1)*L + (VSTEP>>1);
      for(j=0; j<VSTEP; j++, ix+=L){ //hauteur
        for(i=0; i<L-iv*imax; i++){ //largeur
          Velx_interp[ix+i] = Velx_interp[ix-1];
          Vely_interp[ix+i] = Vely_interp[ix-1];
        }
      }
    }
  }

#else //VEL_SLIDE is defined

  /// no interpolation in the vertical planes : we simply copy Velx and Vely into Velx_interp and Vely_interp
#ifdef OPENMP
#pragma omp parallel for private(i,k,ix,ix0)
#endif
  for(k=CLN; k<=CLS; k++){ //profondeur
    ix = (LN+k*DIST_MVT_NS)*HL;
    ix0 = k*HL;
    for(i=0; i<HL; i++, ix++, ix0++){
      Velx_interp[ix] = Velx[ix0];
      Vely_interp[ix] = Vely[ix0];
    }
  }

#endif

#if (DIST_MVT_NS > 1)
  float fz0, fz1;
  int32_t kv, kvmax;

  /// north-south interpolation between all vertical planes
  //LogPrintf("interpolation nord-sud dans les plans intermediaires\n");
#ifdef OPENMP
#pragma omp parallel for private(i,j,k,kv,kvmax,ix,ix0,ix1,fz0,fz1)
#endif
  for(k=CLN; k<=CLS; k++){ /// vertical planes (lattice gas)
    kvmax = (k<CLS) ? DIST_MVT_NS : LNS-CLS*DIST_MVT_NS;
    for(kv=1; kv<kvmax; kv++){ /// intermediate planes (interpolation)
      ix = (LN+k*DIST_MVT_NS+kv)*HL;
      ix0 = (LN+k*DIST_MVT_NS)*HL;
      //ix1 = (LN+(k+1)*DIST_MVT_NS)*HL;
      ix1 = (k<CLS) ? (LN+(k+1)*DIST_MVT_NS)*HL : pbc_mode ? LN*HL : ix0;
      fz1 = (float)kv/kvmax;
      fz0 = 1.0-fz1;
      for(j=0; j<H; j++){
        for(i=0; i<L; i++, ix++, ix0++, ix1++){
          if (rot_map && OutOfSpace(i,LN+k*DIST_MVT_NS+kv)) continue;
          Velx_interp[ix] = Velx_interp[ix0]*fz0 + Velx_interp[ix1]*fz1;
          Vely_interp[ix] = Vely_interp[ix0]*fz0 + Vely_interp[ix1]*fz1;
        }
      }
    }
  }
#endif
}


void dump_mvt_in_out()
{
  static int32_t cpt = 0;
  static int32_t step = 0;
  FILE *fp;

  if (!cpt && !step){
    fp = fopen("MVT_IO.log","w");
    fprintf(fp,"      \t mvt_in    \t mvt_out   \t ratio\n");
    fclose(fp);
  }

  //fichier in-out
  if (step >= VSTEP_TIME){
    fp = fopen("MVT_IO.log","a");
    fprintf(fp,"%04d: \t%09ld \t%09ld \t%f\n", cpt, nb_mvt_in, nb_mvt_out, (nb_mvt_in) ? ((float)nb_mvt_out/nb_mvt_in) : 0);
    fclose(fp);
    step = nb_mvt_in = nb_mvt_out = 0;	//reset
    cpt++;
  }
  else
    step++;
}

void dump_densite()
{
  static int32_t cpt = 0;
  FILE *fp;
  int32_t ix, n;
  //int32_t nb_cel_fluide_dump = 0;

  //recalcul nb_mvt, nb_mvt_sol et nb_cel_fluide
  int32_t nb_mvt = 0; //nombre de particules mobiles
  int32_t nb_mvt_sol = 0; //nombre de particules mobiles localisees dans une cellule solide
  //for (ix=0; ix<CHLD; ix++){
  for (ix=CHL*CLN; ix<CHL*(CLS+1); ix++){
    n = NbMvtCell[CelMvt[ix]];
    nb_mvt += n;
    if ((CelMvt[ix] & MVT_SOLID)) nb_mvt_sol += n;
    //if (!(CelMvt[ix] & MVT_SOLID)) nb_cel_fluide_dump++;
  }
  //densite = (float)nb_mvt / nb_cel_fluide;
  densite = (float)(nb_mvt - nb_mvt_sol) / nb_cel_fluide;

  //fichier densite
  fp = fopen("DENSI.log","a");
  if (!cpt) fprintf(fp,"\tmvt cells \tfluid nodes \tdensity \ttrapped cells\n");
  fprintf(fp,"%04d: \t%09ld \t%09ld \t%f \t%09ld\n", cpt++, (long)nb_mvt, nb_cel_fluide, densite, (long)nb_mvt_sol);
  //fprintf(fp,"%03d : \t%09ld \t%09ld (%09ld) \t%f\n", cpt++, (long)nb_mvt, nb_cel_fluide, nb_cel_fluide_dump, densite);
  fclose(fp);
}

void dump_vel()
{
  static int32_t cpt = 0;
  static int8_t start = 1;
  FILE *fp;

  if ((!Velx_interp) || (!Vely_interp)) {cpt++; return;}

#ifdef STABILITY_ANALYSIS
  //recalcul maxvel et meanvel  (calcul plus precis mais plus long, car on elimine les cellules solides lors du moyennage)
  int32_t ix, ncel;
  float vx, vy, vv;
  maxvel = 0.0;
  meanvel = 0.0;
  //float meanvel2 = 0.0;
  ncel = 0;
  for (ix=0; ix<HLD; ix++){
    vx = Velx_interp[ix];
    vy = Vely_interp[ix];
    vv = vx*vx+vy*vy;
    if (vv>maxvel) maxvel=vv;
    if (Phase[TE[ix].celltype] == FLUID) {meanvel += vv; /*meanvel2 += sqrt(vv);*/ ncel++;}
  }
  maxvel=sqrt(maxvel);
  meanvel=sqrt(meanvel/ncel); //approximation (en toute rigueur, il faut calculer la moyenne des racines carrees => plus de temps de calcul)
  //meanvel2/=ncel; //valeur exacte
#endif

  //fichier VEL.log
  fp = fopen("VEL.log","a");
  if (start){
    start = 0;
    fprintf(fp,"      \t maxvel \t meanvel\n");
  }
  fprintf(fp,"%04d: \t  %04d  \t  %04d\n", cpt++, (int)maxvel, (int)meanvel);
  //fprintf(fp,"%03d : \t  %04d  \t  %04d  \t  %04d\n", cpt++, (int)maxvel, (int)meanvel, (int)meanvel2);
  //fprintf(fp,"nb_mvt = %lu\n", nb_mvt);
  fclose(fp);
}

void dump_mvt(int32_t cpt, int32_t unit)
{
  static int8_t filename[100];
  FILE *fp;

#if 1
  // binary  file
  //sprintf(filename,"HPP_%d_%d_%d.bin",CLNS,CH,CLEO);
  if (unit == UNIT_COMP)
    sprintf(filename,"HPP%04d.bin", cpt);
  else
    sprintf(filename,"HPP%05d_t0.bin", cpt);

  fp = fopen(filename,"w");
  if ( ! fp ){
    ErrPrintf("ERROR: cannot open binary file %s\n", filename);
    exit(-4);
  }


  //while (CelMvt_lock) {LogPrintf("lock3\n"); };

  //CelMvt_lock = 1;

  fwrite(CelMvt, sizeof(MvtField), CHLD, fp);

  //CelMvt_lock = 0;

  fclose(fp);
#endif

  //free(buf);
  /*if (inter){
    int8_t command[100];
    sprintf(command, "gzip %s", nom);
    system(command);
  }*/

  // fichier vitesses
//#if (VSTEP > 1)
#if 1
  if (Velx && Vely){
    if (unit == UNIT_COMP)
      sprintf(filename,"HPP%04d_%d_%d_%d.vel", cpt, CLNS, VH, VL);
    else
      sprintf(filename,"HPP%05d_%d_%d_%d.vel", cpt, CLNS, VH, VL);

    fp = fopen(filename,"w");
    if ( ! fp ){
      ErrPrintf("ERROR: cannot open VEL file %s\n", filename);
      exit(-4);
    }

    //compute_vel();

    fwrite(Velx, sizeof(int), VHLD, fp);
    fwrite(Vely, sizeof(int), VHLD, fp);

    fclose(fp);
  }
#endif
}

void dump_collisions()
{
  FILE *fp;
  int32_t mvt_in, mvt_out;

  fp = fopen("LGCA.log","w");

  //dump des regles de collision entre particules fluides
  fprintf(fp, "          IN            |          OUT\n");
  fprintf(fp, "E  W  B  T EB ET WB WT  |  E  W  B  T EB ET WB WT\n");
  fprintf(fp, "-------------------------------------------------\n");
  for(mvt_in=0; mvt_in<SIZE_MVT_FIELD; mvt_in+=4){
    mvt_out = Collisions[mvt_in];
    fprintf(fp, "%d  %d  %d  %d  %d  %d  %d  %d  |  %d  %d  %d  %d  %d  %d  %d  %d\n",
      (mvt_in & MVT_E)? 1 : 0, (mvt_in & MVT_O)? 1 : 0, (mvt_in & MVT_B)? 1 : 0, (mvt_in & MVT_H)? 1 : 0,
      (mvt_in & MVT_EB)? 1 : 0, (mvt_in & MVT_EH)? 1 : 0, (mvt_in & MVT_OB)? 1 : 0, (mvt_in & MVT_OH)? 1 : 0,
      (mvt_out & MVT_E)? 1 : 0, (mvt_out & MVT_O)? 1 : 0, (mvt_out & MVT_B)? 1 : 0, (mvt_out & MVT_H)? 1 : 0,
      (mvt_out & MVT_EB)? 1 : 0, (mvt_out & MVT_EH)? 1 : 0, (mvt_out & MVT_OB)? 1 : 0, (mvt_out & MVT_OH)? 1 : 0);
  }

#ifdef SOLSLOW
  //dump des regles de collision solide-fluide pour les particules lentes
  //selon la configuration du voisinage
  int32_t slow[4] = {MVT_E, MVT_O, MVT_B, MVT_H};
  int32_t mvt_sol_E, mvt_sol_O, mvt_sol_B, mvt_sol_H, mvt_sol, i;
  fprintf(fp, "\n   SOLID     |      IN      |      OUT\n");
  fprintf(fp, "E  O  B  H   | E  O  B  H   | E  O  B  H\n");
  for (mvt_sol=0; mvt_sol<16; mvt_sol++){
    mvt_sol_E = (mvt_sol & (MVT_E >> 2))? 1 : 0;
    mvt_sol_O = (mvt_sol & (MVT_O >> 2))? 1 : 0;
    mvt_sol_B = (mvt_sol & (MVT_B >> 2))? 1 : 0;
    mvt_sol_H = (mvt_sol & (MVT_H >> 2))? 1 : 0;
    for (i=0; i<4; i++){
      mvt_out = Collisions_SolSlow[mvt_sol][slow[i]];
      //mvt_in = MVT_SOLID | (1 << (2+i));
      //mvt_out = Collisions[mvt_in];
      fprintf(fp, "% d  %d  %d  %d  |  %d  %d  %d  %d  |  %d  %d  %d  %d\n",
              mvt_sol_E, mvt_sol_O, mvt_sol_B, mvt_sol_H,
              (i==0)? 1 : 0, (i==1)? 1 : 0, (i==2)? 1 : 0, (i==3)? 1 : 0,
              (mvt_out & MVT_E)? 1 : 0, (mvt_out & MVT_O)? 1 : 0, (mvt_out & MVT_B)? 1 : 0, (mvt_out & MVT_H)? 1 : 0);
    }
  }
#endif

#ifdef SOLFAST
  //dump des regles de collision solide-fluide pour les particules rapides (en diagonale)
  //selon la configuration du voisinage
  int32_t fast[4] = {MVT_EB, MVT_EH, MVT_OB, MVT_OH};
  int32_t mvt_sol_EB, mvt_sol_OB, mvt_sol_EH, mvt_sol_OH;//, mvt_sol, i;
  fprintf(fp, "\n   SOLID     |      IN      |      OUT\n");
  fprintf(fp, "EB EH OB OH  | EB EH OB OH  | EB EH OB OH\n");
  for (mvt_sol=0; mvt_sol<16; mvt_sol++){
    mvt_sol_EB = (mvt_sol & (MVT_EB >> 6))? 1 : 0;
    mvt_sol_EH = (mvt_sol & (MVT_EH >> 6))? 1 : 0;
    mvt_sol_OB = (mvt_sol & (MVT_OB >> 6))? 1 : 0;
    mvt_sol_OH = (mvt_sol & (MVT_OH >> 6))? 1 : 0;
    for (i=0; i<4; i++){
      mvt_out = Collisions_SolFast[mvt_sol][fast[i]];
      //mvt_in = MVT_SOLID | (1 << (6+i));
      //mvt_out = Collisions[mvt_in];
      fprintf(fp, "% d  %d  %d  %d  |  %d  %d  %d  %d  |  %d  %d  %d  %d\n",
              mvt_sol_EB, mvt_sol_EH, mvt_sol_OB, mvt_sol_OH,
              (i==0)? 1 : 0, (i==1)? 1 : 0, (i==2)? 1 : 0, (i==3)? 1 : 0,
              (mvt_out & MVT_EB)? 1 : 0, (mvt_out & MVT_EH)? 1 : 0, (mvt_out & MVT_OB)? 1 : 0, (mvt_out & MVT_OH)? 1 : 0);
    }
  }
#endif
  fclose(fp);
}



#ifdef DUMP_SIGNATURE
void dump_signature_mvt()
{
  FILE *fp;
  int32_t i;
  uint32_t sig, *aux;

  //calcul de la signature
  sig=0;
  aux = (unsigned int*)CelMvt;
  for(i=0; i<CHLD*sizeof(MvtField)/sizeof(unsigned int); i++, aux++) sig += (*aux); //sig = sig ^ (*aux);

  //dump
  fp = fopen("SIGN_HPP.log","a");
  if ( ! fp ){
    ErrPrintf("erreur ouverture fichier dump signature HPP\n");
    exit(-4);
  }
  fprintf(fp, "%08x : %08x\n", col_iter, sig);
  fclose(fp);

}
#endif

#endif //LGCA
