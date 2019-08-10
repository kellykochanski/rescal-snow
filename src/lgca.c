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


#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <assert.h>
#include <math.h>
//#include <stdint.h>

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

//#undef CYCLAGE_HOR  //pas de cyclage sur le flux (optionnel) -> sinon essoufflement du flux (frottements)

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
extern Cell  *TE;            // la 'terre'
extern int32_t Ncel[];                       // nombre de cellules par type
extern uint8_t opt_vel;       //affichage des vitesses moyennes
extern Vec2 *norm2d;              //normale a la surface
extern int32_t h_ceil;  //epaisseur du plafond
extern int8_t *psol;                // indicateur des plans solides verticaux est-ouest
extern const uint8_t Phase[MAX_CELL]; //phase (fluide ou solide) des types de cellules

MvtField *CelMvt; //cellules en mouvements
MvtField *CelMvt2;  //cellules en mouvements
MvtField Collisions[SIZE_MVT_FIELD];  //tableau des collisions fuilde-fluide
int32_t CelVelx[SIZE_MVT_FIELD];  //precomputed sums of horizontal velocities for a single node
int32_t CelVely[SIZE_MVT_FIELD];  //precomputed sums of vertical velocities for a single node

int32_t *Velx_time = NULL, *Vely_time = NULL; //vitesses locales instantanees
int32_t *Velx = NULL, *Vely = NULL; //vitesses locales moyennees (espace et temps)
int32_t *Veln = NULL; //sliding window for the number of fluid nodes
float *Velx_interp = NULL, *Vely_interp = NULL; //vitesses locales interpolees (espace) et moyennees (temps)
int32_t *Velx_vector = NULL, *Vely_vector = NULL; //sliding vectors of velocities
int32_t *Veln_vector = NULL; //sliding vector for the number of fluid nodes
int32_t CH, CLEO, CLNS, CHL, CHLD, CLN, CLS;     // dimensions du reseau de collisions (CLEO/CLNS = largeur du reseau de collisions en est-ouest/nord-sud)
int32_t VH, VL, VHL, VHLD;  // dimensions du tableau des vitesses
int32_t HLN;    //offset dans le reseau classique
int64_t nb_cel_fluide = 0; //nombre de cellules fluides
int64_t *NbMvt = NULL;  //nombre de particules en mouvement dans chaque couche du reseau
int64_t nb_mvt_in = 0;  //nombre de particules entrant dans le reseau
int64_t nb_mvt_out = 0; //nombre de particules sortant du reseau
#ifdef MVT_REGUL
int64_t *NbMvt0 = NULL; //nombre de particules en mouvement au depart dans chaque couche du reseau
#endif
uint8_t NbMvtCell[SIZE_MVT_FIELD];  //nombre de particules en mouvement dans chaque etat d'une cellule
float densite = 0; //densite moyenne de particules en mouvement
float maxvel = 0.0;  //vitesse interpolee maximale apres stabilisation
float meanvel = 0.0;  //vitesse interpolee moyenne apres stabilisation
int32_t col_iter = 0; //nombre de cycles de collisions
char *nom_fic_mvt = NULL; //nom du fichier VELOC (initial)
int32_t lgca_reset = 1;
int32_t lgca_speedup = 0; //initial speedup of the stabilization of lattice gas
int32_t lgca_mixing = 1; //mixing of the flux (on the left boundary)

void lecture_mvt();
void do_mvt_in(int32_t ix);
void do_regul_mvt();
void do_force_mvt();
void do_mixing();
void compute_vel_interp();
void dump_densite();
void dump_mvt_in_out();
void dump_collisions();


void params_collisions() {
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

void init_collisions() {
  int32_t i;

  /// size of the lattice
  CH = H * NB_MVT_VER;
  CLEO = L * NB_MVT_EO;
  CLNS = (LNS + DIST_MVT_NS - 1) / DIST_MVT_NS;
  CHL = CH * CLEO;

  CLN = 0;
  CLS = CLNS - 1;
  CHLD = CH * CLEO * CLNS;

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
  VH = (CH + VSTEP_H - 1) / VSTEP_H;
  VL = (CLEO + VSTEP_L - 1) / VSTEP_L;
  LogPrintf("VH = %d\n", VH);
  LogPrintf("VL = %d\n", VL);
  VHL = VH * VL;
  VHLD = VHL * CLNS;
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
  if (lgca_mixing) {
    LogPrintf("brouillage flux sur le bord gauche\n");
  }
#ifdef PLAF_REAL
  LogPrintf("rebonds realistes PLAF_REAL\n");
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
  AllocMemory(NbMvt, int64_t, CLNS);
  ResetMemory(NbMvt, int64_t, CLNS);

#ifdef MVT_REGUL
  AllocMemory(NbMvt0, int64_t, CLNS);
  ResetMemory(NbMvt0, int64_t, CLNS);
#endif
  PrintTotalMemory();

  if (nom_fic_mvt) {
    lecture_mvt();
  } else {
    init_mvt();
  }

  // initialisation tableau NbMvtCell[]
  for (i = 0; i < SIZE_MVT_FIELD; i++) {
    NbMvtCell[i] = 0;
    if (i & MVT_E) {
      NbMvtCell[i]++;
    }
    if (i & MVT_O) {
      NbMvtCell[i]++;
    }
    if (i & MVT_B) {
      NbMvtCell[i]++;
    }
    if (i & MVT_H) {
      NbMvtCell[i]++;
    }
    if (i & MVT_EB) {
      NbMvtCell[i]++;
    }
    if (i & MVT_EH) {
      NbMvtCell[i]++;
    }
    if (i & MVT_OB) {
      NbMvtCell[i]++;
    }
    if (i & MVT_OH) {
      NbMvtCell[i]++;
    }
  }

  /// number of gas particules
#ifdef MVT_REGUL
  for (k = CLN; k <= CLS; k++) {
    if (NbMvt[k]) { //TODO
      NbMvt0[k] = NbMvt[k];
    } else //approximation
#if defined(MODEL_DUN) || defined(MODEL_SNO)
      NbMvt0[k] = Ncel[EAUC] * DENSITE / CLNS;
#else
      NbMvt0[k] = Ncel[MOINS] * DENSITE / CLNS;
#endif
    LogPrintf("NbMvt0[%d] = %ld\n", k, NbMvt0[k]);
  }
#endif

  LogPrintf("nb_mvt_in = %" PRId64 "\n", nb_mvt_in);

  // initialisation of Collisions[]
  for (i = 0; i < SIZE_MVT_FIELD; i++) {
    if (i & MVT_OUT) {
      Collisions[i] = MVT_OUT;
    } else if (i & MVT_SOLID) {
      Collisions[i] = MVT_SOLID;
      if (i & MVT_E) {
        Collisions[i] |= MVT_O;
      }
      if (i & MVT_O) {
        Collisions[i] |= MVT_E;
      }
      if (i & MVT_B) {
        Collisions[i] |= MVT_H;
      }
      if (i & MVT_H) {
        Collisions[i] |= MVT_B;
      }

      if (i & MVT_EB) {
        Collisions[i] |= MVT_OH;  //MVT_OH;
      }
      if (i & MVT_EH) {
        Collisions[i] |= MVT_OB;  //MVT_OB;
      }
      if (i & MVT_OB) {
        Collisions[i] |= MVT_EH;  //MVT_EH;
      }
      if (i & MVT_OH) {
        Collisions[i] |= MVT_EB;  //MVT_EB;
      }
    } else {
      Collisions[i] = 0;

      if ((i & MVT_SLOW) == MVT_E) {
        if ((i & MVT_FAST) == MVT_OB) {
          Collisions[i] |= (MVT_O | MVT_EB);
        } else if ((i & MVT_FAST) == MVT_OH) {
          Collisions[i] |= (MVT_O | MVT_EH);
        } else {
          Collisions[i] |= MVT_E;
        }
      } else if ((i & MVT_SLOW) == MVT_O) {
        if ((i & MVT_FAST) == MVT_EB) {
          Collisions[i] |= (MVT_E | MVT_OB);
        } else if ((i & MVT_FAST) == MVT_EH) {
          Collisions[i] |= (MVT_E | MVT_OH);
        } else {
          Collisions[i] |= MVT_O;
        }
      } else if ((i & MVT_SLOW) == MVT_B) {
        if ((i & MVT_FAST) == MVT_EH) {
          Collisions[i] |= (MVT_H | MVT_EB);
        } else if ((i & MVT_FAST) == MVT_OH) {
          Collisions[i] |= (MVT_H | MVT_OB);
        } else {
          Collisions[i] |= MVT_B;
        }
      } else if ((i & MVT_SLOW) == MVT_H) {
        if ((i & MVT_FAST) == MVT_EB) {
          Collisions[i] |= (MVT_B | MVT_EH);
        } else if ((i & MVT_FAST) == MVT_OB) {
          Collisions[i] |= (MVT_B | MVT_OH);
        } else {
          Collisions[i] |= MVT_H;
        }
      } else if ((i & MVT_SLOW) == (MVT_E | MVT_O)) {
        Collisions[i] |= (MVT_H | MVT_B);
      } else if ((i & MVT_SLOW) == (MVT_H | MVT_B)) {
        Collisions[i] |= (MVT_E | MVT_O);
      } else {
        Collisions[i] |= (i & MVT_SLOW);
      }

      if ((i & MVT_FAST) == (MVT_EB | MVT_OH)) {
        Collisions[i] |= (MVT_EH | MVT_OB);
      } else if ((i & MVT_FAST) == (MVT_EH | MVT_OB)) {
        Collisions[i] |= (MVT_EB | MVT_OH);
      } else if (!(Collisions[i] & MVT_FAST)) {
        Collisions[i] |= (i & MVT_FAST);
      }
    }
  }

  if (opt_info) {
    dump_collisions();
  }
  

  // initialization of precomputed sums of velocities
  int32_t *pvx = CelVelx;
  int32_t *pvy = CelVely;
  for (i = 0; i < SIZE_MVT_FIELD; i++, pvx++, pvy++) {
#ifdef VEL_SLIDE
    if (i & MVT_SOLID) {
      continue;
    }
#endif
    if (i & MVT_SLOW) {
      if (i & MVT_E) {
        (*pvx)++;
      }
      if (i & MVT_O) {
        (*pvx)--;
      }
      if (i & MVT_B) {
        (*pvy)++;
      }
      if (i & MVT_H) {
        (*pvy)--;
      }
    }
    if (i & MVT_FAST) {
      if (i & MVT_EB) {
        (*pvx)++;
        (*pvy)++;
      }
      if (i & MVT_EH) {
        (*pvx)++;
        (*pvy)--;
      }
      if (i & MVT_OB) {
        (*pvx)--;
        (*pvy)++;
      }
      if (i & MVT_OH) {
        (*pvx)--;
        (*pvy)--;
      }
    }
  }

#ifdef DUMP_SIGNATURE
  //empreinte gazeuse
  if (opt_info) {
    dump_signature_mvt();
  }
#endif

}


// initialize one site of lattice gas
void init_mvt_cell(MvtField *mvt) {
  static MvtField mvt_slow[] = { MVT_E, MVT_O, MVT_B, MVT_H };
  static MvtField mvt_fast[] = { MVT_EB, MVT_EH, MVT_OB, MVT_OH };
  int32_t alea;
  alea = drand48() * 4;
  if (alea < 4) {
    *mvt = mvt_slow[alea];
    nb_mvt_in++;
  }
  alea = drand48() * 8;
  if (alea < 4) {
    *mvt |= mvt_fast[alea];
    nb_mvt_in++;
  }
}

// initialisation du tableau MvtCel des particules de gaz sur reseau
// initialization of the MvtCel array of gas particles on a network
void init_mvt() {
  int32_t i, j, k, ic;
  Cell *aux;
  MvtField *aux_mvt;

  LogPrintf("reset lattice gas\n");

  if (rot_map) {
    /// set the north-south boundaries (we have to keep an extra vertical plan on north and south sides for the interpolation)
    CLN = 0;
    for (k = 0; psol[LN + k * DIST_MVT_NS]; k++) {
      CLN = k;
    }
    for (; !psol[LN + k * DIST_MVT_NS]; k++);
    CLS = k;
    assert(CLS < CLNS);
  }

  // initialisation of CelMvt[]
  nb_cel_fluide = 0;
  for (k = CLN; k <= CLS; k++) { // depth
    aux = TE + HLN + k * HL * DIST_MVT_NS;
    aux_mvt = CelMvt + k * HL * NB_MVT_EO;
    for (j = 0; j < H; j++) { // height
      for (i = 0; i < L; i++, aux++) { // width
        for (ic = 0; ic < NB_MVT_EO; ic++, aux_mvt++) {
          if (Phase[aux->celltype] == FLUID) {
            *aux_mvt = 0;
            nb_cel_fluide++;
            init_mvt_cell(aux_mvt);
          } else if (Phase[aux->celltype] == SOLID) {
            *aux_mvt = MVT_SOLID;
          }
        }
      }
    }
  }

  if (rot_map) {
    out_of_space_mvt(1);
  }

  //speed up the stabilization of the flow
  if (lgca_speedup) {
    LogPrintf("speedup of the lattice gas: %d\n", lgca_speedup);
    for (i = 0; i < lgca_speedup; i++) {
      do_force_mvt();
    }
  }

  if (Velx_interp) {
    ResetMemory(Velx_interp, float, HLD);
    ResetMemory(Vely_interp, float, HLD);
  }
}

// set the lattice gas out of the rotating space, using periodic boundary conditions
void out_of_space_mvt(int32_t reset_mvt) {
  int32_t i, j, k, ic, ck, typ;
  MvtField *aux_mvt;
  Pos2 cp;

  if (reset_mvt) {
    LogPrintf("initialization of the lattice gas outside the rotating space\n");
  }
  for (ck = CLN; ck <= CLS; ck++) { // depth
    k = LN + ck * DIST_MVT_NS;
    for (j = 0; j < H; j++) { // height
      aux_mvt = CelMvt + j * CLEO * NB_MVT_VER + ck * CHL;
      for (i = 0; i < L; i++) { // width
        if (OutOfSpace(i, k)) {
          if (pbc_mode) {
            cp = RotMapPos(i, k);
            typ = CellType(cp.x, j, cp.y);
            for (ic = 0; ic < NB_MVT_EO; ic++, aux_mvt++) {
              if (Phase[typ] == SOLID) {
                SetMask(*aux_mvt, MVT_SOLID);
              } else {
                UnsetMask(*aux_mvt, MVT_SOLID);
                if (reset_mvt) {
                  init_mvt_cell(aux_mvt);
                  nb_cel_fluide++;
                }
              }
            }
          } else {
            // Open boundary conditions: by default we put solid ground and ceiling
            for (ic = 0; ic < NB_MVT_EO; ic++, aux_mvt++) {
              if ((j <= 1) || (j >= H - 2)) {
                SetMask(*aux_mvt, MVT_SOLID);
              } else {
                UnsetMask(*aux_mvt, MVT_SOLID);
                if (reset_mvt) {
                  init_mvt_cell(aux_mvt);
                  nb_cel_fluide++;
                }
              }
            }
          }
        } else {
          aux_mvt += NB_MVT_EO;
        }
      }
    }
  }
}

// modification globale de la terre
void collisions_mod_terre() {
  int32_t i, j, k, ic;
  Cell *aux;
  MvtField *aux_mvt;

  nb_cel_fluide = 0;
  aux = TE + HLN;
  aux_mvt = CelMvt;
  for (k = 0; k < CLNS; k++, aux += HL * (DIST_MVT_NS - 1)) { // depth
    for (j = 0; j < H; j++) { // height
      for (i = 0; i < L; i++, aux++) { // width
        for (ic = 0; ic < NB_MVT_EO; ic++, aux_mvt++) {
          if (Phase[aux->celltype] == SOLID) {
            SetMask(*aux_mvt, MVT_SOLID);
          } else {
            UnsetMask(*aux_mvt, MVT_SOLID);
            nb_cel_fluide++;
          }
        }
      }
    }
  }
}

void lecture_mvt() {
  FILE *fp;
  int32_t ix;

  LogPrintf("lecture fichier HPP : %s\n", nom_fic_mvt);

  fp = fopen(nom_fic_mvt, "r");
  if (! fp) {
    ErrPrintf("ERROR: cannot open HPP file %s\n", nom_fic_mvt);
    exit(-4);
  }

  fread(CelMvt, sizeof(MvtField), CHLD, fp);

  fclose(fp);

  // initialisation nb_cel_fluide
  for (ix = 0; ix < CHLD; ix++) {
    if (!(CelMvt[ix] & MVT_SOLID)) {
      nb_cel_fluide++;
    }
  }

  // moyennage de la vitesse initiale
  compute_vel(0);
}


int32_t mvt(int32_t index) {
  return (CelMvt[index] & MVT_ALLDIR);
}


int32_t mvt_solid(int32_t ii) {
  return (CelMvt[ii] & MVT_SOLID);
}

int32_t mvt_bas(int32_t ii) {
  return (CelMvt[ii] & (MVT_B | MVT_EB | MVT_OB));
}

int32_t mvt_haut(int32_t ii) {
  return (CelMvt[ii] & (MVT_H | MVT_EH | MVT_OH));
}

int32_t mvt_est(int32_t ii) {
  return (CelMvt[ii] & (MVT_E | MVT_EB | MVT_EH));
}

int32_t mvt_ouest(int32_t ii) {
  return (CelMvt[ii] & (MVT_O | MVT_OB | MVT_OH));
}

int32_t check_mvt(int32_t index) {
  int32_t ii;
  Calcule_cix(index, ii);
  return (CelMvt[ii] & MVT_ALLDIR);
}

int32_t check_mvt_solid(int32_t index) {
  int32_t ii;
  Calcule_cix(index, ii);
  return (CelMvt[ii] & MVT_SOLID);
}

int32_t check_mvt_bas(int32_t index) {
  int32_t ii;
  Calcule_cix(index, ii);
  return (CelMvt[ii] & (MVT_B | MVT_EB | MVT_OB));
}

int32_t check_mvt_haut(int32_t index) {
  int32_t ii;
  Calcule_cix(index, ii);
  return (CelMvt[ii] & (MVT_H | MVT_EH | MVT_OH));
}

int32_t check_mvt_est(int32_t index) {
  int32_t ii;
  Calcule_cix(index, ii);
  return (CelMvt[ii] & (MVT_E | MVT_EB | MVT_EH));
}

int32_t check_mvt_ouest(int32_t index) {
  int32_t ii;
  Calcule_cix(index, ii);
  return (CelMvt[ii] & (MVT_O | MVT_OB | MVT_OH));
}

int32_t check_mvt_est_et_bas(int32_t index) {
  return check_mvt_est(index) && check_mvt_bas(index);
}

int32_t check_mvt_EB(int32_t index) {
  int32_t ii;
  Calcule_cix(index, ii);
  return (CelMvt[ii - 10 * CLEO] & (MVT_B | MVT_E | MVT_EB)); //maneuvering
}

int32_t check_no_mvt_est_et_bas(int32_t index) {
  return !check_mvt_est_et_bas(index);
}

int32_t check_mvt_mean_bas(int32_t index) {
#if (VSTEP > 1)
  int32_t cx, cy, cz;
  Calcule_cxyz(index, cx, cy, cz);
  return (Vely[(cx / VSTEP_L) + (cy / VSTEP_H) * VL + cz * VHL] > 0);
#else
  return (Vely[index] > 0);
#endif
}

int32_t check_mvt_mean_est(int32_t index) {
#if (VSTEP > 1)
  int32_t cx, cy, cz;
  Calcule_cxyz(index, cx, cy, cz);
  return (Velx[(cx / VSTEP_L) + (cy / VSTEP_H) * VL + cz * VHL] > 0);
#else
  return (Velx[index] > 0);
#endif
}

int32_t check_no_mvt_mean_bas(int32_t index) {
  return !check_mvt_mean_bas(index);
}

void collisions_modcell(int32_t type, int32_t index) {
  int32_t x, y, z, ii, ic, cz;
  Calcule_xyz(index, x, y, z);
  cz = (z - LN) / DIST_MVT_NS;
  if (z != cz * DIST_MVT_NS + LN) {
    return;
  }
  ii = x * NB_MVT_EO + y * CLEO + cz * CHL;
  for (ic = 0; ic < NB_MVT_EO; ic++, ii++) {
    if (Phase[type] == SOLID) {
      if (!(CelMvt[ii] & MVT_SOLID)) {
        nb_cel_fluide--;
      }
      SetMask(CelMvt[ii], MVT_SOLID);
    } else {
      if (CelMvt[ii] & MVT_SOLID) {
        nb_cel_fluide++;
      }
      UnsetMask(CelMvt[ii], MVT_SOLID);
    }
  }
}


void do_collisions() {
  int32_t i, j, k, ii;
  for (k = CLN; k <= CLS; k++) { // depth
    ii = CLEO * CH * k + CLEO;
    for (j = 1; j < CH - 1; j++) { // height
      ii += NB_MVT_EO;
      for (i = NB_MVT_EO; i < CLEO - NB_MVT_EO; i++, ii++) { // width
        CelMvt2[ii] = Collisions[CelMvt[ii]];
        if ((j <= h_ceil) && (NbMvtCell[CelMvt[ii] & MVT_FAST] == 1)) {
          //rebond realiste au plafond
          MvtField mvt2 = 0;
          switch (CelMvt[ii] & MVT_FAST) {
          case MVT_EH:
            mvt2 = MVT_EB;
            break;
          case MVT_OH:
            mvt2 = MVT_OB;
            break;
          }
          if (mvt2) {
            UnsetMask(CelMvt2[ii], MVT_FAST);
            SetMask(CelMvt2[ii], mvt2);
          }
        }
      }
      ii += NB_MVT_EO;
    }
    ii += CLEO;
  }
  col_iter++;
}

void propagate_mvtcel(int32_t i, int32_t ii) {
  int32_t ii0;

  CelMvt[ii] = 0;

  //SOLID+OUT
  CelMvt[ii] |= (CelMvt2[ii] & (MVT_SOLID | MVT_OUT));

  //EST
  ii0 = ii - 1;
#ifdef CYCLAGE_MVT
  if (i == NB_MVT_EO) {
    ii0 += CLEO - 2 * NB_MVT_EO;
  }
#endif
  CelMvt[ii] |= (CelMvt2[ii0] & MVT_E);

  //OUEST
  ii0 = ii + 1;
#ifdef CYCLAGE_MVT
  if (i == CLEO - NB_MVT_EO - 1) {
    ii0 -= CLEO - 2 * NB_MVT_EO;
  }
#endif
  CelMvt[ii] |= (CelMvt2[ii0] & MVT_O);

  //BAS
  ii0 = ii - CLEO;
  CelMvt[ii] |= (CelMvt2[ii0] & MVT_B);

  //HAUT
  ii0 = ii + CLEO;
  CelMvt[ii] |= (CelMvt2[ii0] & MVT_H);

  //EST+BAS
  ii0 = ii - CLEO - 1;
#ifdef CYCLAGE_MVT
  if (i == NB_MVT_EO) {
    ii0 += CLEO - 2 * NB_MVT_EO;
  }
#endif
  CelMvt[ii] |= (CelMvt2[ii0] & MVT_EB);

  //EST+HAUT
  ii0 = ii + CLEO - 1;
#ifdef CYCLAGE_MVT
  if (i == NB_MVT_EO) {
    ii0 += CLEO - 2 * NB_MVT_EO;
  }
#endif
  CelMvt[ii] |= (CelMvt2[ii0] & MVT_EH);

  //OUEST+BAS
  ii0 = ii - CLEO + 1;
#ifdef CYCLAGE_MVT
  if (i == CLEO - NB_MVT_EO - 1) {
    ii0 -= CLEO - 2 * NB_MVT_EO;
  }
#endif
  CelMvt[ii] |= (CelMvt2[ii0] & MVT_OB);

  //OUEST+HAUT
  ii0 = ii + CLEO + 1;
#ifdef CYCLAGE_MVT
  if (i == CLEO - NB_MVT_EO - 1) {
    ii0 -= CLEO - 2 * NB_MVT_EO;
  }
#endif
  CelMvt[ii] |= (CelMvt2[ii0] & MVT_OH);

}

void do_propagations() {
  int32_t i, j, k, ii;

  for (k = CLN; k <= CLS; k++) { // profondeur
    for (j = 1; j < CH - 1; j++) { // hauteur
      i = NB_MVT_EO;
      ii = CLEO * CH * k + CLEO * j + i;
      while (i < CLEO - NB_MVT_EO) { //largeur
        propagate_mvtcel(i++, ii++);
      }
    }
  }

  //Retroaction transport/matiere

#ifdef MVT_REGUL
  do_regul_mvt();
#endif
#if defined(MVT_HFOR) || defined(MVT_VFOR) || defined(MVT_GFOR) //forcage flux
  do_force_mvt();
#endif

  if (lgca_mixing) {
    do_mixing();
  }

#ifdef DUMP_SIGNATURE
  //empreinte gazeuse
  if (opt_info && !(col_iter & 0xf)) {
    dump_signature_mvt();  //1 cycle sur 16
  }
#endif

  if (opt_info) {
    dump_mvt_in_out();
  }
}

void do_mvt_in(int32_t ix) {
  float aleat = drand48();
  int32_t k = (int)(ix / CHL);
  //injection d'une particule lente
  if (aleat < 0.5 * DENSITE)
    if (!(CelMvt[ix] & (MVT_SOLID | MVT_OUT | MVT_E | MVT_O))) {
      CelMvt[ix] |= MVT_E;
      NbMvt[k]++;
      nb_mvt_in++;
    }

  //injection d'une particule rapide
  if (aleat < 0.25 * DENSITE) {
    if (!(CelMvt[ix] & (MVT_SOLID | MVT_OUT | MVT_EB | MVT_OH))) {
      CelMvt[ix] |= MVT_EB;
      NbMvt[k]++;
      nb_mvt_in++;
    }
  } else if (aleat < 0.5 * DENSITE) {
    if (!(CelMvt[ix] & (MVT_SOLID | MVT_OUT | MVT_EH | MVT_OB))) {
      CelMvt[ix] |= MVT_EH;
      NbMvt[k]++;
      nb_mvt_in++;
    }
  }
}

#ifdef MVT_REGUL
void do_regul_mvt() {
  int32_t k, ix, y, cpt;

  for (k = CLN; k <= CLS; k++) {
    //on injecte des particules sur le bord gauche
    cpt = 0;
    while ((NbMvt[k] < NbMvt0[k]) && (cpt < CH * 2/*CH/2*/)) {
      y = CH * drand48();
      ix = NB_MVT_EO + y * CLEO + k * CHL;
      do_mvt_in(ix);
      cpt++;
    }

    if (NbMvt[k] < NbMvt0[k] * 0.9) LogPrintf("NbMvt[%d] = %ld\n", k, NbMvt[k])
    }
}
#endif // MVT_REGUL

void do_force_mvt() {
  int32_t ix;
  // injection flux maximum sur le bord gauche
#ifdef MVT_HFOR
  // forcage flux ouest-est sur la couche superieure
  for (z = CLN; z <= CLS; z++) {
    for (x = 0; x < CLEO; x++) {
      ix = x + CLEO * NB_MVT_VER * (h_ceil + 1) + z * CHL;
      if ((CelMvt[ix] & (MVT_E | MVT_O)) == MVT_O) {
        UnsetMask(CelMvt[ix], MVT_O);
        SetMask(CelMvt[ix], MVT_E);
      }
      if ((CelMvt[ix] & (MVT_EH | MVT_OH)) == MVT_OH) {
        UnsetMask(CelMvt[ix], MVT_OH);
        SetMask(CelMvt[ix], MVT_EH);
      }
      if ((CelMvt[ix] & (MVT_EB | MVT_OB)) == MVT_OB) {
        UnsetMask(CelMvt[ix], MVT_OB);
        SetMask(CelMvt[ix], MVT_EB);
      }
    }
  }
#endif // MVT_HFOR

#ifdef MVT_VFOR
  // forcage partiel du flux sur le bord gauche
  int32_t cpt_vfor = 0;

  for (z = CLN; z <= CLS; z++) {
    for (x = 0; x < 1; x++) {
      ix = z * CHL + 5 * CLEO + NB_MVT_EO + x;
      for (y = 0; y < CH; y++, ix += CLEO) {
        //forcage
        if (drand48() < intensite_vfor) {
          if ((CelMvt[ix] & (MVT_E | MVT_O)) == MVT_O) {
            UnsetMask(CelMvt[ix], MVT_O);
            SetMask(CelMvt[ix], MVT_E);
            cpt_vfor++;
          }

          if ((CelMvt[ix] & (MVT_EH | MVT_OH)) == MVT_OH) {
            UnsetMask(CelMvt[ix], MVT_OH);
            SetMask(CelMvt[ix], MVT_EH);
            cpt_vfor++;
          }

          if ((CelMvt[ix] & (MVT_EB | MVT_OB)) == MVT_OB) {
            UnsetMask(CelMvt[ix], MVT_OB);
            SetMask(CelMvt[ix], MVT_EB);
            cpt_vfor++;
          }
        }
      }
    }
  }
#endif //MVT_VFOR

#ifdef MVT_GFOR
  // forcage global du flux
  int32_t cpt_for = 0;
  int32_t nb_mvt_for = nb_cel_fluide * intensite_gfor;

  while (nb_mvt_for > 0) {
    ix = CLN * CHL + drand48() * ((CLS - CLN + 1) * CHL);
    if (CelMvt[ix] & MVT_SOLID) {
      continue;
    }
    if ((CelMvt[ix] & (MVT_E | MVT_O)) == MVT_O) {
      UnsetMask(CelMvt[ix], MVT_O);
      SetMask(CelMvt[ix], MVT_E);
      cpt_for++;
    }

    if ((CelMvt[ix] & (MVT_EH | MVT_OH)) == MVT_OH) {
      UnsetMask(CelMvt[ix], MVT_OH);
      SetMask(CelMvt[ix], MVT_EH);
      cpt_for++;
    }

    if ((CelMvt[ix] & (MVT_EB | MVT_OB)) == MVT_OB) {
      UnsetMask(CelMvt[ix], MVT_OB);
      SetMask(CelMvt[ix], MVT_EB);
      cpt_for++;
    }
    nb_mvt_for--;
  }
#endif // MVT_GFOR
}

void do_mixing() {
  int32_t i, y, z, y1, y2;
  int32_t cix, cix1, cix2, mvt;

  for (z = CLN; z <= CLS; z++) {
    for (i = 0; i < (CH / 2); i++) {

      //brouillage directionnel (bord gauche)
      y = 5 + (CH - 10) * drand48();
      cix = z * CHL + y * CLEO + NB_MVT_EO;
      mvt = CelMvt[cix] & (MVT_EB | MVT_EH);
      if (mvt == MVT_EB) {
        UnsetMask(CelMvt[cix], MVT_EB);
        SetMask(CelMvt[cix], MVT_EH);
      } else if (mvt == MVT_EH) {
        UnsetMask(CelMvt[cix], MVT_EH);
        SetMask(CelMvt[cix], MVT_EB);
      }

      //brouillage positionnel (bord gauche)
      y1 = 5 + (CH - 10) * drand48();
      y2 = 5 + (CH - 10) * drand48();
      cix1 = z * CHL + y1 * CLEO + NB_MVT_EO;
      cix2 = z * CHL + y2 * CLEO + NB_MVT_EO;
      if (!mvt_solid(cix1) && !mvt_solid(cix2)) {
        //permutation
        mvt = CelMvt[cix1];
        CelMvt[cix1] = CelMvt[cix2];
        CelMvt[cix2] = mvt;
      }
    }
  }
}
//#endif


void compute_vel(int8_t flag_interp) {
  static int32_t cpt = 0;
  int32_t i, j, k, ix, cpt1;
  int32_t *pvx, *pvy;
  int32_t vel_size;
  float fluid_ratio;

#ifndef VEL_SLIDE
  vel_size = VHLD;
#else
  vel_size = CHLD;
#endif

  if (!Velx) {
    AllocMemoryPrint("Velx", Velx, int, vel_size);
    ResetMemory(Velx, int, vel_size);
  }
  if (!Vely) {
    AllocMemoryPrint("Vely", Vely, int, vel_size);
    ResetMemory(Vely, int, vel_size);
  }
  if (!Velx_time) {
    if (VSTEP_TIME > 1) {
      AllocMemoryPrint("Velx_time", Velx_time, int, vel_size * VSTEP_TIME);
      ResetMemory(Velx_time, int, vel_size * VSTEP_TIME);
    } else {
      Velx_time = Velx;
    }
  }
  if (!Vely_time) {
    if (VSTEP_TIME > 1) {
      AllocMemoryPrint("Vely_time", Vely_time, int, vel_size * VSTEP_TIME);
      ResetMemory(Vely_time, int, vel_size * VSTEP_TIME);
      PrintTotalMemory();
    } else {
      Vely_time = Vely;
    }
  }

#ifdef VEL_SLIDE
  /// sliding window
  if (!Veln) {
    AllocMemoryPrint("Veln", Veln, int, vel_size);
    ResetMemory(Veln, int, vel_size);
  }
  /// sliding vectors
  if (!Velx_vector) {
    AllocMemoryPrint("Velx_vector", Velx_vector, int, CHLD);
    ResetMemory(Velx_vector, int, CHLD);
  }
  if (!Vely_vector) {
    AllocMemoryPrint("Vely_vector", Vely_vector, int, CHLD);
    ResetMemory(Vely_vector, int, CHLD);
  }
  if (!Veln_vector) {
    AllocMemoryPrint("Veln_vector", Veln_vector, int, CHLD);
    ResetMemory(Veln_vector, int, CHLD);
  }
#endif

#ifndef VEL_SLIDE
  /// averaging over neighboring space by using fixed windows
  int32_t iv, jv, imax, jmax;
  int32_t nb_mvt_fluide;

  for (k = CLN; k <= CLS; k++) { //profondeur
    pvx = Velx_time + VHLD * cpt + VHL * k;
    pvy = Vely_time + VHLD * cpt + VHL * k;
    for (jv = 0; jv < VH; jv++) { //hauteur
      for (iv = 0; iv < VL; iv++, pvx++, pvy++) { //largeur
        // moyenne des vitesses dans un rectangle de taille VSTEP_H x VSTEP_L
        ix = k * CHL + jv * CLEO * VSTEP_H + iv * VSTEP_L;
        if (jv < VH - 1) {
          jmax = VSTEP_H;
        } else {
          jmax = CH - jv * VSTEP_H;
        }
        if (iv < VL - 1) {
          imax = VSTEP_L;
        } else {
          imax = CLEO - iv * VSTEP_L;
        }
        *pvx = *pvy = 0;
        nb_mvt_fluide = 0;
        for (j = 0; j < jmax; j++, ix += CLEO - imax) { //hauteur
          for (i = 0; i < imax; i++, ix++) { //largeur
            *pvx += CelVelx[CelMvt[ix]];
            *pvy += CelVely[CelMvt[ix]];

            if (!(CelMvt[ix] & MVT_SOLID)) {
              nb_mvt_fluide++;
            }
          }
        }
        //moyennage sur les cellules fluides uniquement
        fluid_ratio = (float)nb_mvt_fluide / (imax * jmax);
        if (fluid_ratio < 0.2) {
          fluid_ratio = 0.2;
        }
        *pvx = (*pvx) / fluid_ratio;
        *pvy = (*pvy) / fluid_ratio;
      }
    }
  }
#else // if def VEL_SLIDE
  /// averaging over neighboring space by using sliding windows

  int32_t *pwx, *pwy, *pwn, *pvn;
  int32_t *pwx0, *pwy0, *pwx1, *pwy1;
  int32_t iw, jw;
  int32_t i0 = VSTEP_L / 2;
  int32_t j0 = VSTEP_H / 2;
  int32_t i1 = VSTEP_L - i0;

  pvx = pvy = pvn = NULL;

  for (k = CLN; k <= CLS; k++) { //depth
    pvx = Velx_vector + CHL * k;
    pvy = Vely_vector + CHL * k;
    pvn = Veln_vector + CHL * k;
    for (j = 0; j < CH - VSTEP_H; j++) { //vertical slide
      //sliding window of size (VSTEP_H x VSTEP_L)
      pwx = Velx_time + CHLD * cpt + CHL * k + CLEO * (j + j0) + i0;
      pwy = Vely_time + CHLD * cpt + CHL * k + CLEO * (j + j0) + i0;
      pwn = Veln + CHL * k + CLEO * (j + j0) + i0;
      for (i = 0; i < CLEO - VSTEP_L; i++, pvx++, pvy++, pvn++, pwx++, pwy++, pwn++) { //horizontal slide
        //averaging on the sliding window
        if (j == 0) {
          if (i == 0) {
            //upper left window of the plane
            ix = k * CHL;
            *pwx = *pwy = *pwn = 0;
            for (jw = 0; jw < VSTEP_H; jw++, ix += CLEO - VSTEP_L, pvx += CLEO, pvy += CLEO, pvn += CLEO) { //height
              //sliding vectors of size VSTEP_L
              *pvx = *pvy = *pvn = 0;
              for (iw = 0; iw < VSTEP_L; iw++, ix++) { //width
                if (CelMvt[ix] & MVT_SOLID) {
                  continue;
                }
                *pvx += CelVelx[CelMvt[ix]];
                *pvy += CelVely[CelMvt[ix]];
                (*pvn)++;
              }
              *pwx += *pvx;
              *pwy += *pvy;
              *pwn += *pvn;
            }
          } else {
            //upper window
            ix = k * CHL + (VSTEP_L - 1 + i);
            *pwx = *pwy = *pwn = 0;
            for (jw = 0; jw < VSTEP_H; jw++, ix += CLEO) { //height
              pvx = Velx_vector + CHL * k + CLEO * jw + i;
              pvy = Vely_vector + CHL * k + CLEO * jw + i;
              pvn = Veln_vector + CHL * k + CLEO * jw + i;

              //slide vectors to the right
              *pvx = *(pvx - 1) - CelVelx[CelMvt[ix - VSTEP_L]] + CelVelx[CelMvt[ix]];
              *pvy = *(pvy - 1) - CelVely[CelMvt[ix - VSTEP_L]] + CelVely[CelMvt[ix]];
              *pvn = *(pvn - 1) + (CelMvt[ix - VSTEP_L] & MVT_SOLID) - (CelMvt[ix] & MVT_SOLID);
              *pwx += *pvx;
              *pwy += *pvy;
              *pwn += *pvn;
            }
          }
        } else {
          if (i == 0) {
            //left window
            ix = k * CHL + (VSTEP_H - 1 + j) * CLEO;

            //new sliding vectors
            pvx = Velx_vector + CHL * k + CLEO * (VSTEP_H - 1 + j);
            pvy = Vely_vector + CHL * k + CLEO * (VSTEP_H - 1 + j);
            pvn = Veln_vector + CHL * k + CLEO * (VSTEP_H - 1 + j);
            *pvx = *pvy = *pvn = 0;

            for (iw = 0; iw < VSTEP_L; iw++, ix++) { //width
              if (CelMvt[ix] & MVT_SOLID) {
                continue;
              }
              *pvx += CelVelx[CelMvt[ix]];
              *pvy += CelVely[CelMvt[ix]];
              (*pvn)++;
            }
          } else {
            //sliding window
            ix = k * CHL + (VSTEP_H - 1 + j) * CLEO + (VSTEP_L - 1 + i);

            //slide vectors to the right
            *pvx = *(pvx - 1) - CelVelx[CelMvt[ix - VSTEP_L]] + CelVelx[CelMvt[ix]];
            *pvy = *(pvy - 1) - CelVely[CelMvt[ix - VSTEP_L]] + CelVely[CelMvt[ix]];
            *pvn = *(pvn - 1) + (CelMvt[ix - VSTEP_L] & MVT_SOLID) - (CelMvt[ix] & MVT_SOLID);
          }

          //slide down the window
          *pwx = *(pwx - CLEO) - *(pvx - CLEO * VSTEP_H) + *pvx;
          *pwy = *(pwy - CLEO) - *(pvy - CLEO * VSTEP_H) + *pvy;
          *pwn = *(pwn - CLEO) - *(pvn - CLEO * VSTEP_H) + *pvn;

          //full computation
        }
        if (*pwn > VSTEP_L * VSTEP_H) {
          ErrPrintf("ERROR: compute_vel - i=%d  j=%d  k=%d  pwn=%d\n", i, j, k, *pwn);
          exit(-3);
        }
      }
    }
  }

  /// normalize with the number of fluid cells in sliding windows
  for (k = CLN; k <= CLS; k++) {
    for (j = 0; j < CH - VSTEP_H; j++) {
      pwx = Velx_time + CHLD * cpt + CHL * k + CLEO * (j + j0) + i0;
      pwy = Vely_time + CHLD * cpt + CHL * k + CLEO * (j + j0) + i0;
      pwn = Veln + CHL * k + CLEO * (j + j0) + i0;
      for (i = 0; i < CLEO - VSTEP_L; i++, pwx++, pwy++, pwn++) {
        fluid_ratio = (float)(*pwn) / (VSTEP_L * VSTEP_H);
        if (fluid_ratio < 0.2) {
          fluid_ratio = 0.2;
        }
        if (fluid_ratio < 1.0) {
          *pwx = (*pwx) / fluid_ratio;
          *pwy = (*pwy) / fluid_ratio;
        }
      }
    }
  }

  /// extend the values up to the boundaries
  /// east-west
  for (k = CLN; k <= CLS; k++) {
    for (j = 0; j < CH - VSTEP_H; j++) {
      //left boundary
      pwx = Velx_time + CHLD * cpt + CHL * k + CLEO * (j + j0);
      pwy = Vely_time + CHLD * cpt + CHL * k + CLEO * (j + j0);
      pwx0 = pwx + i0;
      pwy0 = pwy + i0;
      for (i = 0; i < i0; i++) {
        *(pwx + i) = *pwx0;
        *(pwy + i) = *pwy0;
      }
      //rigth boundary
      pwx += CLEO - 1;
      pwy += CLEO - 1;
      pwx1 = pwx - i1;
      pwy1 = pwy - i1;
      for (i = 0; i < i1; i++) {
        *(pwx - i) = *pwx1;
        *(pwy - i) = *pwy1;
      }
    }
  }
  /// top
  for (k = CLN; k <= CLS; k++) {
    for (i = 0; i < CLEO; i++) {
      pwx = Velx_time + CHLD * cpt + CHL * k + i;
      pwy = Vely_time + CHLD * cpt + CHL * k + i;
      pwx0 = pwx + CLEO * j0;
      pwy0 = pwy + CLEO * j0;
      for (j = 0; j < j0; j++) {
        *(pwx + j * CLEO) = *pwx0;
        *(pwy + j * CLEO) = *pwy0;
      }
    }
  }
#endif //VEL_SLIDE


  /// averaging over time by using a circular buffer
  cpt1 = cpt + 1;
  if (cpt1 == VSTEP_TIME) {
    cpt1 = 0;
  }
  if (VSTEP_TIME > 1) {
    for (i = 0; i < vel_size; i++) {
      Velx[i] += Velx_time[cpt * vel_size + i] - Velx_time[cpt1 * vel_size + i];
      Vely[i] += Vely_time[cpt * vel_size + i] - Vely_time[cpt1 * vel_size + i];
    }
  }

  /// computing mean and max velocity values over all space
  float vel_max = 0;
  int32_t vel_size_2d = vel_size / CLNS;
  int32_t vx, vy;
  float vel, vel_sum;
  vel_sum = 0;
  maxvel = 0;
  for (k = CLN; k <= CLS; k++) {
    for (i = 0; i < vel_size_2d; i++) {
      vx = Velx[i + k * vel_size_2d];
      vy = Vely[i + k * vel_size_2d];
      vel = sqrt(vx * vx + vy * vy);
      vel_sum += vel;
      if (vel_max < vel) {
        vel_max = vel;
      }
    }
    if (maxvel < vel_max) {
      maxvel = vel_max;
    }
  }
  meanvel = vel_sum / (vel_size_2d * (CLS - CLN + 1)); //VHLD;

  if (flag_interp) {
    compute_vel_interp();
  }

  cpt = cpt1;
}

/// interpolation of the velocity
void compute_vel_interp() {
  int32_t i, j, k;
  int32_t ix, ix0, ix1;

  if (!Velx_interp) {
    AllocMemoryPrint("Velx_interp", Velx_interp, float, HLD);
    ResetMemory(Velx_interp, float, HLD);
  }
  if (!Vely_interp) {
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
  int32_t imax = VSTEP / NB_MVT_EO;

  /// interpolation in each vertical plane
  //LogPrintf("interpolation dans chaque plan de collisions\n");
  for (k = CLN; k <= CLS; k++) { //profondeur
    for (jv = 0; jv < VH - 1; jv++) { //hauteur
      pvx = Velx + k * VHL + jv * VL;
      pvy = Vely + k * VHL + jv * VL;
      Velx00 = Velx01 = Velx10 = Velx11 = 0;
      Vely00 = Vely01 = Vely10 = Vely11 = 0;
      for (iv = 0; iv < VL - 2; iv++, pvx++, pvy++) { //largeur
        // interpolation des vitesses dans un rectangle de taille (imax)xVSTEP
        ix = (LN + k * DIST_MVT_NS) * HL + jv * L * VSTEP + iv * imax; // + (VSTEP>>1)*L + (VSTEP>>1);
        Velx00 = *pvx;
        Velx10 = *(pvx + VL);
        Vely00 = *pvy;
        Vely10 = *(pvy + VL);
        if (imax > 1) {
          Velx01 = *(pvx + 1);
          Velx11 = *(pvx + VL + 1);
          Vely01 = *(pvy + 1);
          Vely11 = *(pvy + VL + 1);
        }
        for (j = 0; j < VSTEP; j++, ix += L - imax) { //hauteur
          if (imax > 1) { //interpolation entre les 4 sommets du rectangle englobant
            for (i = 0; i < imax; i++, ix++) { //largeur
              a00 = (VSTEP - j) * (imax - i);
              a01 = (VSTEP - j) * i;
              a10 = j * (imax - i);
              a11 = j * i;
              Velx_interp[ix] = (float)(Velx00 * a00 + Velx01 * a01 + Velx10 * a10 + Velx11 * a11) / (VSTEP * imax);
              Vely_interp[ix] = (float)(Vely00 * a00 + Vely01 * a01 + Vely10 * a10 + Vely11 * a11) / (VSTEP * imax);
            }
          } else { //interpolation verticale entre 2 sommets
            Velx_interp[ix] = (float)(Velx00 * (VSTEP - j) + Velx10 * j) / VSTEP;
            Vely_interp[ix] = (float)(Vely00 * (VSTEP - j) + Vely10 * j) / VSTEP;
            ix += imax;
          }
        }
      }

      // interpolation a cote du bord droit
      // on repete la derniere valeur interpolee
      // c'est plus facile et on reduit les effets de bord ...
      ix = (LN + k * DIST_MVT_NS) * HL + jv * L * VSTEP + iv * imax; // + (VSTEP>>1)*L + (VSTEP>>1);
      for (j = 0; j < VSTEP; j++, ix += L) { //hauteur
        for (i = 0; i < L - iv * imax; i++) { //largeur
          Velx_interp[ix + i] = Velx_interp[ix - 1];
          Vely_interp[ix + i] = Vely_interp[ix - 1];
        }
      }
    }
  }

#else //VEL_SLIDE is defined

  /// no interpolation in the vertical planes : we simply copy Velx and Vely into Velx_interp and Vely_interp
  for (k = CLN; k <= CLS; k++) { //profondeur
    ix = (LN + k * DIST_MVT_NS) * HL;
    ix0 = k * HL;
    for (i = 0; i < HL; i++, ix++, ix0++) {
      Velx_interp[ix] = Velx[ix0];
      Vely_interp[ix] = Vely[ix0];
    }
  }

#endif // VEL_SLIDE

#if (DIST_MVT_NS > 1)
  float fz0, fz1;
  int32_t kv, kvmax;

  /// north-south interpolation between all vertical planes
  //LogPrintf("interpolation nord-sud dans les plans intermediaires\n");
  for (k = CLN; k <= CLS; k++) { /// vertical planes (lattice gas)
    kvmax = (k < CLS) ? DIST_MVT_NS : LNS - CLS * DIST_MVT_NS;
    for (kv = 1; kv < kvmax; kv++) { /// intermediate planes (interpolation)
      ix = (LN + k * DIST_MVT_NS + kv) * HL;
      ix0 = (LN + k * DIST_MVT_NS) * HL;
      ix1 = (k < CLS) ? (LN + (k + 1) * DIST_MVT_NS) * HL : pbc_mode ? LN * HL : ix0;
      fz1 = (float)kv / kvmax;
      fz0 = 1.0 - fz1;
      for (j = 0; j < H; j++) {
        for (i = 0; i < L; i++, ix++, ix0++, ix1++) {
          if (rot_map && OutOfSpace(i, LN + k * DIST_MVT_NS + kv)) {
            continue;
          }
          Velx_interp[ix] = Velx_interp[ix0] * fz0 + Velx_interp[ix1] * fz1;
          Vely_interp[ix] = Vely_interp[ix0] * fz0 + Vely_interp[ix1] * fz1;
        }
      }
    }
  }
#endif // DIST_MVT_NS
}


void dump_mvt_in_out() {
  static int32_t cpt = 0;
  static int32_t step = 0;
  char current_output[128];

  if (!cpt && !step){
    output_write("MVT_IO","      \t mvt_in    \t mvt_out   \t ratio\n");
  }

  if (step >= VSTEP_TIME){
    sprintf(current_output, "%04d: \t%09" PRId64" \t%09" PRId64 " \t%f\n", cpt, nb_mvt_in, nb_mvt_out, (nb_mvt_in) ? ((float)nb_mvt_out/nb_mvt_in) : 0);
    output_write("MVT_IO", current_output);
    step = nb_mvt_in = nb_mvt_out = 0;	//reset
    cpt++;
  } else {
    step++;
  }
}

void dump_densite(){
// Write numbers of moving cells, in solid and fluid, to DENSITE.log
  static int32_t cpt = 0; // count number of times function has been called == number of lines written to file
  int32_t ix, n;
  char output[256];

  //recalcul nb_mvt, nb_mvt_sol et nb_cel_fluide
  int32_t nb_mvt = 0; //nombre de particules mobiles // number of mobile particles
  int32_t nb_mvt_sol = 0; //nombre de particules mobiles localisees dans une cellule solide
			  //number of mobile particles inside a solid cell
  for (ix = CHL * CLN; ix < CHL * (CLS + 1); ix++) {
    n = NbMvtCell[CelMvt[ix]];
    nb_mvt += n;
    if ((CelMvt[ix] & MVT_SOLID)) {
      nb_mvt_sol += n;
    }
  }
  // Calculate density
  densite = (float)(nb_mvt - nb_mvt_sol) / nb_cel_fluide;

  sprintf(output, "%04d: \t%09d \t%09" PRId64 " \t%f \t%09ld\n", cpt++, nb_mvt, nb_cel_fluide, densite, (long)nb_mvt_sol);
  output_write("DENSITE", output);
}

void dump_vel() {
  static int32_t cpt = 0;
  char output[128];

  if ((!Velx_interp) || (!Vely_interp)) {
    cpt++;
    return;
  }

  //fichier VEL.log
  sprintf(output,"%04d: \t  %04d  \t  %04d\n", cpt++, (int)maxvel, (int)meanvel);
  output_write("VEL", output);
}

void dump_mvt(int32_t cpt, int32_t unit) {
  static char filename[100];
  FILE *fp;

  // binary  file
  if (unit == UNIT_COMP) {
    sprintf(filename, "HPP%04d.bin", cpt);
  } else {
    sprintf(filename, "HPP%05d_t0.bin", cpt);
  }

  fp = fopen(filename, "w");
  if (! fp) {
    ErrPrintf("ERROR: cannot open binary file %s\n", filename);
    exit(-4);
  }

  fwrite(CelMvt, sizeof(MvtField), CHLD, fp);

  fclose(fp);

  if (Velx && Vely) {
    if (unit == UNIT_COMP) {
      sprintf(filename, "HPP%04d_%d_%d_%d.vel", cpt, CLNS, VH, VL);
    } else {
      sprintf(filename, "HPP%05d_%d_%d_%d.vel", cpt, CLNS, VH, VL);
    }

    fp = fopen(filename, "w");
    if (! fp) {
      ErrPrintf("ERROR: cannot open VEL file %s\n", filename);
      exit(-4);
    }

    fwrite(Velx, sizeof(int), VHLD, fp);
    fwrite(Vely, sizeof(int), VHLD, fp);

    fclose(fp);
  }
}

void dump_collisions()
{
  int32_t mvt_in, mvt_out;
  char current_output[512]; //must hold sprintf strings of 16 ints

  //dump des regles de collision entre particules fluides
  //dump the rules for collisions between fluid particles
  current_output[0] = '\0';
  strcat(current_output, "          IN            |          OUT\n");
  strcat(current_output, "E  W  B  T EB ET WB WT  |  E  W  B  T EB ET WB WT\n");
  strcat(current_output, "-------------------------------------------------\n");
  output_write("LGCA", current_output);

  for(mvt_in=0; mvt_in<SIZE_MVT_FIELD; mvt_in+=4)
  {
    mvt_out = Collisions[mvt_in];
    sprintf(current_output, "%d  %d  %d  %d  %d  %d  %d  %d  |  %d  %d  %d  %d  %d  %d  %d  %d\n",
      (mvt_in & MVT_E)? 1 : 0, (mvt_in & MVT_O)? 1 : 0, (mvt_in & MVT_B)? 1 : 0, (mvt_in & MVT_H)? 1 : 0,
      (mvt_in & MVT_EB)? 1 : 0, (mvt_in & MVT_EH)? 1 : 0, (mvt_in & MVT_OB)? 1 : 0, (mvt_in & MVT_OH)? 1 : 0,
      (mvt_out & MVT_E)? 1 : 0, (mvt_out & MVT_O)? 1 : 0, (mvt_out & MVT_B)? 1 : 0, (mvt_out & MVT_H)? 1 : 0,
      (mvt_out & MVT_EB)? 1 : 0, (mvt_out & MVT_EH)? 1 : 0, (mvt_out & MVT_OB)? 1 : 0, (mvt_out & MVT_OH)? 1 : 0);
    output_write("LGCA", current_output);
  }
/*
#ifdef SOLSLOW
  //dump des regles de collision solide-fluide pour les particules lentes
  //selon la configuration du voisinage
  int32_t slow[4] = {MVT_E, MVT_O, MVT_B, MVT_H};
  int32_t mvt_sol_E, mvt_sol_O, mvt_sol_B, mvt_sol_H, mvt_sol, i;
  current_output[0] = '\0';
  strcat(current_output, "\n   SOLID     |      IN      |      OUT\n");
  strcat(current_output, "E  O  B  H   | E  O  B  H   | E  O  B  H\n");
  output_write("LGCA", current_output);
  for (mvt_sol=0; mvt_sol<16; mvt_sol++)
  {
    mvt_sol_E = (mvt_sol & (MVT_E >> 2))? 1 : 0;
    mvt_sol_O = (mvt_sol & (MVT_O >> 2))? 1 : 0;
    mvt_sol_B = (mvt_sol & (MVT_B >> 2))? 1 : 0;
    mvt_sol_H = (mvt_sol & (MVT_H >> 2))? 1 : 0;
    for (i=0; i<4; i++)
    {
      mvt_out = Collisions_SolSlow[mvt_sol][slow[i]];
      //mvt_in = MVT_SOLID | (1 << (2+i));
      //mvt_out = Collisions[mvt_in];
      sprintf(current_output, "% d  %d  %d  %d  |  %d  %d  %d  %d  |  %d  %d  %d  %d\n",
              mvt_sol_E, mvt_sol_O, mvt_sol_B, mvt_sol_H,
              (i==0)? 1 : 0, (i==1)? 1 : 0, (i==2)? 1 : 0, (i==3)? 1 : 0,
              (mvt_out & MVT_E)? 1 : 0, (mvt_out & MVT_O)? 1 : 0, (mvt_out & MVT_B)? 1 : 0, (mvt_out & MVT_H)? 1 : 0);
      output_write("LGCA", current_output);
    }
  }
#endif
#ifdef SOLFAST
  //dump des regles de collision solide-fluide pour les particules rapides (en diagonale)
  //selon la configuration du voisinage
  int32_t fast[4] = {MVT_EB, MVT_EH, MVT_OB, MVT_OH};
  int32_t mvt_sol_EB, mvt_sol_OB, mvt_sol_EH, mvt_sol_OH;//, mvt_sol, i;
  current_output = '\0';
  strcat(current_output, "\n   SOLID     |      IN      |      OUT\n");
  strcat(current_output, "EB EH OB OH  | EB EH OB OH  | EB EH OB OH\n");
  output_write("LGCA", current_output);
  for (mvt_sol=0; mvt_sol<16; mvt_sol++)
  {
    mvt_sol_EB = (mvt_sol & (MVT_EB >> 6))? 1 : 0;
    mvt_sol_EH = (mvt_sol & (MVT_EH >> 6))? 1 : 0;
    mvt_sol_OB = (mvt_sol & (MVT_OB >> 6))? 1 : 0;
    mvt_sol_OH = (mvt_sol & (MVT_OH >> 6))? 1 : 0;
    for (i=0; i<4; i++)
    {
      mvt_out = Collisions_SolFast[mvt_sol][fast[i]];
      sprintf(current_output, "% d  %d  %d  %d  |  %d  %d  %d  %d  |  %d  %d  %d  %d\n",
              mvt_sol_EB, mvt_sol_EH, mvt_sol_OB, mvt_sol_OH,
              (i==0)? 1 : 0, (i==1)? 1 : 0, (i==2)? 1 : 0, (i==3)? 1 : 0,
              (mvt_out & MVT_EB)? 1 : 0, (mvt_out & MVT_EH)? 1 : 0, (mvt_out & MVT_OB)? 1 : 0, (mvt_out & MVT_OH)? 1 : 0);
      output_write("LGCA", current_output);
    }
  }
#endif
 */
}

#ifdef DUMP_SIGNATURE
void dump_signature_mvt()
{
  char current_output[128];
  size_t i;
  uint32_t sig, *aux;

  //calcul de la signature
  sig = 0;
  aux = (unsigned int*)CelMvt;
  for (i = 0; i < CHLD * sizeof(MvtField) / sizeof(unsigned int); i++, aux++) {
    sig += (*aux);  //sig = sig ^ (*aux);
  }

  //dump
  sprintf(current_output, "%08x : %08x\n", col_iter, sig);
  output_write("SIGN_HPP", current_output);

}
#endif

#endif //LGCA
