/* ReSCAL - 3D space processes
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
#include <math.h>
#include <assert.h>
#include <stdint.h>

#include "defs.h"
#include "macros.h"
#include "format.h"
#include "space.h"
#include "doublets.h"
#include "cells.h"
#include "surface.h"
#include "lgca.h"
#include "transitions.h"
#include "view.h"
#include "simul.h"
#include "trace.h"

extern double csp_time;                  // temps reel simule
extern int32_t pbc_mode;  //periodic boundary conditions
extern int32_t use_lgca;

extern int32_t ava_norm;
extern int32_t rot_mode;

int32_t L = 0, H = 0, D = 0, HL = 0, HLD = 0; // dimensions of the cellular space
int32_t L_bounds = 1, D_bounds = 1; //thickness of lateral boundaries
int32_t LN, LS, LEO, LNS, HLN;    //couloir est-ouest (limite nord, limite sud, largeur nord-sud, ...)
int32_t Hd2, Ld2, Dd2;
Cell  *TE = NULL;          // la 'terre'
char *bin_filename = NULL; //nom du fichier BIN
char *csp_filename = NULL; //nom du fichier CSP
int32_t L_rot = 0, D_rot = 0; // dimensions of the rotating space

#ifdef REFDB_PTR
RefDoublets *RefDB = NULL;     // references des cellules de la terre vers les doublets actifs
#else
RefDoublets_Type *RefDB_Type = NULL;     // references des cellules de la terre vers les doublets actifs
RefDoublets_Ind *RefDB_Ind = NULL;     // references des cellules de la terre vers les doublets actifs
#endif
int32_t graine = 273676;           // graine pour la generation des nombres aleatoires
int32_t h_ceil = 0;               //epaisseur du plafond
int32_t h_floor = 0;               //epaisseur du sol
int32_t nb_pv = 0;                // nombre de plans verticaux (non-DUM)
char *psol = NULL;            // indicateur des plans solides verticaux est-ouest
char *csol = NULL;            // indicateur des colonnes solides
float csp_angle = 0.0;        //angle resultant de toutes les rotations
char *rot_map = NULL;        // periodic mapping of the rotating space
Pos2 *rot_map_pos = NULL;        // periodic mapping of the rotating space
Pos2 *rot_map_pos0 = NULL;        // old periodic mapping of the rotating space

void init_terre() {
}

void lock_csp(int32_t log_flag) {
  (void)log_flag;
#ifdef CSP_MUTEX
#endif
}

void unlock_csp(int32_t log_flag) {
  (void)log_flag; //SUPPRESS: unused warning
#ifdef CSP_MUTEX
#endif
}

void cree_terre() {
  int32_t i, j, k;
  Cell *aux;
  char *filename;

  if (pbc_mode) {
    LogPrintf("cyclage horizontal\n");
  }

#ifdef CELL_TIME
  LogPrintf("datation des cellules (mode 'sedimento')\n");
#endif

  lock_csp(0);

  filename = csp_filename ? csp_filename : bin_filename;
  if (csp_filename) {
    read_csp_header(filename);
  }

#ifdef CYCLAGE_HOR
  if (rot_mode == ROT_MODE_FULL) {
    L_rot = L;
    D_rot = D;
    float diag_length = sqrtf(L * L + D * D);
    LogPrintf("diag_length = %f\n", diag_length);
    int32_t margin = 5;
#ifdef LGCA
    if (use_lgca) {
      margin = DIST_MVT_NS + 1;
    }
#endif // LGCA
    LogPrintf("margin = %d\n", margin);
    L_bounds += (int)ceilf((diag_length - L) / 2) + margin;
    D_bounds += (int)ceilf((diag_length - D) / 2) + margin;
    LogPrintf("L_bounds = %d\n", L_bounds);
    LogPrintf("D_bounds = %d\n", D_bounds);
  }
#endif

  H += 2;
  L += 2 * L_bounds;
  if (!D) {
    D = 1;
  }
  D += 2 * D_bounds;
  HL = H * L;
  HLD = HL * D;
  csp_set_bounds(L_bounds, D_bounds, 1);
  LEO = L - 2;
  LNS = D - 2;

  LogPrintf("size of cellular space with bounds (HxLxD): %d x %d x %d\n", H, L, D);

  AllocMemoryPrint("TE", TE, Cell, HLD);
  ResetMemory(TE, Cell, HLD);

  #ifdef REFDB_PTR
  AllocMemoryPrint("RefDB", RefDB, RefDoublets, HLD);
  #else
  AllocMemoryPrint("RefDB_Type", RefDB_Type, RefDoublets_Type, HLD);
  AllocMemoryPrint("RefDB_Ind", RefDB_Ind, RefDoublets_Ind, HLD);
  #endif
  PrintTotalMemory();

  srand48(graine);
  LogPrintf("random seed: %d\n", graine);

  // le centre de la terre
  Hd2 = (int) H / 2;
  Ld2 = (int) L / 2;
  Dd2 = (int) D / 2;

  aux = TE;

#if defined(MODEL_DUN) || defined(MODEL_SNO)
  for (k = 0; k < D; k++) { // profondeur
    for (j = 0; j < H; j++) { // hauteur
      for (i = 0; i < L; i++, aux++) { //largeur
        if ((i < L_bounds) || (j == 0) || (k < D_bounds) || (i >= L - L_bounds) || (j == H - 1) || (k >= D - D_bounds)) {
          aux->celltype = BORD;
        } else if ((j == 1) || (j == H - 2)) {
          aux->celltype   = DUM;
        } else if (j > 0.9 * H) {
          aux->celltype   = GR;
        } else {
          aux->celltype   = EAUC;
        }
      }
    }
  }
#else
  uint8_t type = (BORD >= 2) ? 0 : 2;
  for (k = 0; k < D; k++) { // profondeur
    for (j = 0; j < H; j++) { // hauteur
      for (i = 0; i < L; i++, aux++) { //largeur
        if ((i == 0) || (j == 0) || (k == 0) || (i == L - 1) || (j == H - 1) || (k == D - 1)) {
          aux->celltype = BORD;
        } else {
          aux->celltype = type + (int) floor(drand48() * 2);
        }
      }
    }
  }
#endif

  #ifdef REFDB_PTR
  memset(RefDB, 0, sizeof(RefDoublets) * HLD);
  #else
  memset(RefDB_Type, -1, sizeof(RefDoublets_Type) * HLD);
  memset(RefDB_Ind, -1, sizeof(RefDB_Ind) * HLD);
  #endif
  if (filename) {
    read_csp(filename);
  }

  //signal_csp_ready();
  unlock_csp(0);

#ifdef CYCLAGE_HOR
  if (rot_mode == ROT_MODE_FULL) {
    rotate_map(0);
  }
#endif
}

void apocalypse() {
  LogPrintf("destruction terre !\n");
#ifdef REFDB_PTR
  if (RefDB) free(RefDB);
#else
  if (RefDB_Type) free(RefDB_Type);
  if (RefDB_Ind) free(RefDB_Ind);
#endif
  if (TE) {
    free(TE);
  }
}


void cellule_terre(uint8_t type, int32_t ix) {
  //changer l'etat d'une cellule de la terre
  elimine_doublet_est(ix);
  elimine_doublet_ouest(ix);
  elimine_doublet_bas(ix);
  elimine_doublet_haut(ix);
  elimine_doublet_sud(ix);
  elimine_doublet_nord(ix);

  modifie_cellule(type, ix);

  ajoute_doublet_est(ix);
  ajoute_doublet_ouest(ix);
  ajoute_doublet_bas(ix);
  ajoute_doublet_haut(ix);
  ajoute_doublet_sud(ix);
  ajoute_doublet_nord(ix);
}

/// Horizontal cycling in case of periodic boundary conditions
int32_t hcycle(int32_t ix) {
  int32_t x, y, z;

  if (TE[ix].celltype != BORD) {
    //WarnPrintf("WARNING: hcycle - it is not a boundary\n");
    return ix;
  }

  Calcule_xyz(ix, x, y, z);

  if (!rot_map) {
    if (x == 0) {
      ix += L - 2;
    } else if (x == L - 1) {
      ix -= L - 2;
    } else if (z == 0) {
      ix += (D - 2) * HL;
    } else if (z == D - 1) {
      ix -= (D - 2) * HL;
    }
  } else {
    Pos2 *cp;
    cp = rot_map_pos + (x + z * L);
    if (rot_map[x + z * L] >= 0) {
      ErrPrintf("ERROR: hcycle - it is not a boundary (%d, %d, %d, %d, %d, %d)\n", x, y, z, rot_map[x + z * L], cp->x, cp->y);
      exit(-1);
    } //debug
    ix = cp->x + y * L + cp->y * HL;
    if (rot_map[cp->x + cp->y * L] < 0) {
      ErrPrintf("ERROR: hcycle - it is a boundary (%d, %d, %d, %d)\n", cp->x, y, cp->y, rot_map[cp->x + cp->y * L]);
      exit(-1);
    } //debug
  }

  return ix;
}

int32_t get_cell_east(int32_t ix) {
  int32_t ix2;
  ix2 = ix + 1;

#ifdef CYCLAGE_HOR
  if (pbc_mode && (TE[ix2].celltype == BORD)) {
    if (!rot_map) {
      ix2 -= (L - 2);
    } else {
      ix2 = hcycle(ix2);
    }
  }
#endif

  return ix2;
}

int32_t get_cell_west(int32_t ix) {
  int32_t ix2;
  ix2 = ix - 1;

#ifdef CYCLAGE_HOR
  if (pbc_mode && (TE[ix2].celltype == BORD)) {
    if (!rot_map) {
      ix2 += (L - 2);
    } else {
      ix2 = hcycle(ix2);
    }
  }
#endif

  return ix2;
}

int32_t get_cell_south(int32_t ix) {
  int32_t ix2;
  ix2 = ix + HL;

#ifdef CYCLAGE_HOR
  if (pbc_mode && (TE[ix2].celltype == BORD)) {
    if (!rot_map) {
      ix2 -= (D - 2) * HL;
    } else {
      ix2 = hcycle(ix2);
    }
  }
#endif

  return ix2;
}

int32_t get_cell_north(int32_t ix) {
  int32_t ix2;
  ix2 = ix - HL;

#ifdef CYCLAGE_HOR
  if (pbc_mode && (TE[ix2].celltype == BORD)) {
    if (!rot_map) {
      ix2 += (D - 2) * HL;
    } else {
      ix2 = hcycle(ix2);
    }
  }
#endif

  return ix2;
}

int32_t get_cell_down(int32_t ix) {
  return ix + L;
}

int32_t get_cell_up(int32_t ix) {
  return ix - L;
}

int32_t get_cell_dir(int32_t ix, char dir) {
  int32_t ix2=0;
  if (dir == EST) {
    ix2 = get_cell_east(ix);
  } else if (dir == OUEST) {
    ix2 = get_cell_west(ix);
  } else if (dir == SUD) {
    ix2 = get_cell_south(ix);
  } else if (dir == NORD) {
    ix2 = get_cell_north(ix);
  } else if (dir == BAS) {
    ix2 = get_cell_down(ix);
  } else if (dir == HAUT) {
    ix2 = get_cell_up(ix);
  } else {
    assert(0); //Panic    
  }
  return ix2;
}

#ifdef CELL_DATA
//permuter les donnees de 2 cellules
void swap_cell_data(int32_t ix, int32_t ix2) {
  Cell aux;
  char flag1 = 1;
  char flag2 = 1;

  //LogPrintf("offset = %d\n", offset);

  aux = TE[ix];
  //LogPrintf("%08x %08x\n", (&TE[ix]), (unsigned char*)(TE+ix)+offset)

  // les cellules BORD ne doivent pas etre modifiees
  if (TE[ix].celltype == BORD) {
    flag1 = 0;
  }
  if (TE[ix2].celltype == BORD) {
    flag2 = 0;
  }

#ifdef IN
  // les cellules IN ne doivent pas etre modifiees
  if (TE[ix].celltype == IN) {
    flag1 = 0;
  }
  if (TE[ix2].celltype == IN) {
    flag2 = 0;
  }
#endif

#ifdef OUT
  // les cellules OUT ne doivent pas etre modifiees
  if (TE[ix].celltype == OUT) {
    flag1 = 0;
  }
  if (TE[ix2].celltype == OUT) {
    flag2 = 0;
  }
#endif

  if (flag1) {
    const celltype_t te_ix_type = TE[ix].celltype; //Back up cell type; we are about to overwrite it.
    TE[ix] = TE[ix2];             //Copy cell data and overwrite type
    TE[ix].celltype = te_ix_type; //Restore cell type
    //The above is a clearer way of doing the following:
    // memcpy((void *)(TE + ix) + CELL_TYPE_SIZE, (void *)(TE + ix2) + CELL_TYPE_SIZE, CELL_DATA_SIZE);  //TE[ix].celldata <- TE[ix2].celldata
  }
  if (flag2) {
    const celltype_t te_ix2_type = TE[ix2].celltype; //Back up cell type; we are about to overwrite it.
    TE[ix2] = aux;                   //Copy cell data and overwrite type
    TE[ix2].celltype = te_ix2_type; //Restore cell type
    //The above is a clearer way of doing the following:
    // memcpy((void *)(TE + ix2) + CELL_TYPE_SIZE, (void *)(&aux) + CELL_TYPE_SIZE, CELL_DATA_SIZE);  //TE[ix2].celldata <- aux.celldata
  }
}

void update_inout_data(int32_t ix, int32_t ix2) {
  (void)ix;  //SUPPRESS: unused warning
  (void)ix2; //SUPPRESS: unused warning
}

//copier les donnees de la cellule ix dans la cellule ix2
void copy_cell_data(int32_t ix, int32_t ix2) {
  const celltype_t te_ix2_type = TE[ix2].celltype; //Back up cell type; we are about to overwrite it.
  TE[ix2] = TE[ix];               //Copy cell data and overwrite type
  TE[ix2].celltype = te_ix2_type; //Restore cell type
  //The above is a clearer way of doing the following:
  // memcpy((void *)(TE + ix2) + CELL_TYPE_SIZE, (void *)(TE + ix) + CELL_TYPE_SIZE, CELL_DATA_SIZE); //TE[ix2].celldata <- aux.celldata
}

#endif  //CELL_DATA

// callback de controle en fonction de l'altitude
int32_t check_alti(int32_t ix, void *data) {
  (void)data; //SUPPRESS: unused warning
  static int32_t start = 1;
  if (start) {
    LogPrintf("controle du flux entrant en fonction de l'altitude\n");
    start = 0;
  }
  //calcul de la position (x,y,z)
  int32_t x, y, z;
  Calcule_xyz(ix, x, y, z);
  (void)x; //SUPPRESS: unused warning

  //calcul probabiliste suivant y
  float alea = drand48();
  float var = 1.0 - (float)y / H; //probabilite plus forte si altitude elevee (y faible)
  int32_t chk_ok = (alea < var);
  return chk_ok;
}

int32_t check_cell_dir(int32_t ix, void *data) {
  DataCheck *pdc = data;
  int32_t dir = pdc->char_data1;
  int32_t cell_type = pdc->char_data2;
  int32_t ix2;

  ix2 = get_cell_dir(ix, dir);

  return TE[ix2].celltype == cell_type;
}

int32_t check_no_cell_dir(int32_t ix, void *data) {
  return !check_cell_dir(ix, data);
}

#if defined(MODEL_DUN) || defined(MODEL_SNO)
// callback de controle en fonction de la cellule a l'ouest
int32_t check_grain_seul(int32_t index, void *data) {
  (void)data; //SUPPRESS: unused check
  return TE[index - 1].celltype == EAUC;
}
#endif

void translation(int32_t dx, int32_t dz) {
  static Cell *auxmap = NULL;
  static float dx_sum = 0;
  static float dz_sum = 0;
  Cell *aux, *aux2;
  int32_t i, j, k, i0, k0;

  if (!dx && !dz) {
    return;
  }

  if (!auxmap) {
    AllocMemoryPrint("auxmap (translation)", auxmap, Cell, HLD);
  }

  /// computing the resultant translation from start by using current orientation
  float c = cos(csp_angle * PI / 180.0);
  float s = sin(csp_angle * PI / 180.0);
  dx_sum += dx * c + dz * s;
  dz_sum += -dx * s + dz * c;
  LogPrintf("translation: ( %d , %d ) - orientation: %.3f - resultant: ( %.3f , %.3f ) - time: %f\n", dx, dz, csp_angle, dx_sum, dz_sum, csp_time);

  aux = TE;
  aux2 = auxmap;

  for (k = 0; k < D; k++) { // profondeur
    for (j = 0; j < H; j++) { // hauteur
      for (i = 0; i < L; i++, aux++, aux2++) { //largeur
        if ((aux->celltype != BORD) && (aux->celltype != DUM) && (aux->celltype != IN) && (aux->celltype != OUT)) {
          i0 = i - dx;
          k0 = k - dz;
          if (pbc_mode) {
            while (i0 < 1) {
              i0 += L - 2;
            }
            while (i0 >= L - 1) {
              i0 -= L - 2;
            }
            while (k0 < 1) {
              k0 += D - 2;
            }
            while (k0 >= L - 1) {
              k0 -= D - 2;
            }
            *aux2 = *(TE + i0 + j * L + k0 * HL);
          } else {
            if ((i0 >= 1) && (i0 < L - 1) && (k0 >= 1) && (k0 < D - 1)) {
              *aux2 = *(TE + i0 + j * L + k0 * HL);
            }
#if defined(MODEL_DUN) || defined(MODEL_SNO)
            else {
              aux2->celltype = EAUC;
            }
#endif
          }
        } else {
          *aux2   = *aux;
        }
      }
    }
  }
  memcpy(TE, auxmap, HLD * sizeof(Cell));

  //recalcul du nombre de cellules de chaque etat
  if (!pbc_mode) {
    init_Ncel();
  }

  //recalcul des tableaux de doublets actifs
  init_db_pos();

#ifdef ALTI
  calcule_alti(ALTI, ALTI_MODE_BAS);
  if (ava_norm) {
    calcule_normales();
  }
#endif

#ifdef LGCA
  //mise-a-jour de la matrice de gaz-sur-reseau
  if (use_lgca) {
    collisions_mod_terre();
  }
#endif
}

/// set a line of points in a horizontal LxD map (for the periodic boundaries of a rotating space)
void set_line(char *hmap, float x1, float y1, float x2, float y2, char c) {
  float sdx, sdy, dxabs, dyabs;
  int32_t i, ix1, iy1, px, py, i2;
  float slope;
  float dx = x2 - x1;
  float dy = y2 - y1;

  assert(hmap);
  dxabs = fabsf(dx);
  dyabs = fabsf(dy);
  sdx = sgn(dx);
  sdy = sgn(dy);
  i = 0;
  ix1 = Round(x1);
  iy1 = Round(y1);
  if (dxabs >= dyabs) {
    slope = (float)dy / (float)dx;
    i2 = Round(dx);
    while (1) {
      py = Round(slope * i);
      hmap[ix1 + i + (iy1 + py)*L] = c;
      if (i == i2) {
        break;
      }
      i += sdx;
    }
  } else {
    slope = (float)dx / (float)dy;
    i2 = Round(dy);
    while (1) {
      px = Round(slope * i);
      hmap[ix1 + px + (iy1 + i)*L] = c;
      if (i == i2) {
        break;
      }
      i += sdy;
    }
  }
}

void verify_rotating_map() {
  int32_t i, j, k;
  int32_t bmin, bmax;
  Pos2 cp, cp2;
  int32_t ix, ix1, ix2;
  for (i = 0; i < L; i++) {
    bmin = -1;
    bmax = -1;
    /// check initialized values
    for (k = 0; k < D; k++) {
      if (rot_map[i + k * L] == 0) {
        ErrPrintf("ERROR: verify_rotating_map - value 0 (%d,%d)\n", i, k);
        exit(-2);
      }
      if (rot_map[i + k * L] == 1) {
        if (bmin == -1) {
          bmin = k;
        }
        if (bmax < k) {
          bmax = k;
        }
      }
    }
    /// check north-south mapping
    if (bmin >= 0) {
      k = bmin - 1;
      LogPrintf("verify_rotating_map: i=%d bmin=%d bmax=%d rot_map1=%d rot_map2=%d\n", i, bmin, bmax, rot_map[i + (bmin - 1)*L], rot_map[i + (bmax + 1)*L]);
      cp = rot_map_pos[i + k * L];
      cp2 = rot_map_pos[cp.x + (cp.y + 1) * L];
      if ((i != cp2.x) || (bmin != cp2.y)) {
        ErrPrintf("ERROR: verify_rotating_map (%d %d %d %d)\n", i, cp2.x, bmin, cp2.y);
      }
      for (j = 0; j < H; j++) {
        ix1 = i + j * L + k * HL;
        ix2 = hcycle(ix1) + HL;
        ix = hcycle(ix2) - HL;
        if (ix != ix1) {
          ErrPrintf("ERROR: verify_rotating_map (ix=%d ix1=%d)\n", ix, ix1);
        }
      }
    }
  }

  LogPrintf("verify_rotating_map: done\n");
}


void translate_rotation_map(int32_t di, int32_t dk) {
  int32_t i, k, i2, k2;
  for (k = 0; k < D; k++) {
    for (i = 0; i < L; i++) {
      if (rot_map[i + k * L] > 0) {
        i2 = i + di;
        k2 = k + dk;
        if ((i2 >= 0) && (i2 < L) && (k2 >= 0) && (k2 < D)) {
          rot_map_pos[i2 + k2 * L].x = i;
          rot_map_pos[i2 + k2 * L].y = k;
          rot_map[i2 + k2 * L] = -rot_map[i + k * L];
        }
      }
    }
  }
}


void rotate_map(float angle) {
  int32_t i, k;

  //LogPrintf("rotate map - angle=%f\n", angle);

  /// creation of the maps
  if (!rot_map) {
    AllocMemoryPrint("rot_map", rot_map, char, L * D);
  }
  ResetMemory(rot_map, char, L * D);

  if (!rot_map_pos) {
    AllocMemoryPrint("rot_map_pos", rot_map_pos, Pos2, L * D);
    ResetMemory(rot_map_pos, Pos2, L * D);
    AllocMemoryPrint("rot_map_pos0", rot_map_pos0, Pos2, L * D);
    ResetMemory(rot_map_pos0, Pos2, L * D);
  } else {
    memcpy(rot_map_pos0, rot_map_pos, sizeof(Pos2)*L * D);
  }

  /// build the inner boundaries of the rotating space
  float x1 = -0.4999 * L_rot;
  float x2 = -x1;
  float y1 = -0.4999 * D_rot;
  float y2 = -y1;
  float cx = L / 2.0 - 0.5;
  float cy = D / 2.0 - 0.5;
  float co = cos(angle * PI / 180.0);
  float si = sin(angle * PI / 180.0);
  /// north-west corner
  int32_t i00 = Round(cx + co * x1 - si * y1);
  int32_t k00 = Round(cy + si * x1 + co * y1);
  /// north-east corner
  int32_t i10 = Round(cx + co * x2 - si * y1);
  int32_t k10 = Round(cy + si * x2 + co * y1);
  /// south-west corner
  int32_t i01 = Round(cx + co * x1 - si * y2);
  int32_t k01 = Round(cy + si * x1 + co * y2);
  /// inner boundary edges
  int32_t di1 = i10 - i00;
  int32_t dk1 = k10 - k00;
  int32_t di2 = i01 - i00; //-dk1;
  int32_t dk2 = k01 - k00; //di1;
  /// translate south-east corner
  int32_t i11 = i10 + di2; //roundf(di0 + co*di11 - si*dk11);
  int32_t k11 = k10 + dk2; //roundf(dk0 + si*di11 + co*dk11);

  /// set the edges (with value 1)
  set_line(rot_map, i00, k00, i01, k01, 1);
  set_line(rot_map, i00, k00, i10, k10, 1);
  set_line(rot_map, i01, k01, i11, k11, 1);
  set_line(rot_map, i10, k10, i11, k11, 1);

  /// fill the map of the rotating space (with value 2)
  int32_t bmin, bmax;
  for (k = 0; k < D; k++) {
    // locate the inner boundaries
    bmin = -1;
    bmax = -1;
    for (i = 0; i < L; i++) {
      if (rot_map[i + k * L] > 0) {
        if (bmin == -1) {
          bmin = i;
        }
        if (bmax < i) {
          bmax = i;
        }
      }
    }
    // fill the interval
    for (i = bmin + 1; i < bmax; i++) {
      if (rot_map[i + k * L] == 0) {
        rot_map[i + k * L] = 2;
      }
    }
  }

  if (abs(di1) > abs(dk1)) {
    di1 += sgn(di1);
    dk2 += sgn(dk2);
  } else {
    di2 += sgn(di2);
    dk1 += sgn(dk1);
  }
  if (pbc_mode) {
    /// set the periodic mapping on the outside of rotating space (and put negative values in rotation map)
    float rx = sqrtf(L * L + D * D) / L_rot;
    float ry = sqrtf(L * L + D * D) / D_rot;
    int32_t imax = (int)ceilf((rx - 1.0) / 2);
    int32_t kmax = (int)ceilf((ry - 1.0) / 2);
    //LogPrintf("rx=%f, ry=%f\n", rx, ry);
    //LogPrintf("imax=%d, kmax=%d\n", imax, kmax);
    for (k = -kmax; k <= kmax; k++) {
      for (i = -imax; i <= imax; i++) {
        if ((i != 0) || (k != 0)) {
          translate_rotation_map(i * di1 + k * di2, i * dk1 + k * dk2);
        }
      }
    }
  }
}


void rotation(float angle, char mode, char flags) {
  static Cell *csp_tmp = NULL;
  static Cell *TE0 = NULL;
  Cell *aux, *aux2;
  float alpha;
  float co;
  float si;
  float di0, dk0;
  float di, dk;
  int32_t i0, k0;
  int32_t i, j, k;

  if (!csp_tmp) {
    AllocMemoryPrint("csp_tmp (rotation)", csp_tmp, Cell, HLD);
  }

  if (flags & ROT_REORIENT_TEMP) {
    if (!csp_angle) {
      return;
    }
    //LogPrintf("reorientation terre\n");
    angle = -csp_angle;
  } else if (flags & ROT_REORIENT_UNDO) {
    if (!TE0) {
      return;
    }
    //LogPrintf("restauration terre\n");
    angle = csp_angle;
  } else {
    if (!angle) {
      return;
    }
    csp_angle = fmodf(csp_angle + angle, 360.0);
    LogPrintf("rotation angle: %f - orientation: %f - time: %f - mode: %s\n", angle, csp_angle, csp_time, (mode == ROT_MODE_DISK) ? "rotating table" : (mode == ROT_MODE_FULL) ? "full space" : "overlap");
  }

  aux = TE;
  aux2 = csp_tmp;

  di0 = L / 2.0 - 0.5;
  dk0 = D / 2.0 - 0.5;
  alpha = angle * PI / 180.0;
  co = cos(alpha);
  si = sin(alpha);

  /// rotating table
  if ((mode == ROT_MODE_DISK) && !(flags & ROT_REORIENT_UNDO)) {
    float dist = 10;
    float R = (D < L) ? Dd2 - dist : Ld2 - dist;
    float r2 = R * R;
    for (k = 0; k < D; k++) { // profondeur
      for (j = 0; j < H; j++) { // hauteur
        for (i = 0; i < L; i++, aux++, aux2++) { //largeur
          di = i - di0;
          dk = k - dk0;
          if (di * di + dk * dk < r2) { //modif oli
            i0 = Round(di0 + co * di + si * dk);
            k0 = Round(dk0 - si * di + co * dk);
            *aux2 = TE[i0 + j * L + k0 * HL];
          } else {
            *aux2 = *aux;
          }
        }
      }
    }
  }

  /// full rotation
  if (mode == ROT_MODE_FULL) {
    if (flags & ROT_REORIENT_TEMP) {
      rotate_map(0);
    } else {
      rotate_map(csp_angle);
    }
    if (!(flags & ROT_REORIENT_UNDO)) {
      /// fill the rotating map with cells
      Pos2 *cp;
      for (k = 0; k < D; k++) {
        for (j = 0; j < H; j++) {
          for (i = 0; i < L; i++, aux++) {
            if ((j > 0) && (j < H - 1) && (rot_map[i + k * L] > 0)) {
              di = i - di0;
              dk = k - dk0;
              i0 = Round(di0 + co * di + si * dk);
              k0 = Round(dk0 - si * di + co * dk);
              if (TE[i0 + j * L + k0 * HL].celltype != BORD) {
                csp_tmp[i + j * L + k * HL] = TE[i0 + j * L + k0 * HL];
              } else if (pbc_mode) {
                /// use the precomputed periodic mapping of the rotating space before rotation
                cp = rot_map_pos0 + (i0 + k0 * L);
                csp_tmp[i + j * L + k * HL] = TE[cp->x + j * L + cp->y * HL];
                if (csp_tmp[i + j * L + k * HL].celltype == BORD /*&& (j==2)*/) {
                  ErrPrintf("ERROR: full rotation - no cell, i=%d, k=%d, i0=%d, k0=%d, x0=%d, y0=%d\n", i, k, i0, k0, rot_map_pos[i0 + k0 * L].x, rot_map_pos0[i0 + k0 * L].y);
                  ErrPrintf("ERROR: final cell %d\n", csp_tmp[i + j * L + k * HL].celltype);
                  exit(-1);
                }
              } else {
                /// find nearest cell in rotating space, towards the center
                int32_t ix0 = i0 + j * L + k0 * HL;
                //LogPrintf("full rotation - trying to find a cell, i0=%d, k0=%d, ix0=%d\n", i0, k0, ix0);
                float dist = sqrtf(di * di + dk * dk);
                float step_i = -di / dist;
                float step_k = -dk / dist;
                //LogPrintf("di=%f, dk=%f\n",di,dk);
                //LogPrintf("step_i=%f, step_k=%f\n",step_i,step_k);
                int32_t nb;
                for (nb = 1; nb < 10; nb++) {
                  di += step_i;
                  dk += step_k;
                  i0 = Round(di0 + co * di + si * dk);
                  k0 = Round(dk0 - si * di + co * dk);
                  ix0 = i0 + j * L + k0 * HL;
                  if (TE[ix0].celltype != BORD) {
                    break;
                  }
                  //LogPrintf("i0=%d, k0=%d\n", i0, k0);
                }
                //LogPrintf("full rotation - cell found, ix0=%d, nb=%d\n", ix0, nb);
                static int32_t nbmax = 3; /// empirical value
                if (nb > nbmax) {
                  WarnPrintf("WARNING: full rotation - cell found too far, nb=%d, i=%d, k=%d, i0=%d, k0=%d\n", nb, i, k, i0, k0);
                  nbmax = nb;
                  if (nb >= 10) {
                    ErrPrintf("ERROR: full rotation - no cell found, nb=%d, i=%d, k=%d, i0=%d, k0=%d\n", nb, i, k, i0, k0);
                    exit(-1);
                  }
                }
                csp_tmp[i + j * L + k * HL] = TE[ix0];
              }
            } else {
              csp_tmp[i + j * L + k * HL] = TE[0]; //boundary cell
            }
          }
        }
      }
    }
  }

  /// periodic mapping with translations near the edges
  if ((mode == ROT_MODE_OVERLAP) && !(flags & ROT_REORIENT_UNDO)) {
    for (k = 0; k < D; k++) { // profondeur
      for (j = 0; j < H; j++) { // hauteur
        for (i = 0; i < L; i++, aux++, aux2++) { //largeur
          di = i - di0;
          dk = k - dk0;
          if ((aux->celltype != BORD) && (aux->celltype != DUM) && (aux->celltype != IN)) {
            i0 = Round(di0 + co * di + si * dk);
            k0 = Round(dk0 - si * di + co * dk);
#ifdef TRACE_DAR_COL
            char no_trace = 0;
            if ((i0 < 1) || (i0 >= L - 1) || (k0 < 1) || (k0 >= D - 1)) {
              no_trace = 1;
            }
#endif
            while (i0 < 1) {
              i0 += L - 2;
            }
            while (i0 >= L - 1) {
              i0 -= L - 2;
            }
            while (k0 < 1) {
              k0 += D - 2;
            }
            while (k0 >= D - 1) {
              k0 -= D - 2;
            }
            //aux2->celltype = (TE+i0+j*L+k0*HL)->celltype;
            *aux2 = TE[i0 + j * L + k0 * HL];
#ifdef TRACE_PAR_COL
            if (no_trace) {
              aux2->color = 0;
            }
#endif

#if defined(MODEL_DUN) || defined(MODEL_SNO)
            if (aux2->celltype == DUM) {
              aux2->celltype = EAUC;
            }
#endif
          } else
            //aux2->celltype = aux->celltype;
          {
            *aux2 = *aux;
          }
        }
      }
    }
  }

  if (flags & ROT_REORIENT_TEMP) {
    //save address of current cellular space
    TE0 = TE;
    //switch cellular space
    TE = csp_tmp;

#ifdef TRACE_PAR_COL
    //correction of the colors
    trace_par_col_rotation(0.0);
#endif
  } else if (flags & ROT_REORIENT_UNDO) {
    //restore cellular space
    TE = TE0;
    TE0 = NULL;
  } else {
    //LogPrintf("fin rotation\n");
    memcpy(TE, csp_tmp, HLD * sizeof(Cell));

    //recalcul du nombre de cellules de chaque etat
    init_Ncel();

    //unlock csp condition
    //signal_csp_ready();

    //recalcul des tableaux de doublets actifs
    init_db_pos();

#ifdef LGCA
    if (use_lgca) {
      //mise-a-jour de la matrice de gaz-sur-reseau
      collisions_mod_terre();
#ifdef CGV
      reset_grad_vel();
#endif
    }
#endif

#ifdef TRACE_PAR_COL
    //correction of the colors
    trace_par_col_rotation(csp_angle);
#endif
#ifdef TRACE_FLUX
    trace_flux_rotate(csp_angle);
#endif
  }

  if (rot_mode == ROT_MODE_FULL) {
    compute_v_walls();
  }

#ifdef ALTI
  calcule_alti(ALTI, ALTI_MODE_BAS);
  if (ava_norm) {
    calcule_normales();
  }
#endif
}

void rotation90(int32_t n) {
  static Cell *auxmap = NULL;
  Cell *aux, *aux2;
  int32_t i, j, k;
  int32_t i0, k0;
  float a = n * PI / 2;
  float co = cos(a);
  float si = sin(a);
  float di0 = L / 2.0 - 0.5;
  float dk0 = L / 2.0 - 0.5;
  float di, dk;

  if (!auxmap) {
    AllocMemoryPrint("auxmap (rotation90)", auxmap, Cell, HLD);
  }
  ResetMemory(auxmap, Cell, HLD);

  csp_angle = fmodf(csp_angle + n * 90, 360.0);
  LogPrintf("rotation90 angle : %d - orientation : %f\n", n * 90, csp_angle);

  aux = TE;
  aux2 = auxmap;
  for (k = 0; k < D; k++) { // profondeur
    for (j = 0; j < H; j++) { // hauteur
      for (i = 0; i < L; i++, aux++, aux2++) { //largeur
        if ((aux->celltype != BORD) && (aux->celltype != DUM) && (aux->celltype != IN) && (aux->celltype != OUT)) {
          di = i - di0;
          dk = k - dk0;
          i0 = Round(di0 + co * di + si * dk);
          k0 = Round(dk0 - si * di + co * dk);
          //i0 = di0 + co*(float)di - si*(float)dk;
          //k0 = dk0 + si*(float)di + co*(float)dk;
          while (i0 < 1) {
            i0 += L - 2;
          }
          while (i0 >= L - 1) {
            i0 -= L - 2;
          }
          while (k0 < 1) {
            k0 += D - 2;
          }
          while (k0 >= D - 1) {
            k0 -= D - 2;
          }
          *aux2 = TE[i0 + j * L + k0 * HL];
        } else {
          *aux2 = *aux;
        }
      }
    }
  }
  //LogPrintf("fin rotation\n");
  memcpy(TE, auxmap, HLD * sizeof(Cell));

  //recalcul des tableaux de doublets actifs
  init_db_pos();

  //recalcul du nombre de cellules de chaque etat
  init_Ncel();

#ifdef ALTI
  calcule_alti(ALTI, ALTI_MODE_BAS);
  if (ava_norm) {
    calcule_normales();
  }
#endif

#ifdef LGCA
  //mise-a-jour de la matrice de gaz-sur-reseau
  if (use_lgca) {
    collisions_mod_terre();
  }
#endif
}


void conditions_bords() {
  Cell *aux;
  int32_t j, k, k1, k2, ll;

  //ouverture aleatoire d'une partie du bord gauche (sans DUM)
  ll = D / 4;
  k1 = drand48() * (D - ll);
  k2 = k1 + ll;
  LogPrintf("bord ouvert entre %d et %d\n", k1, k2);

  for (k = 1; k < D - 1; k++) // profondeur
    for (j = 1; j < H - 1; j++) { // hauteur
      aux = TE + 1 + j * L + k * HL;
      if ((k >= k1) && (k <= k2) && (j < H - 2)) {
        aux->celltype = (aux + 1)->celltype;  //ouverture du bord gauche
      } else {
        aux->celltype = DUM;
      }
    }

  //recalcul des tableaux de doublets actifs
  init_db_pos();
}

int32_t dummy_h_plan(int32_t y) {
  int32_t i, k, ix;
  int32_t res = 1;

  //test si le plan horizontal y est entierement constitue de cellules DUM|BORD|IN|OUT
  for (k = 0; (k < D) && (res); k++) { // profondeur
    ix = y * L + k * HL;
    for (i = 0; (i < L) && (res); i++, ix++) { //largeur
      res = (TE[ix].celltype == DUM) || (TE[ix].celltype == BORD) || (TE[ix].celltype == IN) || (TE[ix].celltype == OUT);
    }
  }
  //LogPrintf("plan_h_solide %d : res = %d   ix = %d\n", y, res, ix);

  return res;
}

int32_t dummy_column(int32_t x, int32_t z) {
  int32_t ix;

  ix = z * HL + L + x;
  //test si la colonne (x,y) est entierement constitue de cellules DUM
  if (TE[ix].celltype == BORD) {
    return BORD;
  }
  while (TE[ix].celltype == DUM) {
    ix += L;
  }

  return (TE[ix].celltype == BORD) ? DUM : 0;
}

int32_t dummy_v_plan(int32_t z) {
  int32_t x = 1;
  int32_t res = 1;

  //test si le plan vertical z est entierement constitue de cellules DUM|BORD
  while ((x < L - 1) && (res = dummy_column(x, z))) {
    x++;
  }
  return res;
}


void calcule_couloir() {
  int32_t i;
  static char start = 1;

  if (!psol) {
    AllocMemory(psol, char, D);
  }
  //initialisation du tableau psol
  psol[0] = psol[D - 1] = 1;
  for (i = 1; i < D - 1; i++) {
    psol[i] = dummy_v_plan(i);
    if (!psol[i]) {
      nb_pv++;
    }
    //if (psol[i]) LogPrintf("plan solide : z = %d\n", i);
  }

  if (rot_map) {
    LN = 1;
    LS = D - 2;
  } else {
    //limite nord du couloir
    for (i = 0; (i < D) && psol[i]; i++);
    LN = i;

    //limite sud du couloir
    for (i = D - 1; (i > 0) && psol[i]; i--);
    LS = i;
  }

  //largeur du couloir
  LNS = LS - LN + 1;
  HLN = LN * HL;
  if (start) {
    LogPrintf("couloir : LN = %d   LS = %d   LNS = %d   nb_pv = %d\n", LN, LS, LNS, nb_pv);
    if (LNS <= 0) {
      //ErrPrintf("ERROR: corridor not found\n");  exit(-1);
      WarnPrintf("WARNING: dummy space !\n");
      LN = 1;
      LS = D - 2;
      LNS = LS - LN + 1;
      HLN = LN * HL;
    }
    start = 0;
  }
}

//calcul des colonnes en DUM|BORD
void compute_v_walls() {
  int32_t i, j;

  if (!csol) {
    AllocMemory(csol, char, L * D);
    ResetMemory(csol, char, L * D);
  }

  for (j = 0; j < D; j++) {
    for (i = 0; i < L; i++) {
      csol[i + j * L] = dummy_column(i, j);
    }
  }

  calcule_couloir();
}

//calcul epaisseur du plafond en DUM
void compute_h_ceil() {
  int32_t j;
  for (j = 1; (j < H) && dummy_h_plan(j); j++); // hauteur
  h_ceil = j - 1;
  LogPrintf("h_ceil = %d\n", h_ceil);
}

//calcul epaisseur du sol en DUM
void compute_h_floor() {
  int32_t j;
  for (j = H - 2; (j > 0) && dummy_h_plan(j); j--); // hauteur
  h_floor = H - 2 - j;
  LogPrintf("h_floor = %d\n", h_floor);
}

Vec3 compute_mass_center(int32_t type) {
  int32_t i, j, k, ix, n;
  Vec3 mc = {0, 0, 0};
  n = 0;
  ix = 0;
  for (k = 0; k < D; k++) {
    for (j = 0; j < H; j++) {
      for (i = 0; i < L; i++, ix++) {
        if (TE[ix].celltype == type) {
          mc.x += i;
          mc.y += j;
          mc.z += k;
          n++;
        }
      }
    }
  }

  if (n) {
    mc.x /= n;
    mc.y /= n;
    mc.z /= n;
    LogPrintf("mass center: ( %.3f , %.3f , %.3f ) - orientation: %.3f\n", mc.x, mc.y, mc.z, csp_angle);
  }

  return mc;
}

void dump_terre(char dump_type, int32_t cpt, int32_t unit) {
  char filename[512];
  char str[100];
  char *ext;

  ext = (dump_type == DUMP_CSP) ? "csp" : "bin";
  *str = 0;
  if (unit == UNIT_T0) {
    sprintf(filename, "%s%05d_t0%s.%s", MOD_NAME, cpt, str, ext);
  } else {
    sprintf(filename, "%s%04d%s.%s", MOD_NAME, cpt, str, ext);
  }

  LogPrintf("write CSP: %s, csp_time = %f (t0)\n", filename, csp_time);
  write_csp(dump_type, filename);

  compress(filename, 1);
}


#ifdef DUMP_SIGNATURE
void dump_signature(int32_t ii){
  size_t i;
  uint32_t sig, *aux;
  char output[128];

  //calcul de la signature
  sig = 0;
  aux = (unsigned int*)TE;
  for (i = 0; i < HLD * sizeof(Cell) / sizeof(unsigned int); i++, aux++) {
    sig += (*aux);  //sig = sig ^ (*aux);
  }

  //dump
  sprintf(output, "%08x : %08x\n", ii, sig);
  output_write("CELLSPACE_SIGNATURE", output);
}
#endif

