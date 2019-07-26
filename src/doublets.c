/* ReSCAL - Doublets
 *
 * Copyright (C) 2011
 *
 * Author: Olivier Rozier <rozier@ipgp.fr>
 *
 * Code based on dissol program,
 * by Eduardo Sepulveda <edo@espci.fr>
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
#include <stdint.h>

#include "defs.h"
#include "macros.h"
#include "cells.h"
#include "doublets.h"
#include "simul.h" // for output functions
#include "transitions.h"
#include "space.h"

extern Cell   *TE;	                      // notre 'terre'
#ifdef REFDB_PTR
extern RefDoublets *RefDB;                // references des cellules vers les doublets actifs
#else
extern RefDoublets_Type *RefDB_Type;      // references des cellules de la terre vers les doublets actifs
extern RefDoublets_Ind *RefDB_Ind;        // references des cellules de la terre vers les doublets actifs
#endif
extern int32_t   H, L, D, HL, HLD;        // les dimensions de la terre
extern int32_t LN, LS, LNS, HLN;          // couloir est-ouest (limite nord, limite sud, largeur nord-sud, ...)
extern const char *etats[];             // les noms des types de cellules
extern int32_t direction[];
extern int32_t orientation[];
extern int32_t pbc_mode;                  // periodic boundary conditions
extern char *rot_map;                   // periodic mapping of the rotating space

Doublet t_doub[MAX_DB];                   // la description des doublets
int32_t db_inv[3][MAX_CELL][MAX_CELL];    // la map inverse pour le type de doublet, etant
                                          // donnes le type de cell et la direction
const char *classes_db[5] = CLASSES_DB; // noms des classes de doublet
int32_t nb_type_db=0;                     // nombre de types de doublets
int32_t *db_pos[MAX_DB];                  // les tableaux contenant la position des doublets actifs
int32_t Ndb[MAX_DB];                      // nombre de doublets actifs par type
int32_t Ndbmax[MAX_DB];                   // taille des tableaux de position des doublets actifs

int32_t init_doublet(int32_t classe, char etat1, char etat2, char actif)
{
  int32_t i;

  // tester si le doublet existe deja ...
  for (i=0; i<nb_type_db; i++){
    if ((t_doub[i].classe == classe) && (t_doub[i].one == etat1) && (t_doub[i].two == etat2))
      break;
  }

  if (i==nb_type_db)
  {
    // definir nouveau doublet
    t_doub[i].classe = classe;
    t_doub[i].one = etat1;
    t_doub[i].two = etat2;
    t_doub[i].actif = actif;
    nb_type_db++;
  }
  else{
    t_doub[i].actif |= actif;
  }

  return i;
}

void split_db_hor(int32_t db_hor, int32_t *db_eo, int32_t *db_ns)
{
  uint8_t etat1, etat2, actif;

  etat1 = t_doub[db_hor].one;
  etat2 = t_doub[db_hor].two;
  actif = t_doub[db_hor].actif;

  LogPrintf("decouplage du doublet horizontal %d [ %s , %s ] ", db_hor, etats[etat1], etats[etat2]);

  *db_eo = db_inv[EST_OUEST][etat1][etat2];
  if (*db_eo < 0){
    *db_eo = init_doublet(EST_OUEST, etat1, etat2, actif);
    db_inv[EST_OUEST][etat1][etat2] = *db_eo;
  }
  else
    t_doub[*db_eo].actif |= actif;

  *db_ns = db_inv[NORD_SUD][etat1][etat2];
  if (*db_ns < 0){
    *db_ns = init_doublet(NORD_SUD, etat1, etat2, actif);
    db_inv[NORD_SUD][etat1][etat2] = *db_ns;
  }
  else
    t_doub[*db_ns].actif |= actif;


  LogPrintf("en un doublet est-ouest %d et un doublet nord-sud %d\n", *db_eo, *db_ns);

  t_doub[db_hor].actif = 0;
}


void bouclage_hor(int32_t *p_typ, int32_t *p_ind, int32_t dir)
{
  int32_t ix = *p_ind;

  (*p_ind) = hcycle(ix);
  (*p_typ) = db_inv[orientation[dir]][TE[*p_ind].celltype][TE[get_cell_dir(ix, dir)].celltype];
}

#ifdef REFDB_PTR
void translate_db_pos(int32_t type, unsigned int64_t offset)
{
  int32_t i, ix;
  LogPrintf ("translation d'adresses sur les references au tableau db_pos[%d]\n",type);
  for (i=0; i<Ndb[type]; i++){
    ix = db_pos[type][i];
    if (RefDB[ix][BAS] == db_pos[type] + i - offset)
      RefDB[ix][BAS] += offset;
    else if (RefDB[ix][EST] == db_pos[type] + i - offset)
      RefDB[ix][EST] += offset;
    else if (RefDB[ix][SUD] == db_pos[type] + i - offset)
      RefDB[ix][SUD] += offset;
    else{
      ErrPrintf ("erreur : translation d'adresse impossible sur la cellule %d\n",ix);
      exit(-1);
    }
  }

  LogPrintf ("translation d'adresses : ok !\n");
}
#endif

void realloc_db_pos(int32_t type)
{
#ifdef REFDB_PTR
  int32_t *old_adr = db_pos[type];
#endif
  Ndbmax[type] <<= 1; // on double la taille du tableau db_pos[type]
  LogPrintf ("reallocation db_pos[%d] : %d doublets\n",type,Ndbmax[type]);
  ReallocMemoryPrint("db_pos", db_pos[type], int, Ndbmax[type]);
#ifdef REFDB_PTR
  if (old_adr != db_pos[type])  //adresse du tableau db_pos[type] a peut-etre change ?
    translate_db_pos(type, ((unsigned long)db_pos[type] - (unsigned long)old_adr)/sizeof(int));
#endif
}

int32_t type_doublet(int32_t ix, int32_t dir)
{
  return db_inv[orientation[dir]][TE[ix].celltype][TE[get_cell_dir(ix, dir)].celltype];
}

void elimine_doublet(int32_t type, int32_t index, int32_t dir)
{
#ifdef CYCLAGE_HOR
  if (pbc_mode && (type == DB_BORD) && (dir != BAS)){
    bouclage_hor(&type, &index, dir);
  }
#endif

  if ((type >= 0) && t_doub[type].actif){
    int32_t *elim, *der;
#ifdef REFDB_PTR
    elim = RefDB[index][dir];
#else
    int32_t ind_elim = RefDB_Ind[index][dir];
    elim = db_pos[type] + ind_elim;
    if (*elim != index) {
      ErrPrintf("erreur adresse doublet actif dans elimine_doublet(): type=%d ind=%d dir=%d Ndb[type]=%d\n *elim=%d index=%d\n", type, ind_elim, dir, Ndb[type], *elim, index);
      exit(-1);
    }
#endif
    if (elim){
#ifdef REFDB_PTR
      RefDB[index][dir] = NULL;
#else
      RefDB_Type[index][dir] = -1;
      RefDB_Ind[index][dir] = -1;
#endif
      Ndb[type]--;
      if (Ndb[type] < 0){
        extern uint64_t iter;
        ErrPrintf("ERROR: Ndb[%d] = %d\n", type, Ndb[type]);
        ErrPrintf("index = %d\n", index);
        int32_t x,y,z;
        Calcule_xyz(index,x,y,z);
        ErrPrintf("x=%d   y=%d   z=%d\n",x,y,z);
        ErrPrintf("iter=%" PRIu64 "\n", iter);
        exit(-1);
      }
      // on bouche le 'trou' avec le dernier doublet du tableau
      der = db_pos[type] + Ndb[type];
      if (elim != der){
        *elim = *der;
#ifdef REFDB_PTR
        if (dir == BAS)
          RefDB[*der][BAS] = elim; //vertical
        else if (RefDB[*der][EST] == der)
          RefDB[*der][EST] = elim; //horizontal est-ouest
        else
          RefDB[*der][SUD] = elim; //horizontal nord-sud
#else
        if (dir == BAS)
          RefDB_Ind[*der][BAS] = ind_elim; //vertical
        else if ((RefDB_Type[*der][EST]== type) && (RefDB_Ind[*der][EST] == Ndb[type]))
          RefDB_Ind[*der][EST] = ind_elim; //horizontal est-ouest
        else
          RefDB_Ind[*der][SUD] = ind_elim; //horizontal nord-sud
#endif
      }
    }
#ifndef PARALLEL
    else{
      ErrPrintf("WARNING: active doublet not referenced (index = %d, dir = %d)\n", index, dir);
    }
#endif
  }
}

void elimine_doublet_est(int32_t index)
{
#ifdef REFDB_PTR
  int32_t type = db_inv[EST_OUEST][TE[index].celltype][TE[index+1].celltype];
#else
  int32_t type = RefDB_Type[index][EST];
#endif
  elimine_doublet(type, index, EST);
}

void elimine_doublet_ouest(int32_t index)
{
#ifdef REFDB_PTR
  int32_t type = db_inv[EST_OUEST][TE[index-1].celltype][TE[index].celltype];
#else
  int32_t type = RefDB_Type[index-1][EST];
#endif
  elimine_doublet(type, index-1, EST);
}

void elimine_doublet_bas(int32_t index)
{
  if (H <= 3) return;
#ifdef REFDB_PTR
  int32_t type = db_inv[VERTICAL][TE[index].celltype][TE[index+L].celltype];
#else
  int32_t type = RefDB_Type[index][BAS];
#endif
  elimine_doublet(type, index, BAS);
}

void elimine_doublet_haut(int32_t index)
{
  if (H <= 3) return;
#ifdef REFDB_PTR
  int32_t type = db_inv[VERTICAL][TE[index-L].celltype][TE[index].celltype];
#else
  int32_t type = RefDB_Type[index-L][BAS];
#endif
  elimine_doublet(type, index-L, BAS);
}

void elimine_doublet_sud(int32_t index)
{
  if (D <= 3) return;
#ifdef REFDB_PTR
  int32_t type = db_inv[NORD_SUD][TE[index].celltype][TE[index+HL].celltype];
#else
  int32_t type = RefDB_Type[index][SUD];
#endif
  elimine_doublet(type, index, SUD);
}

void elimine_doublet_nord(int32_t index)
{
  if (D <= 3) return;
#ifdef REFDB_PTR
  int32_t type = db_inv[NORD_SUD][TE[index-HL].celltype][TE[index].celltype];
#else
  int32_t type = RefDB_Type[index-HL][SUD];
#endif
  elimine_doublet(type, index-HL, SUD);
}

void ajoute_doublet(int32_t type, int32_t index, int32_t dir)
{
#ifdef CYCLAGE_HOR
  if (pbc_mode && (type == DB_BORD) && (dir != BAS))
    bouclage_hor(&type, &index, dir);
#endif

  if ((type >= 0) && t_doub[type].actif){
    if (Ndb[type] >= Ndbmax[type]) // tableau db_pos pas assez grand ?
      realloc_db_pos(type);
    // on ajoute le doublet en fin de tableau
#ifdef REFDB_PTR
    if (RefDB[index][dir]){
      ErrPrintf("ERROR: doublet reference already exists (index = %d, dir = %d, new type = %d)\n", index, dir, type);
      exit(-1);
    }
    RefDB[index][dir] = db_pos[type] + Ndb[type];
    db_pos[type][Ndb[type]++] = index;
#else
    if (RefDB_Type[index][dir] != -1){
      ErrPrintf("ERROR: doublet reference already exists (index = %d, dir = %d, old type = %d, new type = %d)\n", index, dir, RefDB_Type[index][dir], type);
      exit(-1);
    }
    RefDB_Type[index][dir] = type;
    RefDB_Ind[index][dir] = Ndb[type];
    db_pos[type][Ndb[type]++] = index;
#endif
  }
}


void ajoute_doublet_est(int32_t index)
{
  int32_t type = db_inv[EST_OUEST][TE[index].celltype][TE[index+1].celltype];
  ajoute_doublet(type, index, EST);
}

void ajoute_doublet_ouest(int32_t index)
{
  int32_t type = db_inv[EST_OUEST][TE[index-1].celltype][TE[index].celltype];
  ajoute_doublet(type, index-1, EST);
}

void ajoute_doublet_bas(int32_t index)
{
  if (H <= 3) return;
  int32_t type = db_inv[VERTICAL][TE[index].celltype][TE[index+L].celltype];
  ajoute_doublet(type, index, BAS);
}

void ajoute_doublet_haut(int32_t index)
{
  if (H <= 3) return;
  int32_t type = db_inv[VERTICAL][TE[index-L].celltype][TE[index].celltype];
  ajoute_doublet(type, index-L, BAS);
}

void ajoute_doublet_sud(int32_t index)
{
  if (D <= 3) return;
  int32_t type = db_inv[NORD_SUD][TE[index].celltype][TE[index+HL].celltype];
  ajoute_doublet(type, index, SUD);
}

void ajoute_doublet_nord(int32_t index)
{
  if (D <= 3) return;
  int32_t type = db_inv[NORD_SUD][TE[index-HL].celltype][TE[index].celltype];
  ajoute_doublet(type, index-HL, SUD);
}

void init_db_inv()
{
  int32_t i, n;
  int32_t cl, etat1, etat2, db_eo, db_ns;
  int32_t split;
  char flag_conflict;

  // initialisation du tableau db_inv[][][]
  memset(db_inv, DB_INACTIF, sizeof(db_inv));

  // doublets avec une direction bien definie
  for (i=0; i<nb_type_db; i++){
    cl = t_doub[i].classe;
    if ((cl != HORIZONTAL) /*&& t_doub[i].actif*/){
      etat1 = t_doub[i].one;
      etat2 = t_doub[i].two;
      if (db_inv[cl][etat1][etat2] >= 0){
        ErrPrintf("ERROR : duplicated doublet [ %s , %s , %s ] -> %d, %d\n", classes_db[cl], etats[etat1], etats[etat2], db_inv[cl][etat1][etat2], i);

      }
      db_inv[cl][etat1][etat2] = i;
    }
  }

  // cas particulier des doublets dits 'horizontaux'
  n = 0;
  do {
    split = 0;
    LogPrintf("decouplage des doublets horizontaux : passe %d\n", n);
    for (i=0; i<nb_type_db; i++){
      cl = t_doub[i].classe;
      if (cl == HORIZONTAL && t_doub[i].actif){
        etat1 = t_doub[i].one;
        etat2 = t_doub[i].two;
        db_eo = db_inv[EST_OUEST][etat1][etat2];
        db_ns = db_inv[NORD_SUD][etat1][etat2];
        flag_conflict = ((db_eo >= 0) && t_doub[db_eo].actif) || ((db_ns >= 0) && t_doub[db_ns].actif);
        if (flag_conflict){
          split |= split_trans_hor(i);
        }
      }
    }
    n++;
  }
  while (split); //autant de passes que necessaire

  for (i=0; i<nb_type_db; i++){
    cl = t_doub[i].classe;
    if (cl == HORIZONTAL){
      etat1 = t_doub[i].one;
      etat2 = t_doub[i].two;
      db_eo = db_inv[EST_OUEST][etat1][etat2];
      db_ns = db_inv[NORD_SUD][etat1][etat2];
      flag_conflict = ((db_eo >= 0) && t_doub[db_eo].actif) || ((db_ns >= 0) && t_doub[db_ns].actif);
      if (!flag_conflict){
        db_inv[EST_OUEST][etat1][etat2] = i;
        db_inv[NORD_SUD][etat1][etat2] = i;
      }
    }
  }

#ifdef CYCLAGE_HOR
  if (pbc_mode){
    for (i=0; i<MAX_CELL; i++){
      db_inv[EST_OUEST][i][BORD] = DB_BORD;
      db_inv[EST_OUEST][BORD][i] = DB_BORD;
      db_inv[NORD_SUD][i][BORD] = DB_BORD;
      db_inv[NORD_SUD][BORD][i] = DB_BORD;
      db_inv[VERTICAL][i][BORD] = DB_BORD;
      db_inv[VERTICAL][BORD][i] = DB_BORD;
    }
  }
#endif
}

void init_db_pos()
{
  static char first = 1;
  int32_t i,j,k, ix, td, tot;
  Cell *t,*dr, *ba, *de;
  // allocations pour les tableaux de positions db_pos[][]
  tot = 0;
  for(i = 0; i < nb_type_db; i++){
    Ndb[i] = 0;
    if (first){
      Ndbmax[i] = 2*L*D*t_doub[i].actif;
      if (t_doub[i].actif){
        AllocMemory(db_pos[i], int, Ndbmax[i]);
      }
    }
    tot += Ndbmax[i];
  }
  if (first){
    LogPrintf("allocation db_pos : %lu\n", tot*sizeof(int));
    PrintTotalMemory();
    first = 0;
  }

  // initialisation des doublets actifs
#ifdef REFDB_PTR
  memset(RefDB, 0, sizeof(RefDoublets) * HLD);
#else
  memset(RefDB_Type, DB_INACTIF, sizeof(RefDoublets_Type) * HLD);
  memset(RefDB_Ind, -1, sizeof(RefDoublets_Ind) * HLD);
#endif
  t = TE;
  dr = t+1; ba = t+L; de = t+HL;
  ix = 0;
  for(k = 0; k < D; k++){
    for(j = 0; j < H; j++){
      for(i = 0; i < L; i++, t++, dr++, ba++, de++, ix++){
#ifdef CYCLAGE_HOR
        if (pbc_mode){
          if (!rot_map){
            if ((k == 0) || (i == 0)){
#ifndef REFDB_PTR
              RefDB_Type[ix][EST] = DB_BORD;
              RefDB_Type[ix][SUD] = DB_BORD;
#endif
              continue;
            }
          }
          else if (OutOfSpace(i,k)){
#ifndef REFDB_PTR
            RefDB_Type[ix][EST] = DB_BORD;
            RefDB_Type[ix][SUD] = DB_BORD;
#endif
            continue;
          }
        }
#endif
        if ((k < D-1) && (j < H-1) && (i < L-1)) // on ne sort pas de la terre !!
        {
          if (j > 0){
            td = db_inv[EST_OUEST][t->celltype][dr->celltype];
            ajoute_doublet(td, ix, EST);

            if (LNS > 1){
              td = db_inv[NORD_SUD][t->celltype][de->celltype];
              ajoute_doublet(td, ix, SUD);
            }
          }

          td = db_inv[VERTICAL][t->celltype][ba->celltype];
          ajoute_doublet(td, ix, BAS);
        }
      }
    }
  }
}

void fin_db_pos()
{
  int32_t i;
  LogPrintf("liberation db_pos\n");
  for(i = 0; i < nb_type_db; i++){
    if (t_doub[i].actif) FreeMemory(db_pos[i], int, Ndbmax[i]);
  }
}


void dump_doublets(){
//Write information about doublets to DOUBLETS.log
  char current_output[256];

  output_write("DOUBLETS", "\n# NUMBER OF ACTIVE DOUBLETS\n");
  output_write("DOUBLETS", "\n     ");
  output_write("DOUBLETS", "\n# DOUBLETS\n");
  sprintf(current_output, "\nNB_TYPE_DOUBLET = %d\n",nb_type_db);
  output_write("DOUBLETS", current_output);
  for(int i=0; i<nb_type_db; i++){
    sprintf(current_output,"DB(%2d): %s, [%s, %s] %s\n", i, classes_db[t_doub[i].classe], etats[t_doub[i].one], etats[t_doub[i].two], t_doub[i].actif?"active":"");
    output_write("DOUBLETS", current_output);
  }
}

#ifdef INFO_DBL
void dump_db_info(){
//Write more info about doublets to DOUBLETS.log
  static int32_t cpt = 0; // counts calls to this function
  int32_t tot,totmax;
  char current_output[256];

  // Things to do the first time this function is called
  if (!cpt){
    for (int32_t i=0; i<nb_type_db; i++){
      if (t_doub[i].actif){
        sprintf(current_output,"     DB(%2d)", i);
        output_write("DOUBLETS", current_output);
      }
    }
  }

  // Things to do every time
  sprintf(current_output,"\n%04d:",cpt);
  output_write("DOUBLETS", current_output);
  for (int32_t i=tot=totmax=0; i<nb_type_db; i++)
  {
    if (t_doub[i].actif){
      sprintf(current_output," %10d",Ndb[i]);
      output_write("DOUBLETS", current_output);
      tot += Ndb[i];
      totmax += Ndbmax[i];
    }
  }
  cpt++;
}
#endif
