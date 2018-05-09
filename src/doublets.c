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
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */


#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <assert.h>
#include "defs.h"
#include "macros.h"
#include "cells.h"
#include "doublets.h"
#include "transitions.h"
#include "space.h"
#ifdef PARALLEL
#include "synchro.h"
#endif

extern Cell   *TE;	        // notre 'terre'
#ifdef REFDB_PTR
extern RefDoublets *RefDB;      // references des cellules vers les doublets actifs
#else
extern RefDoublets_Type *RefDB_Type;       // references des cellules de la terre vers les doublets actifs
extern RefDoublets_Ind *RefDB_Ind;       // references des cellules de la terre vers les doublets actifs
#endif
extern int   H, L, D, HL, HLD;     // les dimensions de la terre
extern int LN, LS, LNS, HLN;    //couloir est-ouest (limite nord, limite sud, largeur nord-sud, ...)
extern const char *etats[];     // les noms des types de cellules
extern int direction[];
extern int orientation[];
extern int pbc_mode;  //periodic boundary conditions
extern char *rot_map;        // periodic mapping of the rotating space
#ifdef PARALLEL
extern int mode_par;    //mode parallele
#endif

Doublet t_doub[MAX_DB];               // la description des doublets
int db_inv[3][MAX_CELL][MAX_CELL];    // la map inverse pour le type de doublet, etant
                                      // donnes le type de cell et la direction
const char *classes_db[5] = CLASSES_DB; //noms des classes de doublet
int nb_type_db=0;                 // nombre de types de doublets
int *db_pos[MAX_DB];            // les tableaux contenant la position des doublets actifs
int Ndb[MAX_DB];                // nombre de doublets actifs par type
int Ndbmax[MAX_DB];             // taille des tableaux de position des doublets actifs

int init_doublet(int classe, char etat1, char etat2, char actif)
{
  int i;

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

void split_db_hor(int db_hor, int *db_eo, int *db_ns)
{
  unsigned char etat1, etat2, actif;

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


void bouclage_hor(int *p_typ, int *p_ind, int dir)
{
  int ix = *p_ind;

  (*p_ind) = hcycle(ix);
  (*p_typ) = db_inv[orientation[dir]][TE[*p_ind].celltype][TE[get_cell_dir(ix, dir)].celltype];
}

#ifdef REFDB_PTR

void translate_db_pos(int type, unsigned long offset)
{
  int i, ix;

  LogPrintf ("translation d'adresses sur les references au tableau db_pos[%d]\n",type);
#if 0
  extern int arch;

  if (arch == W64){
    LogPrintf ("offset = 0x%016lx\n", offset);
    LogPrintf ("nouvelle adresse = 0x%016x\n", (unsigned long)db_pos[type]);
  }
  else{
    LogPrintf ("offset = 0x%08lx\n", offset);
    LogPrintf ("nouvelle adresse = 0x%08x\n", (unsigned int)db_pos[type]);
  }
#endif

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

void realloc_db_pos(int type)
{
#ifdef REFDB_PTR
  int *old_adr = db_pos[type];
#endif

  //Ndbmax[type] += (Ndbmax[type]>>1); // on augmente de 50% la taille du tableau db_pos[type]
  Ndbmax[type] <<= 1; // on double la taille du tableau db_pos[type]

  LogPrintf ("reallocation db_pos[%d] : %d doublets\n",type,Ndbmax[type]);

  ReallocMemoryPrint("db_pos", db_pos[type], int, Ndbmax[type]);

#ifdef REFDB_PTR
  if (old_adr != db_pos[type])  //adresse du tableau db_pos[type] a peut-etre change ?
    translate_db_pos(type, ((unsigned long)db_pos[type] - (unsigned long)old_adr)/sizeof(int));
#endif
}

int type_doublet(int ix, int dir)
{
  return db_inv[orientation[dir]][TE[ix].celltype][TE[get_cell_dir(ix, dir)].celltype];
}

void elimine_doublet(int type, int index, int dir)
{
#ifdef CYCLAGE_HOR
  if (pbc_mode && (type == DB_BORD) && (dir != BAS)){
    bouclage_hor(&type, &index, dir);
  }
#endif

  if ((type >= 0) && t_doub[type].actif){
    int *elim, *der;
#ifdef REFDB_PTR
    elim = RefDB[index][dir];
#else
    int ind_elim = RefDB_Ind[index][dir];
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
      //if (type==11){LogPrintf("Ndb[%d] = %d\n", type, Ndb[type]);}
      if (Ndb[type] < 0){
        extern unsigned long int iter;
        ErrPrintf("ERROR: Ndb[%d] = %d\n", type, Ndb[type]);
        ErrPrintf("index = %d\n", index);
        int x,y,z;
        Calcule_xyz(index,x,y,z);
        ErrPrintf("x=%d   y=%d   z=%d\n",x,y,z);
        ErrPrintf("iter=%ld\n", iter);
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

void elimine_doublet_est(int index)
{
#ifdef REFDB_PTR
  int type = db_inv[EST_OUEST][TE[index].celltype][TE[index+1].celltype];
#else
  int type = RefDB_Type[index][EST];
#endif
  elimine_doublet(type, index, EST);
}

void elimine_doublet_ouest(int index)
{
#ifdef REFDB_PTR
  int type = db_inv[EST_OUEST][TE[index-1].celltype][TE[index].celltype];
#else
  int type = RefDB_Type[index-1][EST];
#endif
  elimine_doublet(type, index-1, EST);
}

void elimine_doublet_bas(int index)
{
  if (H <= 3) return;
#ifdef REFDB_PTR
  int type = db_inv[VERTICAL][TE[index].celltype][TE[index+L].celltype];
#else
  int type = RefDB_Type[index][BAS];
#endif
  elimine_doublet(type, index, BAS);
}

void elimine_doublet_haut(int index)
{
  if (H <= 3) return;
#ifdef REFDB_PTR
  int type = db_inv[VERTICAL][TE[index-L].celltype][TE[index].celltype];
#else
  int type = RefDB_Type[index-L][BAS];
#endif
  elimine_doublet(type, index-L, BAS);
}

void elimine_doublet_sud(int index)
{
  if (D <= 3) return;
#ifdef REFDB_PTR
  int type = db_inv[NORD_SUD][TE[index].celltype][TE[index+HL].celltype];
#else
  int type = RefDB_Type[index][SUD];
#endif
  elimine_doublet(type, index, SUD);
}

void elimine_doublet_nord(int index)
{
  if (D <= 3) return;
#ifdef REFDB_PTR
  int type = db_inv[NORD_SUD][TE[index-HL].celltype][TE[index].celltype];
#else
  int type = RefDB_Type[index-HL][SUD];
#endif
  elimine_doublet(type, index-HL, SUD);
}

void ajoute_doublet(int type, int index, int dir)
{
#ifdef CYCLAGE_HOR
  if (pbc_mode && (type == DB_BORD) && (dir != BAS))
    bouclage_hor(&type, &index, dir);
#endif

#ifdef PARALLEL
  if (mode_par && (type == DB_TUNNEL)){
    int tc = TE[index].celltype;  //t_doub[type].one;
    //int tc2 = t_doub[type].two;
    //LogPrintf("TUNNEL : tc=%d  \n", tc);
    if (dir == EST){
      if (tc != TUNNEL)
        synchro_tunnel_est(tc, index);
      else
        synchro_tunnel_ouest(TE[index+1].celltype, index+1);
    }
    else if (dir == SUD){
      if (tc != TUNNEL)
        synchro_tunnel_sud(tc, index);
      else
        synchro_tunnel_nord(TE[index+HL].celltype, index+HL);
    }
  }
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


void ajoute_doublet_est(int index)
{
  int type = db_inv[EST_OUEST][TE[index].celltype][TE[index+1].celltype];
  ajoute_doublet(type, index, EST);
}

void ajoute_doublet_ouest(int index)
{
  int type = db_inv[EST_OUEST][TE[index-1].celltype][TE[index].celltype];
  ajoute_doublet(type, index-1, EST);
}

void ajoute_doublet_bas(int index)
{
  if (H <= 3) return;
  int type = db_inv[VERTICAL][TE[index].celltype][TE[index+L].celltype];
  ajoute_doublet(type, index, BAS);
}

void ajoute_doublet_haut(int index)
{
  if (H <= 3) return;
  int type = db_inv[VERTICAL][TE[index-L].celltype][TE[index].celltype];
  ajoute_doublet(type, index-L, BAS);
}

void ajoute_doublet_sud(int index)
{
  if (D <= 3) return;
  int type = db_inv[NORD_SUD][TE[index].celltype][TE[index+HL].celltype];
  ajoute_doublet(type, index, SUD);
}

void ajoute_doublet_nord(int index)
{
  if (D <= 3) return;
  int type = db_inv[NORD_SUD][TE[index-HL].celltype][TE[index].celltype];
  ajoute_doublet(type, index-HL, SUD);
}

void init_db_inv()
{
  int i, n;
  int cl, etat1, etat2, db_eo, db_ns;
  int split;
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
        //if (((db_eo >= 0) ) || ((db_ns >= 0) )){
        if (flag_conflict){
          //split_db_hor(i);
          //LogPrintf("init_db_inv: etat1=%d, etat2=%d, db_eo=%d db_ns=%d\n", etat1, etat2, db_eo, db_ns);
          split |= split_trans_hor(i);
          //split = 0;
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

#ifdef PARALLEL
  // cas particulier des tunnels
  for (i=0; i<MAX_CELL; i++){
    db_inv[EST_OUEST][i][TUNNEL] = DB_TUNNEL;
    db_inv[EST_OUEST][TUNNEL][i] = DB_TUNNEL;
    db_inv[NORD_SUD][i][TUNNEL] = DB_TUNNEL;
    db_inv[NORD_SUD][TUNNEL][i] = DB_TUNNEL;
    db_inv[VERTICAL][i][TUNNEL] = DB_TUNNEL;
    db_inv[VERTICAL][TUNNEL][i] = DB_TUNNEL;
  }
#endif
}

void init_db_pos()
{
  static char first = 1;
  unsigned int i,j,k, ix, td, tot;
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
    //pos[i] = (int*) malloc(nmax[i]*sizeof(int) );
    //LogPrintf ("Allocation de pos[%d] : %ld\n", i, nmax[i]*sizeof(int));
    //LogPrintf("adresse de pos[%d] : %08x\n",i,pos[i]);
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
  //t = TE;
  t = TE; // + LN*HL;
  dr = t+1; ba = t+L; de = t+HL;
  ix = 0;
  //ix = LN*HL;

  //for(k = 0; k < L; k++){
  for(k = 0; k < D; k++){
    for(j = 0; j < H; j++){
      for(i = 0; i < L; i++, t++, dr++, ba++, de++, ix++){
#ifdef CYCLAGE_HOR
        //if (!pbc_mode || ((k > 0) && (i > 0))) //pour eviter la surcharge aux bords !!
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
#ifdef PARALLEL
        if ((dr->celltype != TUNNEL) && (ba->celltype != TUNNEL) && (de->celltype != TUNNEL)) //le doublet n'est pas dans une zone de recouvrement asservie
#endif
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
  //LogPrintf("init_db_pos : fin\n");
}

void fin_db_pos()
{
  int i;
  LogPrintf("liberation db_pos\n");
  for(i = 0; i < nb_type_db; i++){
    //if (db_pos[i]) free(db_pos[i]);
    if (t_doub[i].actif) FreeMemory(db_pos[i], int, Ndbmax[i]);
  }
}


void dump_doublets()
{
  FILE *fp;
  int i;

  fp = fopen("DBL.log","w");
  if ( ! fp ){
	  ErrPrintf("erreur ouverture fichier dump DBL.log\n");
	  exit(-4);
  }

  fprintf(fp,"\n# DOUBLETS\n");
  fprintf(fp,"\nNB_TYPE_DOUBLET = %d\n",nb_type_db);
  for(i=0; i<nb_type_db; i++){
    fprintf(fp,"DB(%2d): %s, [%s, %s] %s\n", i, classes_db[t_doub[i].classe], etats[t_doub[i].one], etats[t_doub[i].two], t_doub[i].actif?"active":"");
    //fprintf(fp,"DB(%2d): %s, [%s, %s] %s | %d\n", i, classes_db[t_doub[i].classe], etats[t_doub[i].one], etats[t_doub[i].two], t_doub[i].actif?"active":"", (t_doub[i].classe <= 2)?db_inv[t_doub[i].classe][t_doub[i].one][t_doub[i].two]:-1);  //DEBUG
  }

  fclose(fp);
}

#ifdef INFO_DBL
void dump_db_info()
{
  static int cpt = 0;
  FILE *fp;
  int i,tot,totmax;

  fp = fopen("DBL.log","a");

  if (!cpt){
    fprintf(fp,"\n# NUMBER OF ACTIVE DOUBLETS\n");
    fprintf(fp,"\n     ");
    for (i=0; i<nb_type_db; i++){
      if (t_doub[i].actif) fprintf(fp,"     DB(%2d)", i);
    }
  }

  fprintf(fp,"\n%04d:",cpt);
  for (i=tot=totmax=0; i<nb_type_db; i++)
  {
    if (t_doub[i].actif){
      fprintf(fp," %10d",Ndb[i]);
      tot += Ndb[i];
      totmax += Ndbmax[i];
    }
  }
  //fprintf(fp,"     (total = %d, max = %d, alloc = %lu)",tot,totmax,totmax*sizeof(int)); //debug info

  fclose(fp);

  cpt++;
}
#endif
