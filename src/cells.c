/* ReSCAL - Cells
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
#include <string.h>
#include <stdint.h>

#include "defs.h"
#include "macros.h"
#include "cells.h"
#include "space.h"
#include "surface.h"
#include "trace.h"

#ifdef LGCA
#include "lgca.h"
#endif

extern Cell   *TE;	      // notre 'terre'
extern int32_t HLD;       // les dimensions de la terre
//extern int32_t H, L, HL, HLD;       // les dimensions de la terre
//extern RefCellules *RefCel;      // references des cellules vers les cellules actives
extern const uint8_t Phase[MAX_CELL];	//phase (fluide ou solide) des types de cellules
#ifdef CELL_TIME
extern double csp_time;
#endif
#ifdef ALTI
extern int32_t LN;
extern int16_t *alti;
#endif
extern int32_t use_lgca;

const int8_t *etats[MAX_CELL] = ETATS;  // les noms des types de cellules

//int8_t t_cell[MAX_CELL];        // cellules actives
//int32_t *cel_pos[MAX_CELL];                   // les tableaux contenant la position des cellules actives
//int32_t Ncelmax[MAX_CELL];                    // taille des tableaux de position des cellules actives
int32_t Ncel[MAX_CELL];                       // nombre de cellules par type

#ifdef PARALLEL
int32_t Ncel_par[MAX_CELL];           // nombre total de cellules par type sur l'ensemble des process
extern int32_t mode_par;              // mode parallele
extern uint8_t pserv;       // flag process serveur
#endif


void init_Ncel()
{
  int32_t i;

  for (i=0; i<MAX_CELL; i++)
    Ncel[i]=0;

  for(i=0; i<HLD; i++){
    Ncel[TE[i].celltype]++;
    //ajoute_cellule(TE[i].celltype);
  }
  //LogPrintf ("Ncel[MOINS] = %d\n",Ncel[MOINS]);
}

void init_cellule(Cell cel, int32_t index)
{
  TE[index] = cel;
#ifdef CELL_TIME
  if (cel.celltype == CELL_TIME)
    TE[index].celltime = (int)csp_time;
  else
    TE[index].celltime = 0;
#endif
#ifdef LGCA
  if (use_lgca) collisions_modcell(cel.celltype, index);
#endif
}

void modifie_cellule(int32_t type, int32_t index)
{
  int32_t old_type;

  old_type = TE[index].celltype;

  Ncel[old_type]--;

  TE[index].celltype = type;

  Ncel[type]++;

#ifdef CELL_TIME
  if (type == CELL_TIME)
    TE[index].celltime = (int)csp_time;
  else
    TE[index].celltime = 0;
#endif

#ifdef ALTI
  //if ((ALTI == old_type) || (ALTI == type)){  //recalcul de l'altitude
  if (Phase[old_type] != Phase[type]){  //recalcul de l'altitude
    modif_alti_cel(index, type);
  }
#endif
#ifdef LGCA
  if (use_lgca) collisions_modcell(type, index);
#endif
}

/*
void echange_cellules(int32_t ix, int32_t ix2)
{
  Cell aux;

  aux = TE[ix];
  TE[ix] = TE[ix2];
  TE[ix2] = aux;
}
*/


//deplace la cellule ix en ix2
void deplace_cellule(int32_t ix, int32_t ix2)
{
  TE[ix2] = TE[ix];
#ifdef CELL_TIME
  if (TE[ix].celltype == CELL_TIME)
    TE[ix2].celltime = (int)csp_time;
  else
    TE[ix2].celltime = 0;
#endif
#ifdef LGCA
  if (use_lgca) collisions_modcell(TE[ix].celltype, ix2);
#endif
}


#ifdef MODEL_ICB
void verifier_Ncel_MOINS()
{
  int32_t i,N;

  for(i=N=0; i<HLD; i++)
    if (TE[i].celltype == MOINS) N++;

  if (N!=Ncel[MOINS]){
    ErrPrintf ("ERROR: %d cells MOINS, Ncel[MOINS] = %d \n", N, Ncel[MOINS]);
    //exit(-1);
  }
  else{
    LogPrintf ("verification Ncel[MOINS] ok\n");
  }
}
#endif


void log_cell()
{
  FILE *fp;
  int32_t i, nb;
  int8_t *filename = "CELL.log";
  static int32_t start = 1;
  static int32_t cpt = 0;

  if (start){
    fp = fopen(filename,"w");
    if ( ! fp ){
      ErrPrintf("ERROR: cannot open file CELL.log\n");
      exit(-4);
    }

    fprintf(fp,"\n# CELL STATES\n");

    nb = 0;
    for(i=0; i<MAX_CELL; i++){
      if (etats[i]) nb++;
    }
    fprintf(fp,"\nNB_STATES = %d\n", nb);

    for(i=0; i<MAX_CELL; i++){
      if (etats[i]){
#ifdef PHASES
       int8_t *str_phase = (Phase[i]==SOLID)?"<SOLID>":"<FLUID>";
#else
       int8_t *str_phase = NULL;
#endif
        fprintf(fp,"ST(%d): %s %s\n", i, etats[i], str_phase);
      }
    }

    fprintf(fp,"\n     ");
    int8_t str[30];
    for(i=0; i<MAX_CELL; i++){
      if (etats[i]){
        str[0] = 0;
        strncat(str, etats[i], 10);
        fprintf(fp, "%11s", str);
      }
    }

    fprintf(fp,"\n");
    fclose(fp);
    start = 0;
  }//start
#ifdef INFO_CEL
  else{
    //wait_csp_ready();

    fp = fopen(filename,"a");

    fprintf(fp,"%04d:",cpt++);
    for (i=0; i<MAX_CELL; i++){
      if (etats[i]){
        fprintf(fp,"  %09d", Ncel[i]);
      }
    }
    fprintf(fp,"\n");

    fclose(fp);

#ifdef PARALLEL
  //TODO : comprendre pourquoi cela ne marche pas bien en mode MPI !?
  /*  if (mode_par){
      synchro_Ncel(); //consolidation des donnees globales
      if (pserv){
        sprintf(nom,"INFO_CEL_PAR.log");
        fp = fopen(nom,"a");
        fprintf(fp,"%04d :",cpt-1);
        for (i=0; i<MAX_CELL; i++){
          if (Ncel_par[i]) fprintf(fp,"\tn[%d]=%d",i,Ncel_par[i]);
        }
        fprintf(fp,"\n");
        fclose(fp);
      }
  }*/
#endif  //PARALLEL
  }
#endif //INFO_CEL
}

