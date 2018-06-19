/* ReSCAL - genesis entry
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
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <stdint.h>

#include "defs.h"
#include "macros.h"
#include "genesis.h"
#include "cells.h"
#include "param.h"
#include "format.h"


int32_t prog=PROG_TOOL;
int32_t  H, L, D, HL, HLD;        // size of cellular space
Cell *TE=NULL;          // cellular space
double csp_time=0.0;
//int32_t  node = -1;            // node number (with PARALLEL option)
int8_t *model_name = NULL; //name of model
int8_t *bin_filename = NULL; //name of BIN file
int8_t *csp_filename = NULL; //name of CSP file
int8_t opt_bin = 0;   //raw binary format
int32_t graine = 68374; //random seed

int8_t *csp_template_str = NULL; //value of CSP_template parameter
CSP_Template csp_template; //CSP template
CSP_Template t_templates[MAX_TEMPLATES]; //available templates
int32_t nb_templates = 0; //number of templates

int8_t *boundary_str = NULL;
int32_t boundary = BC_PERIODIC;  //boundary conditions

double drand48();
void  srand48();

float* load_surface(int8_t *, int, int);

void init_template(int32_t type, int8_t *name, int8_t *desc, int32_t nb_args, ...)
{
  va_list vl;
  int32_t i;

  t_templates[nb_templates].type = type;
  t_templates[nb_templates].name = name;
  t_templates[nb_templates].desc = desc;
  t_templates[nb_templates].nb_args = nb_args;

  AllocMemory(t_templates[nb_templates].args, float, nb_args);
  va_start(vl, nb_args);
  for (i=0; i<nb_args; i++){
    t_templates[nb_templates].args[i] = (float) va_arg(vl, double);
    //LogPrintf("t_templates[%d].args[%d] = %f\n", nb_templates, i, t_templates[nb_templates].args[i]);
  }
  va_end(vl);
  nb_templates++;
}

// TEMPLATES USER'S GUIDE
// ----------------------
//
// A CSP template is a predefined cellular space configuration, with optional arguments.
// The user gives to the parameter Csp_template a value of the form:
// TEMPLATE_NAME, or TEMPLATE_NAME(arg1, ...).
//
// Example: LAYER(5) will generate a sand layer of 5 cells height
//
// Here we explain how to create a new CSP template:
// 1) add a new type in CSP_TEMPLATES.
// 2) define the new template in init_template_list() by adding a line of the form
//    init_template(template type, "TEMPLATE NAME", "Template usage description", number of optional arguments, default value for argument 1, default value for argument 2, ...);
//    Note: default values of arguments MUST be floats (in case of an integer, just add a zero fractional part).
// 3) in genesis() function, write the implementation for the new CSP template.
//

//available CSP template types
#if defined(MODEL_DUN) || defined(MODEL_SNO)
enum CSP_TEMPLATES {CSP_CUSTOM, CSP_LAYER, CSP_LAYER_COL, CSP_BLOCK, CSP_CYLINDER, CSP_CONE, CSP_RCONE, CSP_CONE2, CSP_CONE3, CSP_CONE5, CSP_RCONE5, CSP_RWALL, CSP_WAVES_2D, CSP_WAVY_NS_LAYER, CSP_TRIANGLES, CSP_SRC_DISK, CSP_SRC_DISK_CEIL, CSP_SMILEY, CSP_FORSTEP};
#elif defined(MODEL_RIV)
enum CSP_TEMPLATES {CSP_CUSTOM, CSP_SLOPE};
#else
enum CSP_TEMPLATES {CSP_CUSTOM};
#endif

//initialization of available templates
void init_template_list()
{
  init_template(CSP_CUSTOM, "CUSTOM", "CUSTOM:\t\t\tno template", 0);
#if defined(MODEL_DUN) || defined(MODEL_SNO)
  init_template(CSP_LAYER, "LAYER", "LAYER(h=1.0):\t\t\tsand layer of height <h>", 1, 1.0);
#ifdef CELL_COLOR
  init_template(CSP_LAYER_COL, "LAYER_COL", "LAYER_COL(h=1.0, d=0.5):\tbinary sand layer (two colors of grain) of height <h> and density <d> for white color", 2, 1.0, 0.5);
#endif
  init_template(CSP_BLOCK, "BLOCK", "BLOCK(xmin, xmax, ymin, ymax, h):\tsand block of height <h>", 5, 30.0, 50.0, 50.0, 150.0, 20.0);
  init_template(CSP_CYLINDER, "CYLINDER", "CYLINDER(h=20, w=100):\tcylinder of height <h> and width <w>", 2, 20.0, 100.0);
  //init_template(CSP_CONE, "CONE", "CONE(h=30, w=100):\t\tcone of height <h> and width <w>", 2, 30.0, 100.0);
  init_template(CSP_CONE, "CONE", "CONE(h=30, w=100, x=L/2, y=D/2):\t\tcone of height <h> and width <w> centered on (x,y)", 4, 30.0, 100.0, 0.0, 0.0);
  //init_template(CSP_RCONE, "RCONE", "RCONE(rh=0.66, rw=0.4):\tcone of relative height <rh>*H and relative width <rw>*Min(L,D)", 2, 0.66, 0.4);
  init_template(CSP_RCONE, "RCONE", "RCONE(rh=0.66, rw=0.4, rx=0.5, ry=0.5):\tcone of relative height <rh>*H, relative width <rw>*Min(L,D) and centered on relative position (<rx>*L,<ry>*D)", 4, 0.66, 0.4, 0.5, 0.5);
#ifdef CELL_COLOR
  init_template(CSP_CONE2, "CONE2", "CONE2:\t\t2 cones with different colors and sizes", 0);
#endif
  init_template(CSP_CONE3, "CONE3", "CONE3(h=40, w=126):\t\t3 cones of height <h> and width <w>", 2, 40.0, 126.0);
  init_template(CSP_CONE5, "CONE5", "CONE5(h=40, w=126):\t\t5 cones of height <h> and width <w>", 2, 40.0, 126.0);
  init_template(CSP_RCONE5, "RCONE5", "RCONE5(rh=0.66, rw=0.2):\t5 cones of relative height <rh>*H and relative width <rw>*Min(L,D)", 2, 0.66, 0.2);
  init_template(CSP_RWALL, "RWALL", "RWALL(rh=0.5, rx=0.3):\twall of relative height <rh>*H and relative position <rx>*L", 2, 0.5, 0.3);
  init_template(CSP_WAVES_2D, "WAVES_2D", "WAVES_2D(per_min=10, per_max=60, ufd=0, amp=4, mh=30):\t2d sand waves with increasing periods from <per_min> to <per_max>, amplitude <amp> and mean height <mh> for stability analysis\n\t\t\t\tIt has uniform frequency distribution when <ufd> flag is set", 5, 10.0, 60.0, 0.0, 4.0, 30.0);
  init_template(CSP_WAVY_NS_LAYER, "WAVY_NS_LAYER", "WAVY_NS_LAYER(per=40, amp=15, h=5):\tnorth-south sand layer with wavy east boundary of period <per>, amplitude <amp> and height <h> for stability analysis", 3, 40.0, 15.0, 5.0);
  init_template(CSP_TRIANGLES, "TRIANGLES", "TRIANGLES(n=8, h=10, mh=20):\t<n> east-west periodic triangles of local height <h> and mean height <mh>", 3, 8.0, 10.0, 20.0);
  init_template(CSP_SRC_DISK, "SRC_DISK", "SRC_DISK(w=40, x=100, y=D/2):\tcircular source of sand on the ground with width <w> and centered on (x,y)", 3, 40.0, 100.0, 0.0);
  init_template(CSP_SRC_DISK_CEIL, "SRC_DISK_CEIL", "SRC_DISK_CEIL(w=40, x=100, y=D/2):\tcircular source of sand in the ceiling with width <w> and centered on (x,y)", 3, 40.0, 100.0, 0.0);
  init_template(CSP_SMILEY, "SMILEY", "SMILEY:\t\t\t:-)", 0);
  init_template(CSP_FORSTEP, "FORSTEP", "FORSTEP(h=0.15, hg=0.1, n=2):\tset of integer <n> solid forward-facing steps each <h> cells high, coated with a thickness <hg> of grains", 3, 0.15, 0.1, 2.0);
#elif defined(MODEL_RIV)
  init_template(CSP_SLOPE, "SLOPE", "SLOPE(r=1.0):\t\tsloping ground of ratio <r>", 1, 1.0);
#endif
  //LogPrintf("init_template_list: %d templates\n", nb_templates);
}

int32_t parse_template()
{
  int8_t *ptr;
  int32_t i;
  int32_t err=0;

  //read template name
  ptr = strtok(csp_template_str,"(");
  LogPrintf("csp_template.name = %s\n", ptr);

  //find and copy the full template
  for (i=0; (i < nb_templates) && strcmp(ptr,t_templates[i].name); i++);
  if (i < nb_templates){
    memcpy(&csp_template, &t_templates[i], sizeof(CSP_Template));
    LogPrintf("csp_template.type = %d\n", csp_template.type);
    AllocMemory(csp_template.args, float, csp_template.nb_args);
    memcpy(csp_template.args, t_templates[i].args, csp_template.nb_args*sizeof(float));
  }
  else{
    ErrPrintf("ERROR: Incorrect value for Csp_template, %s\n", csp_template_str);
    exit(-2);
  }

  //read template optional arguments, if any
  ptr = strtok(NULL,",)");
  i = 0;
  while (ptr && (i < csp_template.nb_args)){
    //LogPrintf("csp_template.args[%d] = %s\n", i, ptr);
    csp_template.args[i++] = read_float(ptr, &err);
    if (err){
      ErrPrintf("ERROR: bad argument \"%s\" for template %s\n", ptr, csp_template.name);
      exit(-1);
    }
    ptr = strtok(NULL,",)");
  }
  if (ptr){
    ErrPrintf("ERROR: Too many arguments in Csp_template (%s)\n", ptr);
    exit(-2);
  }
  for (i=0; i < csp_template.nb_args; i++){
    LogPrintf("csp_template.args[%d] = %f\n", i, csp_template.args[i]);
  }

  return i;
}


void genesis_parse()
{
  if (model_name && strcmp(model_name, MOD_NAME)){
    WarnPrintf("ERROR: Wrong model name %s in parameter file (caught in genesis.c). Check model in defs.h and rebuild all.\n", model_name);
    exit(-2);
  }

  if (csp_template_str){
    parse_template();
  }

  if (boundary_str){
    if (!strcmp(boundary_str, "PERIODIC")){
#ifdef CYCLAGE_HOR
      boundary = BC_PERIODIC;
#else
      ErrPrintf("Incorrect value for Boundary : %s, please define CYCLAGE_HOR in defs.h and recompile\n", boundary_str);
      exit(-2);
#endif
    }
    else if (!strcmp(boundary_str, "OPEN")){
      boundary = BC_OPEN;
    }
    else if (!strcmp(boundary_str, "OUT")){
      boundary = BC_OUT;
    }
    else if (!strcmp(boundary_str, "CLOSE")){
      boundary = BC_CLOSE;
    }
    else if (!strcmp(boundary_str, "REINJECTION")){
      boundary = BC_REINJECTION;
    }
    else{
      ErrPrintf("Incorrect value for Boundary : %s\n", boundary_str);
      exit(-2);
    }
    //pbc_mode = (boundary == BC_PERIODIC);
  }
}

void genesis()
{
  int32_t i, j, k, Hd2, Ld2, Dd2, Dd3, Dd4, Hd4, Hd3;
  int32_t Ld3, Ld34, Ld4, Ld8;
  Cell *aux;

  HL = H*L;
  HLD = HL*D;

  // allocation of the cellular space
  AllocMemoryPrint("TE",TE, Cell, HLD);
  ResetMemory(TE, Cell, HLD);

  // relative lengths
  Ld2 = (int) L/2;
  Ld3 = (int) L/3;
  Ld4 = (int) L/4;
  Ld8 = (int) L/8;
  Ld34= (int) 3*L/4;
  Dd2 = (int) D/2;
  Dd3 = (int) D/3;
  Dd4 = (int) D/4;
  Hd2 = (int) H/2;
  Hd3 = (int) H/3;
  Hd4 = (int) H/4;

  aux = TE;


/*****************************************************************************/
/********************************* DUN model *********************************/
/*****************************************************************************/
#if defined(MODEL_DUN) || defined(MODEL_SNO)
  float di, dj, dk;
  float Ldx, Ldy, Ldz;
  float pente = 2.0;
  float h, w, rc, n, hg;
  int32_t xmin, xmax, ymin, ymax;
  float x, y, z;

//// normal mode ////

  //Ldx= Ld2 + 0.5; //Ld4; //10+Hd2;
  Ldx= Ld2-0.5; //Ld2-0.5; //Ld3;//Ld4;
  Ldy= /*0.66**/Dd2-0.5; //Ld4; //(int) L/2;
  Ldz= Hd2;

  int32_t ll1 = (int) L*0.2; //L*0.15; //L*0.5;//L*0.1; //10; //0.1 //0.22; //Ld3
  int32_t ll2 = (int) L*0.2; //0.22; //Ld3
  int32_t l0 = (int) L/8; //L*0.1; //(L-ll1-ll2)*0.4;//25; //0.45
  int32_t l1 = (int) l0+ll1*1.3; //pour 2 triangles
  int32_t l2 = (int)(L-ll1-ll2)*0.25; //0.45
 //int32_t l1 = (int)(L-l0-l2-ll1-ll2);
  int32_t hh = (int) (H/6);//(100/3); //8+(100/3); //(H/3);//(60/2); //H/6; //H/2;
  int32_t H0 = H-2; //H-8
  int32_t hh1 = hh/1.5;
  int32_t hh2 = hh;
  int32_t lc = D; //1; //400; //Ld4; //Ld2; //L/10; //Ld4; //largeur du couloir
  int32_t lcone = lc*0.6; //lc*0.6; //lc/2; //lc*0.6; //lc*2/3;  //Ld2*0.8; //largeur du cone
  int32_t jcone = H0-lcone/2; //sommet du cone
  int32_t hcone = Hd3; //Hd2;//2*Hd3; //H-10; //hauteur du cone
 //int32_t h0 = (int) (H-hh)/2;
 //int32_t h2 = H-h0;
  int32_t rcyl = Ld2-10; //rayon cylindre
  int32_t hgrain = Hd2; //epaisseur de grains
  float periode = 10.0; //periode des ondulations
#ifdef CELL_COLOR
  float density = 0.5; //density of white grains
#endif

  // initialization of the cellular space
  for(k=0; k < D; k++){ // profondeur
    //periode = 10.0+(k-(Dd2-lc/2))/2.0/4.0;
    periode = 10.0+(k-(Dd2-lc/2))/2.0/2.0;
    //periode += 0.5*drand48();
    for(j=0; j < H; j++){ // hauteur
      for(i=0; i < L; i++, aux++){ //largeur
        di = i - Ldx;
        dj = j - Ldz;
        dk = k - Ldy; // + Ld4;
        rc = 0.5*lcone*(j-(H-hcone))/(hcone-1);

        /*** EAUC cells ***/
        aux->celltype = EAUC; //default state

        /*** GR cells ***/
        if (csp_template.type == CSP_LAYER){
          //CSP_LAYER: sand layer
          //format: LAYER(h)
          hh = (int)csp_template.args[0];
          float frac = fmodf(csp_template.args[0],1.0);
          if (j > H-2-hh) aux->celltype = GR;
          if ((j == H-2-hh) && frac){ //h is not an integer
            if (drand48() < frac) aux->celltype = GR;
          }
        }
#ifdef CELL_COLOR
        else if (csp_template.type == CSP_LAYER_COL){
          //CSP_LAYER_COL: binary sand layer
          //format: LAYER_COL(h,d)
          hh = (int)csp_template.args[0];
          density = csp_template.args[1];
          if (density>1){ErrPrintf("ERROR: bad parameter in LAYER_COL template, <d> must be <= 1 (%f)\n", density); exit(-1);}
          float frac = fmodf(csp_template.args[0],1.0);
          if (j > H-2-hh) aux->celltype = GR;
          if ((j == H-2-hh) && frac){ //h is not an integer
            if (drand48() < frac) aux->celltype = GR;
          }
        }
#endif
        else if (csp_template.type == CSP_BLOCK){
          //CSP_BLOCK: sand block
          //format: BLOCK(xmin, xmax, ymin, ymax, h)
          xmin = (int)csp_template.args[0];
          xmax = (int)csp_template.args[1];
          ymin = (int)csp_template.args[2];
          ymax = (int)csp_template.args[3];
          hh = (int)csp_template.args[4];
          if ((i >= xmin) && (i < xmax) && (k >= ymin) && (k < ymax) && (j > H-2-hh)){
            aux->celltype = GR;
            //aux->celltype = GRV;
          }
        }
        else if (csp_template.type == CSP_CYLINDER){
          //CSP_CYLINDER: sand cylinder
          //format: CYLINDER(h, w)
          h = csp_template.args[0];
          w = csp_template.args[1];
          if ((j >= (H-h-1)) && (di*di + dk*dk < 0.25*w*w)) aux->celltype = GR;
        }
        else if (csp_template.type == CSP_CONE){
          //CSP_CONE: sand cone
          //format: CONE(h, w, x, y)
          hcone = csp_template.args[0];
          lcone = csp_template.args[1];
          x = csp_template.args[2];
          y = csp_template.args[3];
          if (x==0) x=Ldx;
          if (y==0) y=Ldy;
          rc = 0.5*lcone*(j-(H-hcone))/(hcone-1);
          if ((j > (H-hcone)) && ((i-x)*(i-x) + (k-y)*(k-y) < rc*rc)) aux->celltype = GR;
        }
        else if (csp_template.type == CSP_RCONE){
          //CSP_RCONE: sand cone with relative dimensions and position
          //format: RCONE(rh, rw, rx, ry)
          hcone = csp_template.args[0]*H;
          lcone = csp_template.args[1]*Min(L,D);
          x = csp_template.args[2]*L - 0.5;
          y = csp_template.args[3]*D - 0.5;
          rc = 0.5*lcone*(j-(H-hcone))/(hcone-1);
          if ((j > (H-hcone)) && ((i-x)*(i-x) + (k-y)*(k-y) < rc*rc)) aux->celltype = GR;
          //if (j > H-3) aux->celltype = GR; //uniform layer
        }
        else if (csp_template.type == CSP_CONE2){
          //CSP_CONE2: 2 sand cones for dune collisions, with different colors
          //format: CONE2
          hcone = H/2;
          float lcone0 = D/2;
          float lcone1 = D/3;
          float rc0 = 0.5*lcone0*(j-(H-hcone))/(hcone-1);
          float rc1 = 0.5*lcone1*(j-(H-hcone))/(hcone-1);
          float di0 = i - Ld3;
          float dk0 = k - D*3/5;
          float di1 = i - L/8;
          float dk1 = k - D/3;
          if ((j > (H-hcone)) && (di0*di0 + dk0*dk0 < rc0*rc0)) aux->celltype = GR;
          if ((j > (H-hcone)) && (di1*di1 + dk1*dk1 < rc1*rc1)){
            aux->celltype = GR;
#ifdef CELL_COLOR
            aux->color = 1;
#endif
          }
        }
        else if (csp_template.type == CSP_CONE3){
          //CSP_CONE3: 3 sand cones with fixed dimensions
          //format: CONE3(h, w)
          hcone = csp_template.args[0];
          lcone = csp_template.args[1];
          rc = 0.5*lcone*(j-(H-hcone))/(hcone-1);
          float di0 = i - Ldx - 0.5*Ldx*cos(0);
          float dk0 = k - Ldy - 0.5*Ldy*sin(0);
          float di2 = i - Ldx - 0.5*Ldx*cos(2*PI/3);
          float dk2 = k - Ldy - 0.5*Ldy*sin(2*PI/3);
          float di4 = i - Ldx - 0.5*Ldx*cos(4*PI/3);
          float dk4 = k - Ldy - 0.5*Ldy*sin(4*PI/3);
          if ((j > (H-hcone)) && (di0*di0 + dk0*dk0 < rc*rc)) aux->celltype = GR;
          if ((j > (H-hcone)) && (di2*di2 + dk2*dk2 < rc*rc)) aux->celltype = GR;
          if ((j > (H-hcone)) && (di4*di4 + dk4*dk4 < rc*rc)) aux->celltype = GR;
        }
        else if ((csp_template.type == CSP_CONE5) || (csp_template.type == CSP_RCONE5)){
          if (csp_template.type == CSP_CONE5){
            //CSP_CONE5: 5 sand cones with fixed dimensions
            //format: CONE5(h, w)
            hcone = csp_template.args[0];
            lcone = csp_template.args[1];
          }
          else{
            //CSP_RCONE5: 5 sand cones with relative dimensions
            //format: RCONE5(rh, rw)
            hcone = csp_template.args[0]*H;
            lcone = csp_template.args[1]*Min(L,D);
          }
          rc = 0.5*lcone*(j-(H-hcone))/(hcone-1);
          float di1 = i - Ldx - 0.6*Ldx*cos(PI/5);
          float dk1 = k - Ldy - 0.6*Ldy*sin(PI/5);
          float di3 = i - Ldx - 0.6*Ldx*cos(3*PI/5);
          float dk3 = k - Ldy - 0.6*Ldy*sin(3*PI/5);
          float di5 = i - Ldx - 0.6*Ldx*cos(5*PI/5);
          float dk5 = k - Ldy - 0.6*Ldy*sin(5*PI/5);
          float di7 = i - Ldx - 0.6*Ldx*cos(7*PI/5);
          float dk7 = k - Ldy - 0.6*Ldy*sin(7*PI/5);
          float di9 = i - Ldx - 0.6*Ldx*cos(9*PI/5);
          float dk9 = k - Ldy - 0.6*Ldy*sin(9*PI/5);
          if ((j > (H-hcone)) && (di1*di1 + dk1*dk1 < rc*rc)) aux->celltype = GR;
          if ((j > (H-hcone)) && (di3*di3 + dk3*dk3 < rc*rc)) aux->celltype = GR;
          if ((j > (H-hcone)) && (di5*di5 + dk5*dk5 < rc*rc)) aux->celltype = GR;
          if ((j > (H-hcone)) && (di7*di7 + dk7*dk7 < rc*rc)) aux->celltype = GR;
          if ((j > (H-hcone)) && (di9*di9 + dk9*dk9 < rc*rc)) aux->celltype = GR;
          //if ((j > (H-hcone)) && (di*di + dk*dk < rc*rc)) aux->celltype = GR; //central cone
          //if (j > H-3) aux->celltype = GR; //uniform layer
        }
        else if (csp_template.type == CSP_RWALL){
          //CSP_RWALL(rh,rx): north-south wall or DUM cells
          int32_t h = csp_template.args[0]*H;
          int32_t x = csp_template.args[1]*L;
          if ((i>=x-1) && (i<=x+1) && (j>=H-h)) aux->celltype = DUM;
        }
        else if (csp_template.type == CSP_WAVES_2D){
          //CSP_WAVES_2D: 2d sand waves with increasing periods for stability analysis
          float per_min = csp_template.args[0];
          float per_max = csp_template.args[1];
          float ufd_flag = csp_template.args[2];
          float amp = csp_template.args[3];
          int32_t mh = (int)csp_template.args[4];
          if (per_min<=0){ErrPrintf("ERROR: bad parameter in WAVES_2D template, <per_min> must be positive (%f)\n", per_min); exit(-1);}
          if (ufd_flag){
            // uniform frequency distribution
            float fmin = 1.0/per_max;
            float fmax = 1.0/per_min;
            float f = fmax + (fmin-fmax)*k/D;
            periode = 1.0/f;
          }
          else
            periode = per_min + (per_max-per_min)*k/D;
          if ((i==0) && (j==0) && (k%2==0)) LogPrintf("k = %d, periode = %f\n", k, periode);
          if (j > (H-mh)-0.5*amp*sin(i*2*PI/periode)) aux->celltype = GR; //variable ondulations
          if ((k & 1) || (j > H0) || (j==0)) aux->celltype = DUM; //ground + ceiling + walls east-west
        }
        else if (csp_template.type == CSP_WAVY_NS_LAYER){
          //CSP_WAVY_NS_LAYER: north-south sand layer with wavy east boundary for stability analysis
          float per = csp_template.args[0];
          float amp = csp_template.args[1];
          int32_t hh = (int)csp_template.args[2];
          int32_t i_wave = (int)(40 + 2 + 0.5*amp*(1+sin(k*2*PI/per)));
          if ((i >= 40) && (i <= i_wave) && (j <= H0) && (j > H0-hh)) aux->celltype = GR; //rangee de sable avec front de forme sinusoidale
        }
        else if (csp_template.type == CSP_TRIANGLES){
          //CSP_TRIANGLES: east-west periodic triangles
          float lx0 = 0.5*L/csp_template.args[0];
          float lx1 = fmodf(i, lx0*2);
          hh=csp_template.args[1];
          int32_t mh=csp_template.args[2];
          if (((lx1<lx0) && (j>H0-mh+hh/2-(hh*lx1/lx0))) || ((lx1>=lx0) && (j>H0-mh+hh/2-(hh*(2*lx0-lx1)/lx0)))) aux->celltype = GR;
        }
        else if ((csp_template.type == CSP_SRC_DISK) || (csp_template.type == CSP_SRC_DISK_CEIL)){
          //CSP_SRC_DISK: circular source of sand
          w = csp_template.args[0];
          x = csp_template.args[1];
          y = csp_template.args[2];
          if (y == 0) y = D*0.5;
          z = (csp_template.type == CSP_SRC_DISK) ? H0+1 : 0 ;
          if ((j == z) && ((i-x)*(i-x) + (k-y)*(k-y) <= w*w/4.0)) aux->celltype = IN;
        }
        else if (csp_template.type == CSP_SMILEY){
          //CSP_SMILEY: :-)
          if (j<H/2 && ((cos(PI*(i-L/2)/L/2)*cos(PI*(k-D/2)/D/2)<0.8 && cos(PI*(i-L/2)/L/2)*cos(PI*(k-D/2)/D/2)>0.75)
             || (cos(2*PI*(i-L/2)/L)<0.4 && sin(PI/4+2*PI*(k-D/2)/D)<0.4 && cos(2*PI*(i-L/2)/L)>0.2 && sin(PI/4+2*PI*(k-D/2)/D)>0.2 && k<L/2)
             || (cos(PI*(i-L/2)/L)*cos(PI*(k-D/2)/D)<0.8 && cos(PI*(i-L/2)/L)*cos(PI*(k-D/2)/D)>0.75&&k>3*D/5)))
               aux->celltype = GR;
        }
        else if (csp_template.type == CSP_FORSTEP){
	  //CSP_FORSTEP: set of forward-facing steps
          h = csp_template.args[0]; // height of each step
	  n = csp_template.args[1]; // number of steps
	  hg = csp_template.args[2]; // thickness of space that is granular
	  float b = h; // thickness of solid layer below steps
	  for( float step = -1.0; step < n; step = step + 1.0 ){
            float m_s;
            m_s = n*h/L;
#ifdef GRV  // if cohesive grain (GRV) defined, make cohesive steps
	    if ( (i>=L*step/n) && (i<L*(step+1.0)/n) && 
	      (j>=(H-b-(h-m_s*(i-L*step/n)))) ) aux->celltype = GRV;
#else       // else make solid unerodible (DUM) steps
	    if ( (i>=L*step/n) && (i<L*(step+1.0)/m) &&
	      (j>=(H-b-(h-m_s*(i-L*step/n)))) ) aux->celltype = DUM; 
#endif
	    if (j<= hg) aux->celltype = GR;
	  }
	}
        else{
          //CSP_CUSTOM
          //if (0) //pas de grain
          //if ((di*di + dk*dk < (j-(H-hcone))*(j-(H-hcone))) && (j > (H-hcone))) //cone de hauteur hcone
          //if ((di*di + dk*dk < 4*(j-2*Hd3)*(j-2*Hd3)) && (j > 2*Hd3)) //cone de hauteur Hd3 et de rayon 2*Hd3
          //if ((di*di + dk*dk < 2.5*(j-Hd3)*(j-Hd3)) && (j > Hd3)) //cone de hauteur 2*Hd3 pour les dunes etoiles
          //if ((di*di + dk*dk < (j-jcone)*(j-jcone)) && (j > Hd2)) //cone tronque
          //if (((di*di + dk*dk < (j-jcone)*(j-jcone)) && (j > jcone)) || (j > H-10)) //cone + couche uniforme
          //if ((j > 12) && (di*di + dk*dk < 400)) // disque
          //if ((j > H*0.65) && (fabs(di) < Ld8) && (fabs(dk) < Ld8)) // parallelepipede rectangle
          //if ((j > Hd2) && ((di+dk)*(di+dk) < 400) && (dk*dk < 400)) // parallelepipede 'tordu'
          //if (j + (i/20)> H-10) //pente EO
          //if (j + k*pente > H-20) //pente NS
          //if (j + i > H) //pente EO
          if (j > H-10) //couche uniforme
          //if (j>H-8-(drand48()-0.90)) //couche aleaoire
          //if ((j > H-4) && (k>=Dd4) && (k<D-Dd4)) //couche uniforme est-ouest
          //if ((j > H-5) && (di*di + dk*dk < rcyl*rcyl)) //couche uniforme circulaire
          //if (/*(i>1) &&*/ (j > H-Hd2*drand48())) //couche non-uniforme
          //if ((j > H-2) || ((j == H-2) && i < L*0.5 && (drand48()>0.75/*0.5*/))) //couche non-uniforme
          //if (j > Hd2-5*cos(10*i*PI/L)) //ondulations
          //if (j > Hd2-2*sin(i*2*PI/periode)) //ondulations variables
          //if (((i>=l0) && (i<l0+1/*l0*/) && (j>H-hh)) || (j>=H-1)) //rectangle
          //if (((i>=l0) && (i<l0+1) && (j>H-hh)) || (j< (H-hh-10)*(Ld2-k)/lc) || (j>=H-1)) //rectangle avec toit incline
          //if (((i>=l0) && (i<l0+ll1) && ((i-l0)-((H0-j/*+drand48()*5*/)*ll1/hh) > 0))
          //      || (j>=H-1)) //triangle
          //if (((i>=l0) && (i<l0+ll1) && ((i-l0)-((H0-j/*+drand48()*5*/)*ll1/hh) > 0))
          //           || ((i>=l0+ll1) && ((i-l0-ll1)-(hh-(H0-j))/5 < 0))
          //           || (j>=H0)) //triangle + triangle inverse
          //if (((i>=l0) && (i<l0+ll1) && ((i-l0)-((H-j)*ll1/hh1) > 0))
          //      || ((i>=l1) && (i<l1+ll2) && ((i-l1)-((H-j)*ll2/hh2) > 0))
          //      || (j>=H-1)) //2 triangles
          /*if (((i>=l0) && (i<=l0+ll1) && ((i-l0)-((H-j)*ll1/hh) > 0)) || (j< H*((k<Ld2) ? Ld2-k : k-Ld2)/Ld2)
            || (j>=H-1)) //triangle avec toit incline*/
                          /*else if (((i>=l0) && (i<=l0+ll1) && ((i-l0)-((H-j)*ll1/hh) > 0)) || (j< Hd2*((k<Ld2) ? Ld2-k : k-Ld2)/ll)
            || (j>=H-1)) //triangle avec toit en biseau*/
          //if (((i>=l0) && (i<=l0+ll1) && ((i-l0)-((H-j)*ll1/hh) > 0)) || (j< (H-hh-10)*(Ld2-k)/lc) || (j>=H-1)) //triangle avec toit incline
          /*if (((i>=l0) && (i<l0+ll1) && ((i-l0)-((H-j)*ll1/hh) > 0))
            || ((i>=l0+ll1) && (i<L-ll2-l2) && (j>=hh))
            || ((i>=L-ll2-l2) && (i<L-l2) && ((L-l2-i)-((H-j)*ll2/hh) > 0))
            || (j>=H-1)) //trapeze*/
          //if (di*di + dj*dj + dk*dk < Hd4*Hd4) // sphere
          //if ((di*di + dj*dj + dk*dk < Hd4*Hd4) || ((di+50)*(di+50) + (dj+30)*(dj+30) + dk*dk < Hd4*Hd4/4)) // 2 spheres
          //if (j >H0-3*(1+cos(10*i*k*PI/L))) //test mh
          //if (j >H0-3*(1+cos(10*i*k*PI/L/D))) //test mh
          //if (j >H0-3*(0+cos(10*i*PI/L)-cos(10*k*PI/D))) // champs de bosses
          //if ((((l1=i%(l0*2))<l0) && (j>H0-Hd3+hh/2-(hh*l1/l0))) || ((l1>=l0) && (j>H0-Hd3+hh/2-(hh*(2*l0-l1)/l0)))) // barres triangulaires periodiques
          {
            aux->celltype = GR;
          }
        }

        /*** DUM cells ***/
        //if (0)
        //if (j > H-2) //plateau
        if ((j > H0) || (j==0)) //plateau + plafond
        //if ((j > H0-round(drand48())) || (j==0)) //plateau rugueux + plafond
        //if ((j > H0) || (j==0) || ((i==Ld3) && (j>=Hd2))) //plateau + plafond + mur
        //if ((j > H0) || (j==0) || (k>=D-5) || (k<5)) //plateau + plafond + couloir
        //if ((j > H0) || (j<4)) //plateau + plafond epais (4)
        //if ((k < Dd2-lc/2) || (k > Dd2+lc/2) || (k & 1) || (j > H0) || (j==0)) //plateau + plafond + couloir + peigne
        //if ((j > H0-5*(1+cos(10*i*PI/L))) || (j==0)) //plateau ondule + plafond
        //if ((j > H0-5*(1+cos(10*i*PI/L))) || (j < 4)) //plateau ondule + plafond epais (4)
        //if (((i == L-1) && (j>H-L))) //couloir + mur
        //if ((j > H-2) && (di*di + dk*dk < Ld2*Ld2)) // disque + plateau + plafond
        //if (((j > hgrain) && (di*di + dk*dk >= rcyl*rcyl)) || (j > H0) || (j==0)) // plateau + trou circulaire + plafond
        //if (j + (i/20)> H-2) //pente EO
        //if ((j + (i/20)> H-2) || (j==0)) //pente EO + plafond
        //if (j + k*pente > H+20) //pente NS
        //if ((j > H0) || (j==0) || ((j > (H-hcone)) && (di*di+dk*dk < rc*rc))) //plateau + cone + plafond
          if (aux->celltype != IN) aux->celltype   = DUM;

        /*** IN cells ***/
        if (boundary == BC_REINJECTION){
          //if ((k == Ld2) && (i == Ld2) && (j == 0)) //source de grains en haut au centre
          if (/*(node == 0) &&*/ (i == 5) && (j == 0)) //source de grains en haut a gauche
          //if ((j == 0) && (di*di + dk*dk < 5*5)) ////source de grains en disque
          //if (((i == 40) || (i == 41)) && (j > H0)) //source de grains nord-sud au sol
          /*float Ldx1 = 100-0.5;//-Ld3*cos(PI/4);
          float Ldy1 = Dd2-0.5;//-Ld3*sin(PI/4);
          if ((j > H0) && ((i-Ldx1)*(i-Ldx1) + (k-Ldy1)*(k-Ldy1) <= 20*20)) ////source de grains en disque au sol
            aux->celltype = IN;*/
          /*float Ldx1 = Ld2-0.5-Ld3*cos(PI/3)+(L/10)*sin(PI/3);
          float Ldy1 = Dd2-0.5-Ld3*sin(PI/3)-(L/10)*cos(PI/3);
          float Ldx2 = Ld2-0.5-Ld3*cos(PI/3)-(L/10)*sin(PI/3);
          float Ldy2 = Dd2-0.5-Ld3*sin(PI/3)+(L/10)*cos(PI/3);
          if ((j > H0) && (((i-Ldx1)*(i-Ldx1) + (k-Ldy1)*(k-Ldy1) <= 25)
                       || ((i-Ldx2)*(i-Ldx2) + (k-Ldy2)*(k-Ldy2) <= 25))) ////2 sources de grains en disque au sol*/
          //if ((j > H0) && ((di*di + dk*dk >= Ld3*Ld3) && ((di+cos(PI/3))*(di+cos(PI/3)) + (dk+sin(PI/3))*(dk+sin(PI/3)) < Ld3*Ld3))) ////source de grains en demi-cercle au sol, sur la table tournante, avec un angle initial de Pi/3
          //if ((j == 0) && ((di*di + dk*dk >= rcyl*rcyl) && ((di+1)*(di+1) + dk*dk < rcyl*rcyl))) ////source de grains en demi-cercle
          //if ((j == 0) && (dk > -5) && (dk < 5) && ((di*di + dk*dk >= rcyl*rcyl) && ((di+1)*(di+1) + dk*dk < rcyl*rcyl))) ////source de grains en arc de cercle
          //if ((j == 0) && (((k >= i) && (i+k <= L) && (di*di + dk*dk >= rcyl*rcyl) && ((di+1)*(di+1) + dk*dk < rcyl*rcyl))
          //     || (((k <= i) || (i+k >= L)) && (i+1==round(Ld2-rcyl/sqrt(2)))))) //source de grains en quart de cercle + segments de droite
            aux->celltype = IN;
        }

        /*** OUT cells ***/
        //if ((j > H0) && (i>=L-5) && (i<L-1)) ////sortie de grains au sol pres du bord est
        //if ((j-1 == hgrain) && (di*di + dk*dk >= rcyl*rcyl)) ////sortie de grains en cercle
        //if ((j-1 == hgrain) && ((di>0) && ((di-1)*(di-1) + (dk-sgn(dk))*(dk-sgn(dk)) >= rcyl*rcyl))) ////sortie de grains en demi-cercle
        //  aux->celltype = OUT;
        /*if ((j==H-1) && (k>=Dd3-2) && (k<Dd3+2)) ////sortie de grains entre les couloirs
          aux->celltype = OUT;
        if ((j==H-1) && (k>=2*Dd3-2) && (k<2*Dd3+2)) ////sortie de grains entre les couloirs
          aux->celltype = OUT;*/

#ifdef CELL_COLOR
        /*** color of cells ***/
        if ((aux->celltype == GR) && (csp_template.type != CSP_CONE2)){
          //aux->color = 1; //unformly colored
          aux->color = (drand48()<density) ? 1 : 0; //randomly colored
          //aux->color = (k > Dd2) ? 1 : 0; //colored area to the north
          //aux->color = (k == Dd2+1) && (j == H-20) && (i == 40) ? 1 : 0; //one colored cell for TRACE_PAR_COL
          //aux->color = (k == Dd2+1) && (i == 40) ? H-1-j : 0; //column of colored cells for TRACE_PAR_COL
        }
#endif
#ifdef USE_VEGETATION
        /*** vegetated cells ***/
        //if ((i>=150) && (j==H0)) aux->celltype = VEG;
        //if (/*(i>=150) &&*/ (j>H0)) aux->celltype = VEG;
#endif
      }
    }
  }


//// precomputed mode (height map) ////
/*
  float* surf;
  int32_t lx = 800; //799;
  int32_t ly = 400; //399;
  int32_t cx = L/2;  //lx/2;
  int32_t cy = D/2; //L*0.4;  //L/3; //L/2;
  int32_t ly1 = (L-ly)/2;
  int32_t ly2 = ly1+ly;
  float alpha = 0; //-PI/6.0;

  surf = load_surface("data/scale_dune_4.txt", lx, ly);
  //surf = load_surface("data/trunc_dune_2.txt", lx, ly);

  for(k=0; k < D; k++){ // profondeur
    for(j=0; j < H; j++){ // hauteur
      for(i=0; i < L; i++, aux++){ //largeur
        // rotation
        int32_t i0 = (i-cx)*cos(alpha) + (k-cy)*sin(alpha) + cx;
        int32_t k0 = (i-cx)*sin(alpha) - (k-cy)*cos(alpha) + cy;
        if (0)
        //if (k>L*0.8)
        //if ((k<L*0.4) || (k>L*0.6))
        //if ((k<Ld2) || (k>Ld2))
          aux->celltype = DUM;
        else if ((k0>=ly1) && (k0<ly2) && (i0>=0) && (i0<lx) && (j > H-surf[i0+(k0-ly1)*lx]))
          aux->celltype = GR;
        else
          //aux->celltype = EAUC;
          //aux->celltype = (drand48() < 0.92) ? EAUC : EAUT;
          //aux->celltype = (j>Hd2) ? EAUC : EAUT;
          aux->celltype = (j>0)||(i>2) ? EAUC : EAUT;
       }
    }
  }
*/
#endif

/*****************************************************************************/
/********************************* AVA model *********************************/
/*****************************************************************************/
#ifdef MODEL_AVA
  int32_t hc = H*0.4;
  float Ldx, Ldy;
  float di, dk;
  float tg = tan(PI/6);

  Ldx= Ld2;//5 + 0.5; //Ld4; //10+Hd2;
  Ldy= Dd2;//5 + 0.5; //Ld4; //(int) L/2;

  for(k=0; k < D; k++){ // profondeur
    for(j=0; j < H; j++){ // hauteur
      for(i=0; i < L; i++, aux++){ //largeur
        di = i - Ldx;
        dk = k - Ldy;
        if (di<0) di = -di;
        if (dk<0) dk = -dk;

        //if (0)
        if (j > H-2) //plateau
        //if ((j > H-2) || ((j >= H - 20) && (di*di + dk*dk < 200))) //plateau et socle
           aux->celltype   = DUM;
        //else if ((i==0) && (k==0) && (j==0)) //source de grains
        //else if ((i==Ld2) && (k==Ld2) && (j==0)) //source ponctuelle de grains
        else if ((di*di + dk*dk <= 25) && (j==0)) //source etendue de grains
        //else if ((((i-Ld3)*(i-Ld3) + dk*dk <= 25) || ((i-2*Ld3)*(i-2*Ld3) + dk*dk <= 25)) && (j==0)) //2 sources etendues de grains
        //else if ((j==0) && (di*di + dk*dk < 55*55) && ((di-1)*(di-1) + (dk-1)*(dk-1) >= 45*45)) //source circulaire de grains
          aux->celltype = IN;
        //else if ((di*di + dk*dk <= 25) && (j==H-2)) //sortie circulaire de grains
        //  aux->celltype = OUT;
        //else if ((j>H/4) && (di*di + dk*dk <= (Ld2-5)*(Ld2-5))) //cylindre
        //else if ((di*di + dk*dk < 2*(j-Hd3)*(j-Hd3)) && (j > Hd3)) //cone de hauteur 2*Hd3
        //else if ((di*di + dk*dk < (j-H+hc)*(j-H+hc)/(tg*tg)) && (j > H-hc)) //cone de hauteur hc et de pente tg
        //    aux->celltype   = DUM;
        else{
          aux->celltype   = AIR;
          //aux->celltype   = GR;
#ifdef CELL_COLOR
          //aux->color = (drand48()<0.5) ? 1 : 0;
#endif
        //if ((j == H-1) && (di*di + dk*dk >= (Ld2-5)*(Ld2-5))) ////sortie de grains en cercle
        //  aux->celltype = OUT;
        }
      }
    }
  }
#endif

/*****************************************************************************/
/********************************* CMB model *********************************/
/*****************************************************************************/
#ifdef MODEL_CMB
  for(k=0; k < D; k++) // profondeur
    for(j=0; j < H; j++) // hauteur
      for(i=0; i < L; i++, aux++){ //largeur
 //Modele classique
        if (j < Hd2)
        //else if ((j < Hd2) || ((j < Hd2+20) && (i > Ld2-10) && (i < Ld2+10)))
          aux->celltype   = PLUS;
        else
          aux->celltype   = MOINS;
        //aux->celltype   = (int) floor(drand48()*3);

/* Modele experimental 1 (2 plaques)*/
/*
        int32_t L1 = L/3;
        int32_t L2 = L*2/3;
        int32_t H1 = H/3;
        int32_t H2 = H*2/3;
        if ((i>L1) && (i<L2) && (k>L1) && (k<L2) && ((j<=H1) || (j>=H2)))
          aux->celltype = PLUS;
        else if ((i==L1) && (k>=L1) && (k<=L2) && ((j<=H1+5) || (j>=H2-5)))
          aux->celltype = DUM;
        else if ((i==L2) && (k>=L1) && (k<=L2) && ((j<=H1) || (j>=H2)))
          aux->celltype = DUM;
        else if ((k==L1) && (i>=L1) && (i<=L2) && ((j<=H1+5) || (j>=H2-5)))
          aux->celltype = DUM;
        else if ((k==L2) && (i>=L1) && (i<=L2) && ((j<=H1) || (j>=H2)))
          aux->celltype = DUM;
        else if (drand48() < 0.5)
          aux->celltype = ZERO;
        else
          aux->celltype = MOINS;
*/
// Modele experimental 2 (sillon diagonal)
 /*
        int32_t D = H/6;
        int32_t L1 = (L/4)-D;
        int32_t L2 = (L/4)+D;
        int32_t L3 = (3*L/4)-D;
        int32_t L4 = (3*L/4)+D;
        int32_t H1 = Hd2-D;
        int32_t H2 = Hd2+D;
        int32_t ii = i+k-Ld2; //fenetre courante
        if (ii<0) ii += L;
        if (ii>=L) ii -= L;
        if (ii<L1)
          aux->celltype = (j<H2) ? PLUS : MOINS;
        else if (ii<L2)
          aux->celltype = (j<H2-(ii-L1)) ? PLUS : MOINS;
        else if (ii<L3)
          aux->celltype = (j<H1) ? PLUS : MOINS;
        else if (ii<L4)
          aux->celltype = (j<H1+(ii-L3)) ? PLUS : MOINS;
        else
          aux->celltype = (j<H2) ? PLUS : MOINS;

  */
/* Modele experimental 3 (1 faille oblique)*/
/*
        int32_t L1 = (L-Hd2)*0.9;
        int32_t L2 = (L-Hd2)*1.1;
        int32_t ii = i+j; // fenetre courante
        if ((ii>L1) && (ii<L2))
          aux->celltype = MOINS;
        else
          aux->celltype = PLUS;
*/

/* Modele experimental 4 (2 fractures auto-affines)*/
/*
        //premiere fracture oblique
        static float* surf1 = NULL;
        static float* surf2 = NULL;
        if (!surf1) surf1 = load_surface("surf256.dat", 256, 256);
        if (!surf2) surf2 = load_surface("surf256b.dat", 256, 256);
        float alpha = 70*PI/180;
        float amp = 1.0;

        int32_t H1 = Hd2-5;
        int32_t H2 = Hd2+5;

        // rotation
        int32_t i0 = (i-Ld2)*cos(alpha) + (j-Hd2)*sin(alpha) + Ld2;
        int32_t j0 = (i-Ld2)*sin(alpha) - (j-Hd2)*cos(alpha) + Hd2;
        if (i0 < 0) i0 += 256;
        if ((j0 > H1+amp*surf1[i0*256+k]) && (j0 < H2-amp*surf2[i0*256+k]))
          aux->celltype = ZERO;
        else
          aux->celltype = PIERRE;

        //deuxieme fracture verticale
        static float* surf3 = NULL;
        static float* surf4 = NULL;
        if (!surf3) surf3 = load_surface("surf256c.dat", 256, 256);
        if (!surf4) surf4 = load_surface("surf256d.dat", 256, 256);
        amp = 1.0;

        int32_t L1 = (L/3)-5;
        int32_t L2 = (L/3)+5;

        if ((k > L1+amp*surf3[i+j*256]) && (k < L2-amp*surf4[i+j*256]))
          aux->celltype = ZERO;
*/
      }
#endif

/*****************************************************************************/
/********************************* ICB model *********************************/
/*****************************************************************************/
#ifdef MODEL_ICB
  for(k=0; k < D; k++) // profondeur
    for(j=0; j < H; j++) // hauteur
      for(i=0; i < L; i++, aux++){ //largeur
        if (j < H*2/3)
          aux->celltype   = (drand48() < 0.5)?ZERO:MOINS;
        else
          aux->celltype   = PLUS;
        //aux->celltype   = (int) floor(drand48()*3);
      }
#endif

/*****************************************************************************/
/********************************* RIV model *********************************/
/*****************************************************************************/
#ifdef MODEL_RIV
    float* surf;
    //int32_t lx = 501;
    //int32_t ly = 529;
    //int32_t ly1 = (L-ly)/2;
    //int32_t ly2 = ly1+ly;
    int32_t lx = 512; //128;
    int32_t ly = 512; //128;
    int32_t ly1 = 0;
    int32_t ly2 = ly-1; //ly*0.48;
    int32_t ll = Ld2;
    float niv = H*0.9;
    //float niv = H-1;
    float hout = H-1;
    float amp = 30.0;

    //modele precalcule (altimap)
    //surf = load_surface("data/guade.txt", lx, ly);
    //surf = load_surface("../../surface/surf128.dat", lx, ly);
    //surf = load_surface("surf512_05.txt", lx, ly);

    for(k=0; k < D; k++){ // profondeur
      for(j=0; j < H; j++){ // hauteur
        for(i=0; i < L; i++, aux++){ //largeur

          /*** EAU cells ***/
          aux->celltype = EAU; //default state

          /*** TERRE cells ***/
          if (csp_template.type == CSP_SLOPE){
            //CSP_SLOPE: sloping ground
            //format: SLOPE(r)
            float r = csp_template.args[0];
            if (j > niv-((D-1-k)*r)) aux->celltype = TERRE;
          }
          else{
            //CSP_CUSTOM
            //if ((k>=ly1) && (k<ly2) && (i<lx) && (j > H-surf[i+(k-ly1)*lx]/10.0))
            //if ((k>=ly1) && (k<ly2) && (i<lx) && (j > H-surf[i+(k-ly1)*lx]*amp))
            //if ((k>=ly1) && (k<ly2) && (i<lx) && (j > (H*0.4)-surf[i+(k-ly1+(L/2))*lx]*2.0))
            //if ((k>=ly1) && (k<ly2) && (i<lx) && (j > niv-((L-1-k)/2)-surf[i+(k-ly1)*lx]*2.0))
            //if ((k>=ly1) && (k<ly2) && (i<lx) && (j > niv-((L-1-k)/2)-surf[i+(k-ly1)*lx]/2.0)) //pente NS avec rugosite
            //if ((k>=ly1) && (k<ly2) && (i<lx) && (j > niv-((L-1-k)/2)-((i<Ld2)?(Ld2-i):(i-Ld2)))) //pente NS avec sillon
            if (j > niv)
            //if (j > niv-((L-1-k)*0.9)) //pente de 0.9
            //if (j > niv-((L-1-k)*0.8)+H*sin(5.0*i*6.28/L)/20) //pente de 0.8 + sinusoide
            //if ((j > niv-(k*0.9)) && (j > niv-((L-1-k)*0.9))) //triangle de pente 0.9
              aux->celltype = TERRE;
              //aux->celltype = (((int) floor(drand48()*100) == 0)) ? PIERRE : TERRE;
          }

          /*** DUM cells ***/
          //if ((i==0) || (k==0) || (i==L-1) || (k==D-1)) aux->celltype = DUM;  //barriere infranchissable tout autour (modele ferme)
            if ((i==0) || (k==0) || (i==L-1) || ((k==D-1) && (j>=hout))) aux->celltype = DUM;  //barriere infranchissable tout autour et jusqu'a une certaine hauteur sur le bord sud
            //if ((i==0) || (i==L-1) || (((k==0) || (k==D-1)) && (j>=hout))) aux->celltype = DUM;  //barriere infranchissable tout autour et jusqu'a une certaine hauteur sur les bords sud et nord
            //if ((i==0) || (i==L-1) || (((k==0) || (k==D-1)) && ((j>=hout) || (i!=Ld2)))) aux->celltype = DUM;  //barriere infranchissable tout autour et jusqu'a une certaine hauteur au centre des bords sud et nord
            //if ((i<Ld2-(ll/2)) || (i>Ld2+(ll/2))) aux->celltype = DUM; //couloir

          /*** IN cells ***/
           if ((j==1) && (i==Ld2) && (k==1)) aux->celltype = IN; //source
         }
       }
     }
#endif

#ifdef MODEL_RIV2
  LogPrintf("modele RIVIERE-2\n");
//#define WAVES
#ifdef WAVES
#define NB_WAV 6
  float offsetx[NB_WAV];
  float offsetz[NB_WAV];
  float alphaw[NB_WAV];
  float Hwave[NB_WAV];
  float Fwave[NB_WAV];
  int32_t iw;
  for (iw=0; iw<NB_WAV; iw++){
    offsetx[iw] = drand48()*10;
    offsetz[iw] = drand48()*10;
    alphaw[iw] = drand48();
    LogPrintf("alphaw[%d] = %f\n", iw, alphaw[iw]);
    //Hwave[iw] = (iw==0)? (float)H*0.2 : Hwave[iw-1]/1.2;
    //Fwave[iw] = (iw==0)? (float)6.28/L : Fwave[iw-1]*1.5;
    Hwave[iw] = (iw < 2)? (float)H*0.2 : Hwave[iw-2]/1.5;
    Fwave[iw] = (iw < 2)? (float)6.28/L : Fwave[iw-2]*2.0;
  }
#endif
  for(k=0; k < D; k++){ // profondeur
    for(i=0; i < L; i++){ //largeur
      //float alpha =  1.0 - (float)(i+k)/(L+L); //pente diagonale
      float alpha =  1.0 - (float)(k)/(D);  //pente nord-sud
      int32_t H1 = 1; //H*0.7; //H*0.1
      int32_t H2 = H*0.8; //H*0.4
      //int32_t H1 = Hd2*1.1;
      //int32_t H2 = Hd2*1.7;
      //float hh = H1*alpha + H2*(1-alpha); //sol en pente
      float hh = H2;//H*0.3;
      //float hh = (i<Ld2)? H*0.5 : H*0.9;//H*0.3; //falaise
      /*int32_t x1 = L/6;
      int32_t x2 = L*5/6;
      int32_t z1 = L/6;
      int32_t z2 = L*5/6;
      if ((i>=x1) && (i<x2) && (k>=z1) && (k<z2)) hh = H*0.1; //rectangle sureleve*/
      //float hh = H*0.2; //H*0.4;
      /*float dd = ((L/2) - i)*((L/2) - i) + ((L/2) - k)*((L/2) - k);
      if (dd < 3000) // bosse centrale;
        hh -= (3000-dd)/120;*/
      /*int32_t ii = i-k+Ld2; //fenetre pour le sillon
      int32_t L1 = L*0.1;
      int32_t L2 = L*0.5;
      if ((ii>L1) && (ii<L2)) //sillon diagonal
        hh += 3;*/
#ifdef WAVES
      hh = H*0.5;
      for (iw=1; iw<NB_WAV; iw++){
        hh += Hwave[iw] * cos(offsetx[iw]+(float)(alphaw[iw]*i+(1-alphaw[iw])*k)*Fwave[iw]) * cos(offsetz[iw]+(float)((1-alphaw[iw])*i-alphaw[iw]*k)*Fwave[iw]); //ondulations de surface
      }
#endif
      //hh += H*0.1*drand48();  //bruit blanc
      //int32_t L34 = L*0.75;
      //if ((i-L34)*(i-L34)+(k-L34)*(k-L34)<100) //encoche circulaire
      //  hh=H;
      //if ((i+k<L/2) || ((L-i+k)<L/2) || ((L-k+i)<L/2))  //carre en diagonal
      //  aux->celltype = DUM;
      //else if (i+k>3*L/2)
      //  aux->celltype = MOINS;
      //else
      float niv = H*0.9;
      for(j=0; j < H; j++){ // hauteur
        aux = TE + k*HL + j*L + i;
        //aux->celltype = (j<hh) ? MOINS : PLUS;  //(j<H-1) ? PLUS : DUM;
        aux->celltype = ((k<D*0.48) && (j > niv-((D-1-k)/2))) ? PLUS : DUM;
        //aux->celltype = (j<hh) ? MOINS : ((int) floor(drand48()*20) > 0) ? PLUS : PIERRE;
        if ((j==0) && (i==Ld2) && (k==1)) aux->celltype = IN;  //source ponctuelle
        //if (((i==0) || (k==0)) && (j>=niv)) aux->celltype = DUM;  //barriere infranchissable jusqu'a une certaine hauteur sur le bord
        if (((i==L-1) || (k==D-1)) && (j>=niv)) aux->celltype = DUM;  //barriere infranchissable jusqu'a une certaine hauteur sur le bord
        if ((i==0) || (k==0) || (i==L-1) || ((k==D-1) && (j>=niv))) aux->celltype = DUM;  //barriere infranchissable tout autour et jusqu'a une certaine hauteur sur le bord sud
        //if (((i==L-1) || (k==L-1)) && (i+k<L+L-10)) aux->celltype = DUM;  //petite ouverture en coin
        //if ((i==0) || (k==0) || (i==L-1) || ((k==L-1) && ((i<Ld2-2) || (i>Ld2+2)))) aux->celltype = DUM;  //petite ouverture au milieu
      }
    }
  }

#endif

#ifdef MODEL_DUM
// pour les tests ...
  for(k=0; k < D; k++) // profondeur
    for(j=0; j < H; j++) // hauteur
      for(i=0; i < L; i++, aux++){ //largeur
        aux->celltype = (i+k & 1) ? MOINS : PLUS;
      }
#endif

/*****************************************************************************/
/********************************* CRY model *********************************/
/*****************************************************************************/
#ifdef MODEL_CRY
  // sol plat et arrivee du gaz
  for(k=0; k < D; k++) // profondeur
    for(j=0; j < H; j++) // hauteur
      for(i=0; i < L; i++, aux++){ //largeur
        int32_t H1 = H*0.4;
        int32_t H2 = H*0.8;
        if ((j==H1) && (i<1) && (k<1))
          aux->celltype = IN;
        else
          aux->celltype = (j<H2) ? AIR : PLUS;
      }
#endif

/*****************************************************************************/
/********************************* DIF model *********************************/
/*****************************************************************************/
#ifdef MODEL_DIF
  for(k=0; k < D; k++) // profondeur
    for(j=0; j < H; j++) // hauteur
      for(i=0; i < L; i++, aux++){ //largeur
        if (i<Ld2)
          aux->celltype = ZERO;
        else
          aux->celltype = ONE;
      }
#endif

/*****************************************************************************/
/********************************* D2G model *********************************/
/*****************************************************************************/
#ifdef MODEL_D2G
  for(k=0; k < D; k++) // profondeur
    for(j=0; j < H; j++) // hauteur
      for(i=0; i < L; i++, aux++){ //largeur
        if ((j==0) && (i==Ld2) && (k==Dd2))
          aux->celltype = IN;
        else if (j>=H-5)
          aux->celltype = TERRE;
        else
          aux->celltype = AIR;
      }
#endif

/*****************************************************************************/
/********************************* LIFE model ********************************/
/*****************************************************************************/
#ifdef MODEL_LIFE
  for(k=0; k < D; k++) // profondeur
    for(j=0; j < H; j++) // hauteur
      for(i=0; i < L; i++, aux++){ //largeur
        //aux->celltype = (drand48()<0.5) ? ALIVE : DEAD; //random
        //aux->celltype = (k & 1) ? ALIVE : DEAD; //horizontal lines
        //if ((i-Ld2)*(i-Ld2) + (k-Dd2)*(k-Dd2) < Ld3*Ld3) //disk
        //if (((i-Ld2)*(i-Ld2) + (k-Dd2)*(k-Dd2) < Ld4*Ld4) || (k & 1)) //disk + horizontal lines
        if ((i==Ld2) || (k & 1)) //1 vertical line + horizontal lines
        //if ((drand48()<0.01) || (k & 1)) //random cells + horizontal lines
          aux->celltype = ALIVE;
        else
          aux->celltype = DEAD;
      }
#endif
}


float* load_surface(int8_t *nom, int32_t lx, int32_t ly)
{
  FILE *f;
  float *altimap, *pt_al;
  int32_t i, j;

  LogPrintf("lecture surface %s\n", nom);

  altimap = malloc(lx*ly*sizeof(float));

  //read ascii file
  f = fopen(nom,"r");
  if (!f){
    ErrPrintf("ERROR : can't open file %s\n", nom);
    exit(-1);
  }
  pt_al = altimap;
  for(j=0; j < ly; j++)
    for(i=0; i < lx; i++, pt_al++)
      fscanf(f, "%f \n", pt_al);

  fclose(f);

  return altimap;
}


void dump_terre(int32_t parx, int32_t pary)
{
  FILE *fp;
  int8_t str[256], *filename;
  Cell *pt;
  int32_t i, j;
  int8_t *buf, *aux, *ext;
  int32_t dump_type;


  if ((parx == 1) && (pary == 1)){
    // check extension
    if (csp_filename){
      filename = csp_filename;
      dump_type = DUMP_CSP;
      ext = filename + strlen(filename) - 4;
      if (strcmp(ext,".csp") && strcmp(ext,".CSP")){
        ErrPrintf("ERROR: invalid extension %s\n", filename);
        exit(-1);
      }
    }
    else if (bin_filename){
      filename = bin_filename;
      dump_type = DUMP_BIN;
      ext = filename + strlen(filename) - 4;
      if (strcmp(ext,".bin") && strcmp(ext,".BIN")){
        ErrPrintf("ERROR: invalid extension %s\n", filename);
        exit(-1);
      }
    }
    else{
      ext = opt_bin ? "bin" : "csp";
      dump_type = opt_bin ? DUMP_BIN : DUMP_CSP;
      sprintf(str,"%s_%d_%d_%d.%s",MOD_NAME,H,L,D,ext);
      filename = str;
    }
    LogPrintf("format %s\n", ext+1);
    csp_set_bounds(0,0,0);
    write_csp(dump_type, filename);
  }
  else{ //mode parallele (raw binary format) //TODO: csp format
    int32_t LX = L / parx;
    int32_t LY = D / pary;
    int32_t tx, ty, tn, k;
    LogPrintf("LX = %d\nLY = %d\n", LX, LY);
    buf = (char*)malloc(sizeof(char)*LX*H);
    for(tn=0, ty=0; ty<pary; ty++){
      for(tx=0; tx<parx; tx++, tn++){
        sprintf(str, "%s_%d_%d_%d_%d.bin", MOD_NAME, H, LX, LY, tn);
        filename = str;
        fp = fopen(filename,"w");
        if (!fp){
          ErrPrintf("erreur creation fichier %s\n", filename);
          exit(-4);
        }
        LogPrintf("creation fichier %s\n", filename);

        for(k=0; k<LY; k++){
          aux = buf;
          pt = TE + (ty*LY+k)*HL + tx*LX;
          for(j=0; j<H; j++, pt+=L){
            for(i=0; i<LX; i++){
              *aux++ = (pt+i)->celltype;
            }
          }
          fwrite(buf, sizeof(char), LX*H, fp);
        }
        fclose(fp);
        //LogPrintf("Data size: %d\n", H*LX*LY);
      }
    }
    free(buf);
  }
}

void dump_par_conf(int32_t parx, int32_t pary)
{
  FILE *fp;
  int32_t i, j, pid;
  int8_t nom[64];

  pid = 0;

  for (j=0; j<pary; j++)
    for (i=0; i<parx; i++, pid++){
      sprintf(nom, "parallel%d.cfg", pid);
      fp = fopen(nom, "w");
      if ( ! fp ){
        ErrPrintf("erreur ouverture fichier parallel.cfg\n");
        exit(-4);
      }
      LogPrintf("creation fichier de configuration %s\n", nom);
      fprintf(fp, "PID %d\n", pid);
      if (!pid) fprintf(fp, "NN %d\n", parx*pary);
      //if (parx > 1){
        if ((i<parx-1) || (boundary == BC_PERIODIC)) fprintf(fp, "EST %d\n", (i<parx-1)? pid+1 : pid+1-parx);
        if ((i>0) || (boundary == BC_PERIODIC)) fprintf(fp, "OUEST %d\n", (i>0)? pid-1 : pid-1+parx);
      //}
      //if (pary > 1){
        if ((j<pary-1) || (boundary == BC_PERIODIC)) fprintf(fp, "SUD %d\n", (j<pary-1)? pid+parx : i);
        if ((j>0) || (boundary == BC_PERIODIC)) fprintf(fp, "NORD %d\n", (j>0)? pid-parx : i+parx*(pary-1));
      //}
      fclose(fp);
    }

}


void params_genesis()
{
  param_family("GENESIS", "Genesis parameters");
  parameter("Csp_template", "CSP template for the generation of the initial cellular space", &csp_template_str, PARAM_STRING,"GENESIS");
#ifdef CYCLAGE_HOR
  parameter("Boundary", "boundary conditions (PERIODIC|OPEN|CLOSE|REINJECTION), PERIODIC by default", &boundary_str, PARAM_STRING,"GENERAL");
#else
  parameter("Boundary", "boundary conditions (OPEN|CLOSE|REINJECTION), REINJECTION by default", &boundary_str, PARAM_STRING,"GENERAL");
#endif
}

void usage()
{
  printf("genesis %s (%s)\n",VER_NUM,VER_DAT);
  printf("%s model\n",MOD_NAME);

  printf("\nusage: \n genesis H L D [OPTIONS]\n");
  printf(" genesis -f PARAMETERS_FILE [OPTIONS]\n\n");

  param_usage();

  printf("\nCSP TEMPLATES (%d)\n", nb_templates);
  for(int32_t i=0; i<nb_templates; i++) printf("  %s\n", t_templates[i].desc);

  printf("\nOPTIONS\n");
  printf("  -bin   \traw binary output\n");
  printf("  -s <n> \trandom seed\n");
#ifdef PARALLEL
  printf("  -par <n>x<m> \tparallel mode, grid <n>x<m>\n");
//  printf("  -node <n> \tnode number <n> in parallel mode\n");
#endif

  exit(-1);
}



int32_t main(int32_t argc, int8_t **argv)
{
  int32_t i, n=1, parx = 1, pary = 1;
  uint8_t opt_par = 0;

#ifdef LOG_FILE
  log_file = fopen("GENESIS.log","w");
#endif

  init_list_params();
  init_template_list();
  params_genesis();

  if (argc < 3) usage();

  LogPrintf("genesis\n");
  printf("%s model\n",MOD_NAME);

  //log the command line options
  LogPrintf("args: ");
  for(i=1; i<argc; i++) LogPrintf("%s ", argv[i]);
  LogPrintf("\n");

  if (!strcmp(argv[n],"-f")){
    int8_t *ext;
    // read parameters
    read_parameters(argv[++n]);
    //assert(filename && strlen(bin_filename)>=4);
    genesis_parse();
  }
  else{
    if (argc < 4) usage();
    H = atoi(argv[n]);
    L = atoi(argv[++n]);
    D = atoi(argv[++n]);
  }

  for(n++; n<argc; n++){
    if (!strcmp(argv[n],"-s") || !strcmp(argv[n],"-g")){
      graine = atoi(argv[++n]);
    }
    else if (!strcmp(argv[n],"-bin")){
      opt_bin = 1;
    }
#ifdef PARALLEL
    else if (!strcmp(argv[n],"-par")){
      sscanf(argv[++n], "%dx%d", &parx, &pary);
      LogPrintf("geom = %s parx = %d \t pary = %d\n", argv[n], parx, pary);
      opt_par = 1;
      L *= parx;
      D *= pary;
    }
/*
    else if (!strcmp(argv[n],"-node")){
      node = atoi(argv[++n]);
    }
*/
#endif
    else{
      ErrPrintf("ERROR: unknown option %s\n", argv[n]);
    }
  }

  LogPrintf("random seed: %d\n", graine);
  srand48(graine);

  if (opt_par) dump_par_conf(parx, pary);

  genesis();

  dump_terre(parx, pary);

  LogPrintf("genesis : done\n");

  return 1;
}



