/*ReSCAL - genesis entry
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


int32_t prog = PROG_TOOL;
int32_t  H, L, D, HL, HLD;        // size of cellular space
Cell *TE = NULL;        // cellular space
double csp_time = 0.0;
//int32_t  node = -1;            // node number (with PARALLEL option)
char *model_name = NULL; //name of model
char *bin_filename = NULL; //name of BIN file
char *csp_filename = NULL; //name of CSP file
char opt_bin = 0;   //raw binary format
int32_t graine = 68374; //random seed

char *csp_template_str = NULL; //value of CSP_template parameter
CSP_Template csp_template; //CSP template
CSP_Template t_templates[MAX_TEMPLATES]; //available templates
int32_t nb_templates = 0; //number of templates

char *boundary_str = NULL;
int32_t boundary = BC_PERIODIC;  //boundary conditions

double drand48();
void  srand48();

float* load_surface(char *, int, int);

//available CSP template types
#if defined(MODEL_DUN) || defined(MODEL_SNO)
enum CSP_TEMPLATES {INPUT_ELEVATION, CSP_CUSTOM, CSP_LAYER, CSP_SNOWFALL,  CSP_LAYER_COL, CSP_BLOCK, CSP_CYLINDER, CSP_CONE, CSP_RCONE, CSP_SNOWCONE, CSP_CONE2, CSP_CONE3, CSP_CONE5, CSP_RCONE5, CSP_RWALL, CSP_WAVES_2D, CSP_WAVY_NS_LAYER, CSP_WAVE, CSP_TRIANGLES, CSP_SRC_DISK, CSP_SRC_DISK_CEIL, CSP_SMILEY, CSP_FORSTEP};
#else
enum CSP_TEMPLATES {CSP_CUSTOM};
#endif

void init_template(int32_t type, char *name, char *desc, int32_t nb_args, ...) {
  va_list vl;
  int32_t i;

  t_templates[nb_templates].type = type;
  t_templates[nb_templates].name = name;
  t_templates[nb_templates].desc = desc;
  t_templates[nb_templates].nb_args = nb_args;
  t_templates[nb_templates].file = NULL;

  AllocMemory(t_templates[nb_templates].args, float, nb_args);
  va_start(vl, nb_args);
  for (i = 0; i < nb_args; i++) {
    if (t_templates[nb_templates].type == INPUT_ELEVATION && i == 0) {
      t_templates[nb_templates].file = (char*) va_arg(vl, const char*);
    } else { 
      t_templates[nb_templates].args[i] = (float) va_arg(vl, double);
      //LogPrintf("t_templates[%d].args[%d] = %f\n", nb_templates, i, t_templates[nb_templates].args[i]);      
    }
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

//initialization of available templates
void init_template_list() {
  init_template(CSP_CUSTOM, "CUSTOM", "CUSTOM:\t\t\tno template", 0);
#if defined(MODEL_DUN) || defined(MODEL_SNO)
  init_template(INPUT_ELEVATION, "INPUT_ELEVATION", "INPUT_ELEVATION(filename with elevation values, 0 or 1 if using injection layer)", 2, "ALTI00000_t0.data", 0);
  init_template(CSP_LAYER, "LAYER", "LAYER(h=1.0):\t\t\tsand layer of height <h>", 1, 1.0);
  init_template(CSP_SNOWFALL, "SNOWFALL", "SNOWFALL(h=1.0):\t\t\tsand layer of height <h> with IN cells along ceiling", 1, 1.0);
  init_template(CSP_SNOWCONE, "SNOWCONE", "SNOWCONE(h=30, w=100, x=L/2, y=D/2, d=h/5):\tcone of height <h> and width <w> centered on (x,y) with a layer depth <d>", 5, 30.0, 100.0, 0.0, 0.0, 0.0);
  init_template(CSP_BLOCK, "BLOCK", "BLOCK(xmin, xmax, ymin, ymax, h):\tsand block of height <h>", 5, 30.0, 50.0, 50.0, 150.0, 20.0);
  init_template(CSP_CYLINDER, "CYLINDER", "CYLINDER(h=20, w=100):\tcylinder of height <h> and width <w>", 2, 20.0, 100.0);
  //init_template(CSP_CONE, "CONE", "CONE(h=30, w=100):\t\tcone of height <h> and width <w>", 2, 30.0, 100.0);
  init_template(CSP_CONE, "CONE", "CONE(h=30, w=100, x=L/2, y=D/2):\t\tcone of height <h> and width <w> centered on (x,y)", 4, 30.0, 100.0, 0.0, 0.0);
  //init_template(CSP_RCONE, "RCONE", "RCONE(rh=0.66, rw=0.4):\tcone of relative height <rh>*H and relative width <rw>*Min(L,D)", 2, 0.66, 0.4);
  init_template(CSP_RCONE, "RCONE", "RCONE(rh=0.66, rw=0.4, rx=0.5, ry=0.5):\tcone of relative height <rh>*H, relative width <rw>*Min(L,D) and centered on relative position (<rx>*L,<ry>*D)", 4, 0.66, 0.4, 0.5, 0.5);
  init_template(CSP_CONE3, "CONE3", "CONE3(h=40, w=126):\t\t3 cones of height <h> and width <w>", 2, 40.0, 126.0);
  init_template(CSP_CONE5, "CONE5", "CONE5(h=40, w=126):\t\t5 cones of height <h> and width <w>", 2, 40.0, 126.0);
  init_template(CSP_RCONE5, "RCONE5", "RCONE5(rh=0.66, rw=0.2):\t5 cones of relative height <rh>*H and relative width <rw>*Min(L,D)", 2, 0.66, 0.2);
  init_template(CSP_RWALL, "RWALL", "RWALL(rh=0.5, rx=0.3):\twall of relative height <rh>*H and relative position <rx>*L", 2, 0.5, 0.3);
  init_template(CSP_WAVES_2D, "WAVES_2D", "WAVES_2D(per_min=10, per_max=60, ufd=0, amp=4, mh=30):\t2d sand waves with increasing periods from <per_min> to <per_max>, amplitude <amp> and mean height <mh> for stability analysis\n\t\t\t\tIt has uniform frequency distribution when <ufd> flag is set", 5, 10.0, 60.0, 0.0, 4.0, 30.0);
  init_template(CSP_WAVY_NS_LAYER, "WAVY_NS_LAYER", "WAVY_NS_LAYER(per=40, amp=15, h=5):\tnorth-south sand layer with wavy east boundary of period <per>, amplitude <amp> and height <h> for stability analysis", 3, 40.0, 15.0, 5.0);
  init_template(CSP_WAVE, "WAVE", "WAVE(height=15, ground=0):\tnorth-south wave of sand/snow with height h on top of hardened ground of thickness hg=0", 2, 15.0, 0.0);
  init_template(CSP_TRIANGLES, "TRIANGLES", "TRIANGLES(n=8, h=10, mh=20):\t<n> east-west periodic triangles of local height <h> and mean height <mh>", 3, 8.0, 10.0, 20.0);
  init_template(CSP_SRC_DISK, "SRC_DISK", "SRC_DISK(w=40, x=100, y=D/2):\tcircular source of sand on the ground with width <w> and centered on (x,y)", 3, 40.0, 100.0, 0.0);
  init_template(CSP_SRC_DISK_CEIL, "SRC_DISK_CEIL", "SRC_DISK_CEIL(w=40, x=100, y=D/2):\tcircular source of sand in the ceiling with width <w> and centered on (x,y)", 3, 40.0, 100.0, 0.0);
  init_template(CSP_SMILEY, "SMILEY", "SMILEY:\t\t\t:-)", 0);
  init_template(CSP_FORSTEP, "FORSTEP", "FORSTEP(h=0.15, hg=1, n=2, p=0):\tset of integer <n> solid forward-facing steps each <h> cells high, coated with a thickness <hg> of grains, perturbed in third dimension if p=1 (default p=0)", 4, 0.15, 0.1, 2.0, 0.0);
#endif
}

int32_t parse_template() {
  char *ptr;
  int32_t i;
  int32_t err = 0;

  //read template name
  ptr = strtok(csp_template_str, "(");
  LogPrintf("csp_template.name = %s\n", ptr);

  //find and copy the full template
  for (i = 0; (i < nb_templates) && strcmp(ptr, t_templates[i].name); i++);
  if (i < nb_templates) {
    memcpy(&csp_template, &t_templates[i], sizeof(CSP_Template));
    LogPrintf("csp_template.type = %d\n", csp_template.type);
    AllocMemory(csp_template.args, float, csp_template.nb_args);
    memcpy(csp_template.args, t_templates[i].args, csp_template.nb_args * sizeof(float));
  } else {
    ErrPrintf("ERROR: Incorrect value for Csp_template, %s\n", csp_template_str);
    exit(-2);
  }

  //read template optional arguments, if any
  ptr = strtok(NULL, ",)");
  i = 0;
  while (ptr && (i < csp_template.nb_args)) {
    //LogPrintf("csp_template.args[%d] = %s\n", i, ptr);
    if (csp_template.type != INPUT_ELEVATION || i > 0) {
      csp_template.args[i++] = read_float(ptr, &err);
    } else { 
      csp_template.file = ptr;
      LogPrintf("File name from input_elevation: %s\n", csp_template.file);
      i++;
    }

    if (err) {
      ErrPrintf("ERROR: bad argument \"%s\" for template %s\n", ptr, csp_template.name);
      exit(-1);
    }
    ptr = strtok(NULL, ",)");
  }
  if (ptr) {
    ErrPrintf("ERROR: Too many arguments in Csp_template (%s)\n", ptr);
    exit(-2);
  }
  for (i = 0; i < csp_template.nb_args; i++) {
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

  if (csp_template_str) {
    parse_template();
  }

  if (boundary_str) {
    if (!strcmp(boundary_str, "PERIODIC")) {
#ifdef CYCLAGE_HOR
      boundary = BC_PERIODIC;
#else
      ErrPrintf("Incorrect value for Boundary : %s, please define CYCLAGE_HOR in defs.h and recompile\n", boundary_str);
      exit(-2);
#endif
    } else if (!strcmp(boundary_str, "OPEN")) {
      boundary = BC_OPEN;
    } else if (!strcmp(boundary_str, "OUT")) {
      boundary = BC_OUT;
    } else if (!strcmp(boundary_str, "CLOSE")) {
      boundary = BC_CLOSE;
    } else if (!strcmp(boundary_str, "REINJECTION")) {
      boundary = BC_REINJECTION;
    } else {
      ErrPrintf("Incorrect value for Boundary : %s\n", boundary_str);
      exit(-2);
    }
    //pbc_mode = (boundary == BC_PERIODIC);
  }
}

void genesis() {
  int32_t i, j, k;
  Cell *aux;

  HL = H * L;
  HLD = HL * D;

  // allocation of the cellular space
  AllocMemoryPrint("TE", TE, Cell, HLD);
  ResetMemory(TE, Cell, HLD);

  aux = TE;


/*****************************************************************************/
/********************************* DUN model *********************************/
/*****************************************************************************/
#if defined(MODEL_DUN) || defined(MODEL_SNO)
  float di, dk;
  float Ldx, Ldy;
  float h, w, rc, n, hg, p;
  int32_t xmin, xmax, ymin, ymax;
  float x, y, z;

  int32_t j_value;
  int32_t j_found;
  int32_t injection_val = 0;

//// normal mode ////

  Ldx = ((int32_t) L / 2) - 0.5; //((int32_t) L / 2)-0.5; //((int32_t) L / 3);//Ld4;
  Ldy = /*0.66**/((int32_t) D / 2) - 0.5; //Ld4; //(int) L/2;
  int32_t hh = (int)(H / 6); //(100/3); //8+(100/3); //(H/3);//(60/2); //H/6; //H/2;
  int32_t H0 = H - 2; //H-8
  int32_t lc = D; //1; //400; //Ld4; //((int32_t) L / 2); //L/10; //Ld4; //largeur du couloir
  int32_t lcone = lc * 0.6; //lc*0.6; //lc/2; //lc*0.6; //lc*2/3;  //((int32_t) L / 2)*0.8; //largeur du cone
  int32_t hcone = ((int32_t) H / 3); //Hd2;//2*((int32_t) H / 3); //H-10; //hauteur du cone
  float periode = 10.0; //periode des ondulations

  FILE* file;		
  int32_t array_index = 0;		
  int32_t input_array[L*D];		

//// read integers from file and store in array ////		

  if (csp_template.type == INPUT_ELEVATION) {
    injection_val = csp_template.args[1];
    file = fopen(csp_template.file, "r");
    j_found = fscanf(file, "%d", &j_value);

    //store elevation values in input_array
    while (j_found == 1) {
      input_array[array_index] = j_value;
      j_found = fscanf(file, " %d", &j_value);
      array_index += 1;
    }
    fclose(file);
  }

  // initialization of the cellular space
  for (k = 0; k < D; k++) { // profondeur
    //periode = 10.0+(k-(((int32_t) D / 2)-lc/2))/2.0/4.0;
    periode = 10.0 + (k - (((int32_t) D / 2) - lc / 2)) / 2.0 / 2.0;
    //periode += 0.5*drand48();
    for (j = 0; j < H; j++) { // hauteur
      for (i = 0; i < L; i++, aux++) { //largeur
        di = i - Ldx;
        dk = k - Ldy; // + Ld4;
        rc = 0.5 * lcone * (j - (H - hcone)) / (hcone - 1);

        /*** EAUC cells ***/
        aux->celltype = EAUC; //default state

        /*** GR cells ***/

        // *** INPUT_ELEVATION CASE: no template (values read in from data file) *** //
        if (csp_template.type == INPUT_ELEVATION) {
            if (input_array[array_index] >= (H / 2)) {
              ErrPrintf("%d height in data file exceeds height (H) specified in param file. ", input_array[array_index]);
              exit(-2);
            } 
            //INPUT_ELEVATION CASE: elevation values read in from text file
            //format: INPUT_ELEVATION(text file name, 0 or 1 injection layer)
            if (injection_val == 1 && j == 1) {
              aux->celltype = IN;
            }
            array_index = i + L*k;
            if ((H - j) < input_array[array_index]) {
              aux->celltype = GR;
            }
	} else if (csp_template.type == CSP_LAYER) {
          //CSP_LAYER: sand layer
          //format: LAYER(h)
          hh = (int)csp_template.args[0];
          float frac = fmodf(csp_template.args[0], 1.0);
          if (j > H - 2 - hh) {
            aux->celltype = GR;
          }
          if ((j == H - 2 - hh) && frac) { //h is not an integer
            if (drand48() < frac) {
              aux->celltype = GR;
            }
          }
	} else if (csp_template.type == CSP_SNOWFALL) {        
	  //CSP_SNOWFALL: sand layer with injection grains at the top
	  //format: SNOWFALL(h)
	  if (j == 1) {
		aux->celltype = IN;
	  }
	  hh = (int)csp_template.args[0];
          float frac = fmodf(csp_template.args[0], 1.0);
          if (j > H - 2 - hh) {
		aux->celltype = GR;
          }
          if ((j == H - 2 - hh) && frac) { //h is not an integer
            if (drand48() < frac) {
		 aux->celltype = GR;
            }
          }	
  } else if (csp_template.type == CSP_SNOWCONE) {
          //CSP_SNOWCONE: snow layer seeded with a singular cone of a specified size
          //format: SNOWCONE(h, w, x, y, d)
          //parameters for the cone
          hcone = csp_template.args[0];
          lcone = csp_template.args[1];
          x = csp_template.args[2];
          y = csp_template.args[3];
          //parameters for the sand layer
          hh = (int)csp_template.args[4];
          float frac = fmodf(csp_template.args[4], 1.0);
          //adds injection grains as the top layer of the space, one layer thick
          if (j == 1) {
                aux->celltype = IN;
          }
          //adds the layer of standard grains along the bottom with the specified depth.
          if (j > H - 2 - hh) {
                aux->celltype = GR;
          }
          if ((j == H - 2 - hh) && frac) { //h is not an integer
            if (drand48() < frac) {
                 aux->celltype = GR;
            }
          }
          //builds the cone of grains
          if (x == 0) {
            x = Ldx;
          }
          if (y == 0) {
            y = Ldy;
          }
          rc = 0.5 * lcone * (j - (H - hcone)) / (hcone - 1);
          if ((j > (H - hcone)) && ((i - x) * (i - x) + (k - y) * (k - y) < rc * rc)) {
            aux->celltype = GR;
          }
	} else if (csp_template.type == CSP_BLOCK) {
          //CSP_BLOCK: sand block
          //format: BLOCK(xmin, xmax, ymin, ymax, h)
          xmin = (int)csp_template.args[0];
          xmax = (int)csp_template.args[1];
          ymin = (int)csp_template.args[2];
          ymax = (int)csp_template.args[3];
          hh = (int)csp_template.args[4];
          if ((i >= xmin) && (i < xmax) && (k >= ymin) && (k < ymax) && (j > H - 2 - hh)) {
            aux->celltype = GR;
            //aux->celltype = GRV;
          }
        } else if (csp_template.type == CSP_CYLINDER) {
          //CSP_CYLINDER: sand cylinder
          //format: CYLINDER(h, w)
          h = csp_template.args[0];
          w = csp_template.args[1];
          if ((j >= (H - h - 1)) && (di * di + dk * dk < 0.25 * w * w)) {
            aux->celltype = GR;
          }
        } else if (csp_template.type == CSP_CONE) {
          //CSP_CONE: sand cone
          //format: CONE(h, w, x, y)
          hcone = csp_template.args[0];
          lcone = csp_template.args[1];
          x = csp_template.args[2];
          y = csp_template.args[3];
          if (x == 0) {
            x = Ldx;
          }
          if (y == 0) {
            y = Ldy;
          }
          rc = 0.5 * lcone * (j - (H - hcone)) / (hcone - 1);
          if ((j > (H - hcone)) && ((i - x) * (i - x) + (k - y) * (k - y) < rc * rc)) {
            aux->celltype = GR;
          }
        } else if (csp_template.type == CSP_RCONE) {
          //CSP_RCONE: sand cone with relative dimensions and position
          //format: RCONE(rh, rw, rx, ry)
          hcone = csp_template.args[0] * H;
          lcone = csp_template.args[1] * Min(L, D);
          x = csp_template.args[2] * L - 0.5;
          y = csp_template.args[3] * D - 0.5;
          rc = 0.5 * lcone * (j - (H - hcone)) / (hcone - 1);
          if ((j > (H - hcone)) && ((i - x) * (i - x) + (k - y) * (k - y) < rc * rc)) {
            aux->celltype = GR;
          }
          //if (j > H-3) aux->celltype = GR; //uniform layer
        } else if (csp_template.type == CSP_CONE2) {
          //CSP_CONE2: 2 sand cones for dune collisions, with different colors
          //format: CONE2
          hcone = H / 2;
          float lcone0 = D / 2;
          float lcone1 = D / 3;
          float rc0 = 0.5 * lcone0 * (j - (H - hcone)) / (hcone - 1);
          float rc1 = 0.5 * lcone1 * (j - (H - hcone)) / (hcone - 1);
          float di0 = i - ((int32_t) L / 3);
          float dk0 = k - D * 3 / 5;
          float di1 = i - L / 8;
          float dk1 = k - D / 3;
          if ((j > (H - hcone)) && (di0 * di0 + dk0 * dk0 < rc0 * rc0)) {
            aux->celltype = GR;
          }
          if ((j > (H - hcone)) && (di1 * di1 + dk1 * dk1 < rc1 * rc1)) {
            aux->celltype = GR;
          }
        } else if (csp_template.type == CSP_CONE3) {
          //CSP_CONE3: 3 sand cones with fixed dimensions
          //format: CONE3(h, w)
          hcone = csp_template.args[0];
          lcone = csp_template.args[1];
          rc = 0.5 * lcone * (j - (H - hcone)) / (hcone - 1);
          float di0 = i - Ldx - 0.5 * Ldx * cos(0);
          float dk0 = k - Ldy - 0.5 * Ldy * sin(0);
          float di2 = i - Ldx - 0.5 * Ldx * cos(2 * PI / 3);
          float dk2 = k - Ldy - 0.5 * Ldy * sin(2 * PI / 3);
          float di4 = i - Ldx - 0.5 * Ldx * cos(4 * PI / 3);
          float dk4 = k - Ldy - 0.5 * Ldy * sin(4 * PI / 3);
          if ((j > (H - hcone)) && (di0 * di0 + dk0 * dk0 < rc * rc)) {
            aux->celltype = GR;
          }
          if ((j > (H - hcone)) && (di2 * di2 + dk2 * dk2 < rc * rc)) {
            aux->celltype = GR;
          }
          if ((j > (H - hcone)) && (di4 * di4 + dk4 * dk4 < rc * rc)) {
            aux->celltype = GR;
          }
        } else if ((csp_template.type == CSP_CONE5) || (csp_template.type == CSP_RCONE5)) {
          if (csp_template.type == CSP_CONE5) {
            //CSP_CONE5: 5 sand cones with fixed dimensions
            //format: CONE5(h, w)
            hcone = csp_template.args[0];
            lcone = csp_template.args[1];
          } else {
            //CSP_RCONE5: 5 sand cones with relative dimensions
            //format: RCONE5(rh, rw)
            hcone = csp_template.args[0] * H;
            lcone = csp_template.args[1] * Min(L, D);
          }
          rc = 0.5 * lcone * (j - (H - hcone)) / (hcone - 1);
          float di1 = i - Ldx - 0.6 * Ldx * cos(PI / 5);
          float dk1 = k - Ldy - 0.6 * Ldy * sin(PI / 5);
          float di3 = i - Ldx - 0.6 * Ldx * cos(3 * PI / 5);
          float dk3 = k - Ldy - 0.6 * Ldy * sin(3 * PI / 5);
          float di5 = i - Ldx - 0.6 * Ldx * cos(5 * PI / 5);
          float dk5 = k - Ldy - 0.6 * Ldy * sin(5 * PI / 5);
          float di7 = i - Ldx - 0.6 * Ldx * cos(7 * PI / 5);
          float dk7 = k - Ldy - 0.6 * Ldy * sin(7 * PI / 5);
          float di9 = i - Ldx - 0.6 * Ldx * cos(9 * PI / 5);
          float dk9 = k - Ldy - 0.6 * Ldy * sin(9 * PI / 5);
          if ((j > (H - hcone)) && (di1 * di1 + dk1 * dk1 < rc * rc)) {
            aux->celltype = GR;
          }
          if ((j > (H - hcone)) && (di3 * di3 + dk3 * dk3 < rc * rc)) {
            aux->celltype = GR;
          }
          if ((j > (H - hcone)) && (di5 * di5 + dk5 * dk5 < rc * rc)) {
            aux->celltype = GR;
          }
          if ((j > (H - hcone)) && (di7 * di7 + dk7 * dk7 < rc * rc)) {
            aux->celltype = GR;
          }
          if ((j > (H - hcone)) && (di9 * di9 + dk9 * dk9 < rc * rc)) {
            aux->celltype = GR;
          }
          //if ((j > (H-hcone)) && (di*di + dk*dk < rc*rc)) aux->celltype = GR; //central cone
          //if (j > H-3) aux->celltype = GR; //uniform layer
        } else if (csp_template.type == CSP_RWALL) {
          //CSP_RWALL(rh,rx): north-south wall or DUM cells
          int32_t h = csp_template.args[0] * H;
          int32_t x = csp_template.args[1] * L;
          if ((i >= x - 1) && (i <= x + 1) && (j >= H - h)) {
            aux->celltype = DUM;
          }
        } else if (csp_template.type == CSP_WAVES_2D) {
          //CSP_WAVES_2D: 2d sand waves with increasing periods for stability analysis
          float per_min = csp_template.args[0];
          float per_max = csp_template.args[1];
          float ufd_flag = csp_template.args[2];
          float amp = csp_template.args[3];
          int32_t mh = (int)csp_template.args[4];
          if (per_min <= 0) {
            ErrPrintf("ERROR: bad parameter in WAVES_2D template, <per_min> must be positive (%f)\n", per_min);
            exit(-1);
          }
          if (ufd_flag) {
            // uniform frequency distribution
            float fmin = 1.0 / per_max;
            float fmax = 1.0 / per_min;
            float f = fmax + (fmin - fmax) * k / D;
            periode = 1.0 / f;
          } else {
            periode = per_min + (per_max - per_min) * k / D;
          }
          if ((i == 0) && (j == 0) && (k % 2 == 0)) {
            LogPrintf("k = %d, periode = %f\n", k, periode);
          }
          if (j > (H - mh) - 0.5 * amp * sin(i * 2 * PI / periode)) {
            aux->celltype = GR;  //variable ondulations
          }
          if ((k & 1) || (j > H0) || (j == 0)) {
            aux->celltype = DUM;  //ground + ceiling + walls east-west
          }
        } else if (csp_template.type == CSP_WAVE) {
            // north-south triangular wave
	    float h = csp_template.args[0];
	    float hg = csp_template.args[1];
	    float i2 = i - L/4.0;
	    // Hardened layer below wave
	    if ( j >= H - hg ) aux -> celltype = GRV;
	    // Downwind half of triangular wave
	    else if ( ( i2 >= 0 ) && ( i2 < 2.0*h )
			    && (j >= H - hg - h*(1.0-2*i2/h)) ) aux -> celltype = GR;
	    // Upwind half of triangular wave
	    else if ( ( i2 <  0 ) && ( i2 >= -2.0*h )
			    && (j >= H - hg - h*(1.0+2*i2/h)) ) aux -> celltype = GR;

        } else if (csp_template.type == CSP_WAVY_NS_LAYER) {
          //CSP_WAVY_NS_LAYER: north-south sand layer with wavy east boundary for stability analysis
          float per = csp_template.args[0];
          float amp = csp_template.args[1];
          int32_t hh = (int)csp_template.args[2];
          int32_t i_wave = (int)(40 + 2 + 0.5 * amp * (1 + sin(k * 2 * PI / per)));
          if ((i >= 40) && (i <= i_wave) && (j <= H0) && (j > H0 - hh)) {
            aux->celltype = GR;  //rangee de sable avec front de forme sinusoidale
          }
        } else if (csp_template.type == CSP_TRIANGLES) {
          //CSP_TRIANGLES: east-west periodic triangles
          float lx0 = 0.5 * L / csp_template.args[0];
          float lx1 = fmodf(i, lx0 * 2);
          hh = csp_template.args[1];
          int32_t mh = csp_template.args[2];
          if (((lx1 < lx0) && (j > H0 - mh + hh / 2 - (hh * lx1 / lx0))) || ((lx1 >= lx0) && (j > H0 - mh + hh / 2 - (hh * (2 * lx0 - lx1) / lx0)))) {
            aux->celltype = GR;
          }
        } else if ((csp_template.type == CSP_SRC_DISK) || (csp_template.type == CSP_SRC_DISK_CEIL)) {
          //CSP_SRC_DISK: circular source of sand
          w = csp_template.args[0];
          x = csp_template.args[1];
          y = csp_template.args[2];
          if (y == 0) {
            y = D * 0.5;
          }
          z = (csp_template.type == CSP_SRC_DISK) ? H0 + 1 : 0 ;
          if ((j == z) && ((i - x) * (i - x) + (k - y) * (k - y) <= w * w / 4.0)) {
            aux->celltype = IN;
          }
        } else if (csp_template.type == CSP_SMILEY) {
          //CSP_SMILEY: :-)
          if (j<H/2 && ((cos(PI*(i-L/2)/L/2)*cos(PI*(k-D/2)/D/2)<0.8 && cos(PI*(i-L/2)/L/2)*cos(PI*(k-D/2)/D/2)>0.75)
             || (cos(2*PI*(i-L/2)/L)<0.4 && sin(PI/4+2*PI*(k-D/2)/D)<0.4 && cos(2*PI*(i-L/2)/L)>0.2 && sin(PI/4+2*PI*(k-D/2)/D)>0.2 && k<L/2)
             || (cos(PI*(i-L/2)/L)*cos(PI*(k-D/2)/D)<0.8 && cos(PI*(i-L/2)/L)*cos(PI*(k-D/2)/D)>0.75&&k>3*D/5)))
               aux->celltype = GR;
        }
        else if (csp_template.type == CSP_FORSTEP){
	  //CSP_FORSTEP: set of forward-facing steps transverse to wind
          h  = csp_template.args[0]; // height of each step
	  n  = csp_template.args[1]; // number of steps
	  hg = csp_template.args[2]; // thickness of space that is granular
	  p  = csp_template.args[3]; // binary - is step perturbed in 3rd dimension?
	  float b  = 2.0*h;      // thickness of solid layer below steps
	  float p0 = p*L/n/20.0; // amplitude of perturbation, if any (0 if p=0)
	  float l  = D/2.0;      // wavelength of perturbation, if any
	  for( float step = -1.0; step < n; step = step + 1.0 ){
            float m_s = n*h/L;
	    float i2  = i + p0*sin(2*M_PI*k/l); //apply perturbation
	    if ( (i2>=L*step/n) 
			    && (i2<L*(step+1.0)/n) 
			    && (j>=(H-b-(h-m_s*(i2-L*step/n))))
	       ) 
#ifdef GRV  // if possible, make cohesive (GRV) steps
		    aux->celltype = GRV;
#else       // else make unerodible (DUM) steps
		    aux->celltype = DUM; 
#endif
	    if ((hg>=1) && (j<= hg)) aux->celltype = GR;
	    if ((hg<1)  && (j==1) && (k%(int) round(1/hg)==0)) aux->celltype = GR;
	  }
	}
        else{
          //CSP_CUSTOM
          if (j > H-10) //couche uniforme
          {
            aux->celltype = GR;
          }
        }

        /*** DUM cells ***/
        if ((j > H0) || (j == 0)) //plateau + plafond
          if (aux->celltype != IN) {
            aux->celltype   = DUM;
          }

        /*** IN cells ***/
        if (boundary == BC_REINJECTION) {
          if ((i == 5) && (j == 0)) { //source de grains en haut a gauche
            aux->celltype = IN;
          }
        }
      }
    }
  }
#endif
}


float* load_surface(char *nom, int32_t lx, int32_t ly) {
  FILE *f;
  float *altimap, *pt_al;
  int32_t i, j;

  LogPrintf("lecture surface %s\n", nom);

  altimap = malloc(lx * ly * sizeof(float));

  //read ascii file
  f = fopen(nom, "r");
  if (!f) {
    ErrPrintf("ERROR : can't open file %s\n", nom);
    exit(-1);
  }
  pt_al = altimap;
  for (j = 0; j < ly; j++)
    for (i = 0; i < lx; i++, pt_al++) {
      fscanf(f, "%f \n", pt_al);
    }

  fclose(f);

  return altimap;
}


void dump_terre(int32_t parx, int32_t pary) {
  FILE *fp;
  char str[256], *filename;
  Cell *pt;
  int32_t i, j;
  char *buf, *aux, *ext;
  int32_t dump_type;


  if ((parx == 1) && (pary == 1)) {
    // check extension
    if (csp_filename) {
      filename = csp_filename;
      dump_type = DUMP_CSP;
      ext = filename + strlen(filename) - 4;
      if (strcmp(ext, ".csp") && strcmp(ext, ".CSP")) {
        ErrPrintf("ERROR: invalid extension %s\n", filename);
        exit(-1);
      }
    } else if (bin_filename) {
      filename = bin_filename;
      dump_type = DUMP_BIN;
      ext = filename + strlen(filename) - 4;
      if (strcmp(ext, ".bin") && strcmp(ext, ".BIN")) {
        ErrPrintf("ERROR: invalid extension %s\n", filename);
        exit(-1);
      }
    } else {
      ext = opt_bin ? "bin" : "csp";
      dump_type = opt_bin ? DUMP_BIN : DUMP_CSP;
      sprintf(str, "%s_%d_%d_%d.%s", MOD_NAME, H, L, D, ext);
      filename = str;
    }
    LogPrintf("format %s\n", ext + 1);
    csp_set_bounds(0, 0, 0);
    write_csp(dump_type, filename);
  } else { //mode parallele (raw binary format) //TODO: csp format
    int32_t LX = L / parx;
    int32_t LY = D / pary;
    int32_t tx, ty, tn, k;
    LogPrintf("LX = %d\nLY = %d\n", LX, LY);
    buf = (char*)malloc(sizeof(char) * LX * H);
    for (tn = 0, ty = 0; ty < pary; ty++) {
      for (tx = 0; tx < parx; tx++, tn++) {
        sprintf(str, "%s_%d_%d_%d_%d.bin", MOD_NAME, H, LX, LY, tn);
        filename = str;
        fp = fopen(filename, "w");
        if (!fp) {
          ErrPrintf("erreur creation fichier %s\n", filename);
          exit(-4);
        }
        LogPrintf("creation fichier %s\n", filename);

        for (k = 0; k < LY; k++) {
          aux = buf;
          pt = TE + (ty * LY + k) * HL + tx * LX;
          for (j = 0; j < H; j++, pt += L) {
            for (i = 0; i < LX; i++) {
              *aux++ = (pt + i)->celltype;
            }
          }
          fwrite(buf, sizeof(char), LX * H, fp);
        }
        fclose(fp);
        //LogPrintf("Data size: %d\n", H*LX*LY);
      }
    }
    free(buf);
  }
}

void dump_par_conf(int32_t parx, int32_t pary) {
  FILE *fp;
  int32_t i, j, pid;
  char nom[64];

  pid = 0;

  for (j = 0; j < pary; j++)
    for (i = 0; i < parx; i++, pid++) {
      sprintf(nom, "parallel%d.cfg", pid);
      fp = fopen(nom, "w");
      if (! fp) {
        ErrPrintf("erreur ouverture fichier parallel.cfg\n");
        exit(-4);
      }
      LogPrintf("creation fichier de configuration %s\n", nom);
      fprintf(fp, "PID %d\n", pid);
      if (!pid) {
        fprintf(fp, "NN %d\n", parx * pary);
      }
      //if (parx > 1){
      if ((i < parx - 1) || (boundary == BC_PERIODIC)) {
        fprintf(fp, "EST %d\n", (i < parx - 1) ? pid + 1 : pid + 1 - parx);
      }
      if ((i > 0) || (boundary == BC_PERIODIC)) {
        fprintf(fp, "OUEST %d\n", (i > 0) ? pid - 1 : pid - 1 + parx);
      }
      //}
      //if (pary > 1){
      if ((j < pary - 1) || (boundary == BC_PERIODIC)) {
        fprintf(fp, "SUD %d\n", (j < pary - 1) ? pid + parx : i);
      }
      if ((j > 0) || (boundary == BC_PERIODIC)) {
        fprintf(fp, "NORD %d\n", (j > 0) ? pid - parx : i + parx * (pary - 1));
      }
      //}
      fclose(fp);
    }

}


void params_genesis() {
  param_family("GENESIS", "Genesis parameters");
  parameter("Csp_template", "CSP template for the generation of the initial cellular space", &csp_template_str, PARAM_STRING, "GENESIS");
#ifdef CYCLAGE_HOR
  parameter("Boundary", "boundary conditions (PERIODIC|OPEN|CLOSE|REINJECTION), PERIODIC by default", &boundary_str, PARAM_STRING, "GENERAL");
#else
  parameter("Boundary", "boundary conditions (OPEN|CLOSE|REINJECTION), REINJECTION by default", &boundary_str, PARAM_STRING, "GENERAL");
#endif
}

void usage() {
  printf("genesis %s (%s)\n", VER_NUM, VER_DAT);
  printf("%s model\n", MOD_NAME);

  printf("\nusage: \n genesis H L D [OPTIONS]\n");
  printf(" genesis -f PARAMETERS_FILE [OPTIONS]\n\n");

  param_usage();

  printf("\nCSP TEMPLATES (%d)\n", nb_templates);
  for (int32_t i = 0; i < nb_templates; i++) {
    printf("  %s\n", t_templates[i].desc);
  }

  printf("\nOPTIONS\n");
  printf("  -bin   \traw binary output\n");
  printf("  -s <n> \trandom seed\n");

  exit(-1);
}



int32_t main(int32_t argc, char **argv) {
  int32_t i, n = 1, parx = 1, pary = 1;
  uint8_t opt_par = 0;

#ifdef LOG_FILE
  char *filename = output_path("GENESIS");
  log_file = fopen(filename,"w");
  free(filename);
#endif

  init_list_params();
  init_template_list();
  params_genesis();

  if (argc < 3) {
    usage();
  }

  LogPrintf("genesis\n");
  printf("%s model\n", MOD_NAME);

  //log the command line options
  LogPrintf("args: ");
  for (i = 1; i < argc; i++) {
    LogPrintf("%s ", argv[i]);
  }
  LogPrintf("\n");

  if (!strcmp(argv[n], "-f")) {
//     char *ext;
    // read parameters
    read_parameters(argv[++n]);
    //assert(filename && strlen(bin_filename)>=4);
    genesis_parse();
  } else {
    if (argc < 4) {
      usage();
    }
    H = atoi(argv[n]);
    L = atoi(argv[++n]);
    D = atoi(argv[++n]);
  }

  for (n++; n < argc; n++) {
    if (!strcmp(argv[n], "-s") || !strcmp(argv[n], "-g")) {
      graine = atoi(argv[++n]);
    } else if (!strcmp(argv[n], "-bin")) {
      opt_bin = 1;
    } else {
      ErrPrintf("ERROR: unknown option %s\n", argv[n]);
    }
  }

  LogPrintf("random seed: %d\n", graine);
  srand48(graine);

  if (opt_par) {
    dump_par_conf(parx, pary);
  }

  genesis();

  dump_terre(parx, pary);

  LogPrintf("genesis : done\n");

  return 1;
}



