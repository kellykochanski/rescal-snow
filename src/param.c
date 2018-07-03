/* ReSCAL - Parsing of parameters
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


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>
#include <stdint.h>

#include "defs.h"
#include "macros.h"
#include "param.h"

#define LG_MAX 1024 /// maximum length of a line in parameters file

extern char *model_name;
extern char *bin_filename;
extern char *csp_filename;

extern int32_t H, L, D;
extern int32_t graine;

int32_t nb_families = 0;  /// number of families
Family *list_families = NULL; /// list of the families

int32_t nb_params = 0;  /// number of de parameters
Parameter *list_params = NULL; /// list of the parameters

/// DECLARATION OF A PARAMETER
/// par_name : name of the parameter
/// par_usage : definition of the parameter
/// par_adr : address of the parameter (in memory)
/// par_type : type of the parameter, possible values are PARAM_STRING | PARAM_INT | PARAM_FLOAT | PARAM_DOUBLE
/// par_family : family of the parameter
void parameter(char *par_nom, char *par_usage, void *par_adr, uint8_t par_type, char *par_family);


void init_list_families() {
  AllocMemory(list_families, Family, NB_PAR_FAM_MAX);
  ResetMemory(list_families, Family, NB_PAR_FAM_MAX);

  param_family("", "Other parameters");
}


void param_family(char *name, char *desc) {
  list_families[nb_families].name = strdup(name);
  list_families[nb_families].description = strdup(desc);
  list_families[nb_families].nb = 0;
  nb_families++;
  assert(nb_families < NB_PAR_FAM_MAX);
}

void display_family(uint8_t family) {
  int32_t j, l;

  assert(family < nb_families);

  if (list_families[family].nb) {
    printf("  %s:\n", list_families[family].description);
    //printf("  %s (%d):\n", list_families[family].description, list_families[family].nb);
    for (j = 0; j < nb_params; j++) {
      if ((list_params[j].family == family) && (list_params[j].is_visible)) {
        l = strlen(list_params[j].name);
        printf("    %s%s %s\n", list_params[j].name, (l <= 2) ? "\t\t" : "\t", list_params[j].usage);
      }
    }
  }
}

uint8_t get_family(char *name) {
  int8_t family = 0;

  if (name) {
    int32_t i = 1;
    while ((i < nb_families) && (strcmp(name, list_families[i].name))) {
      i++;
    }
    if (i < nb_families) {
      family = i;
    }
  }

  //LogPrintf("get_family: name=%s, family=%d\n", name, family);
  return family;
}

void init_list_params() {
  init_list_families();

  AllocMemory(list_params, Parameter, NB_PARAM_MAX);
  ResetMemory(list_params, Parameter, NB_PARAM_MAX);

  param_family("GENERAL", "General parameters");

  parameter("Model", "name of model", &model_name, PARAM_STRING, "GENERAL");

  /// dimensions of cellulat space
  parameter("H", "height", &H, PARAM_INT, "GENERAL");
  parameter("L", "length", &L, PARAM_INT, "GENERAL");
  parameter("D", "depth", &D, PARAM_INT, "GENERAL");

  /// input data file
  parameter("Bin_file", "initial cellular space (binary file)", &bin_filename, PARAM_STRING, "GENERAL");
  parameter("Csp_file", "initial cellular space (csp file)", &csp_filename, PARAM_STRING, "GENERAL");

  /// random seed
  parameter("Seed", "random seed", &graine, PARAM_INT, "GENERAL");
}


void param_usage() {
  int32_t i, j;

  printf("\nFORMAT OF PARAMETERS_FILE\n");
  //printf("\nFORMAT OF PARAMETERS_FILE (for \"%s\" model)\n", MOD_NAME);
  printf("  <parameter> = <value>\n\n");

  /// insert parameters within families
  for (j = 0; j < nb_params; j++) {
    i = get_family(list_params[j].family_name);
    list_params[j].family = i;
    if (list_params[j].is_visible) {
      list_families[i].nb++;
    }
  }

  /// display list of parameters sorted by family
  for (i = 1; i < nb_families; i++) {
    display_family(i);
  }
  display_family(0);
}

void bad_args() {
  ErrPrintf("ERROR: cannot read arguments in command-line\n");

  exit(-4);
}

void bad_params(char *str) {
  ErrPrintf("ERROR: cannot read parameter %s\n", str);

  exit(-4);
}


int32_t read_int(char *s, int32_t *err) {
  if (!isdigit(s[0]) && ((s[0] != '-') || !isdigit(s[1]))) {
    *err = 1;
    return 0;
  }
  return atoi(s);
}

double read_float(char *s, int32_t *err) {
  //LogPrintf("read_float : %s\n",s);
  if (!isdigit(s[0]) && ((s[0] != '-') || !isdigit(s[1]))) {
    *err = 1;
    return 0;
  }
  return atof(s);
}

int32_t read_boolean(char* str, int32_t *err) {
  if (!strcmp(str, "YES") || !strcmp(str, "1")) {
    return 1;
  } else if (!strcmp(str, "NO") || !strcmp(str, "0")) {
    return 0;
  } else {
    *err = 1;
    return 0;
  }
}

void parameter(char *par_name, char *par_usage, void *par_adr, uint8_t par_type, char *par_family) {
  list_params[nb_params].name = par_name;
  list_params[nb_params].usage = par_usage;
  list_params[nb_params].param = par_adr;
  list_params[nb_params].type = par_type;
  list_params[nb_params].is_visible = 1;
  //list_params[nb_params].family = get_family(par_family);
  list_params[nb_params].family_name = par_family;
  nb_params++;
  assert(nb_params < NB_PARAM_MAX);
}

void parameter_hidden(char *par_name, char *par_usage, void *par_adr, uint8_t par_type, char *par_family) {

  list_params[nb_params].name = par_name;
  list_params[nb_params].usage = par_usage;
  list_params[nb_params].param = par_adr;
  list_params[nb_params].type = par_type;
  list_params[nb_params].is_visible = 0;
  //list_params[nb_params].family = get_family(par_family);
  list_params[nb_params].family_name = par_family;
  nb_params++;
  assert(nb_params < NB_PARAM_MAX);
}

void init_val_param(char *par_name, char *par_val) {
  int32_t i = 0;
  int32_t err = 0;

  while ((i < nb_params) && (strcmp(list_params[i].name, par_name))) {
    i++;
  }

  if (i < nb_params) {
    LogPrintf("%s = %s\n", par_name, par_val);
    switch (list_params[i].type) {
    case PARAM_STRING:
      *(char**)list_params[i].param = (char*) malloc(strlen(par_val) + 1);
      strcpy(*(char**)list_params[i].param, par_val);
      break;
    case PARAM_BOOLEAN:
      *(int*)list_params[i].param = read_boolean(par_val, &err);
      break;
    case PARAM_INT:
      *(int*)list_params[i].param = read_int(par_val, &err);
      break;
    case PARAM_FLOAT:
      *(float*)list_params[i].param = (float) read_float(par_val, &err);
      break;
    case PARAM_DOUBLE:
      *(double*)list_params[i].param = read_float(par_val, &err);
      break;
    }
    if (err) {
      ErrPrintf("ERROR: bad value \"%s\" for parameter %s\n", par_val, par_name);
      exit(-1);
    }

    if (list_params[i].init) {
      ErrPrintf("ERROR: %s parameter is duplicated\n", par_name);
      exit(-1);
    }
    list_params[i].init = 1;
    if (!list_params[i].is_visible) {
      WarnPrintf("WARNING: %s\n", list_params[i].usage);
    }
  } else {
    LogPrintf("WARNING: unkown parameter - %s\n", par_name);
    //exit(-4);
  }
}

void read_line(FILE *fp, char *str) {
  char *ptr, *ptrmax, c;

  ptr = str;
  ptrmax = str + LG_MAX;

  /// read complete line with all spaces removed
  while (!feof(fp) && (ptr < ptrmax)) {
    c = fgetc(fp);
    if ((c == '\n') || (c == '\r') || (c == EOF)) {
      break;
    } else if ((c > 32) && (c < 127)) {
      *ptr++ = c;
    }
  }

  if (ptr >= ptrmax) {
    bad_params(str);
  }

  *ptr = 0;
}


void read_parameters(char *param_file) {
  FILE *fp;
  char *par_name;
  char *par_val;
  char str[LG_MAX];
  char *ptr;
  int32_t n = 0;

  LogPrintf("read parameters file : %s\n", param_file);

  fp = fopen(param_file, "r");
  if (!fp) {
    ErrPrintf("ERROR: cannot open file %s\n", param_file);
    exit(-4);
  }

  while (!feof(fp)) {
    *str = 0;
    read_line(fp, str);
    //LogPrintf("str=%s\n", str);
    if (*str == '#') { /// this is a comment
      //LogPrintf("commentaire ...\n");
    } else if (*str == 0) { /// empty line
      //LogPrintf("ligne vide ...\n");
    } else {
      ptr = strchr(str, '=');
      if (ptr) {
        /// we split the string str at '=' character
        *ptr = 0;
        par_name = str;
        par_val = ptr + 1;
        init_val_param(par_name, par_val);
        n++;
      } else {
        bad_params(str);
      }
    }
  }

  fclose(fp);

  LogPrintf("read_parameters : %d parametres lus\n", n);

  if (n < nb_params) {
    int32_t i = 0;
    for (i = 0; i < nb_params; i++)
      if (!list_params[i].init && list_params[i].is_visible) {
        LogPrintf("WARNING: parameter %s not found\n", list_params[i].name);
      }
  }
}

int32_t param_is_set(char *str) {
  int32_t i = 0;
  while ((i < nb_params) && (strcmp(list_params[i].name, str))) {
    i++;
  }

  return list_params[i].init;
}
