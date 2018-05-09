/* ReSCAL - Parsing of parameters
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
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

#define NB_PAR_FAM_MAX 100 /// max number of (parameter) families

#define NB_PARAM_MAX 200 /// max number of parameters


typedef struct{
  char *name;
  char *description;
  int nb;
} Family;

typedef struct{
  char *name;
  char *usage;
  void *param;
  unsigned char type;
  unsigned char init;
  unsigned char is_visible;
  unsigned char family;
  char *family_name;
} Parameter;

enum PARAM_TYPES {PARAM_STRING, PARAM_BOOLEAN, PARAM_INT, PARAM_FLOAT, PARAM_DOUBLE};

void init_list_params();
void param_usage();
void param_family(char *name, char *desc);
void parameter(char *par_nom, char *par_usage, void *par_adr, unsigned char par_type, char *par_family);
void parameter_hidden(char *par_nom, char *par_usage, void *par_adr, unsigned char par_type, char *par_family);
int read_int(char *s, int *err);
double read_float(char *s, int *err);
void read_parameters(char *param_file);
int param_is_set(char *str);

