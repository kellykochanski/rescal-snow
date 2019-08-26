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

#include <stdint.h>

typedef struct _template {
  int32_t type;
  char *name;
  char *desc;
  int32_t nb_args;
  float *args;
  char *file;
} CSP_Template;

#define MAX_TEMPLATES 20 //max number of templates

enum BOUNDARY_CONDITIONS {BC_PERIODIC, BC_OPEN, BC_OUT, BC_CLOSE, BC_REINJECTION};

