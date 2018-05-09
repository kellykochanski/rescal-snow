/* ReSCAL - CSP Format
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

#define CSP_MAGIC_NUM "\212CSP"

enum DUMP_TYPES {DUMP_CSP, DUMP_BIN, DUMP_PNG, DUMP_JPEG, DUMP_LOG};

void csp_set_warning(unsigned char check);
void csp_set_bounds(int xb, int yb, int zb);
void compress(char *filename, char force);
int read_csp_header(char *filename);
void read_csp(char *filename);
void write_csp(char dump_type, char *filename);
