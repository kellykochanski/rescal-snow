/* ReSCAL - Cells tracing
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



#define MODE_TRACE defined(TRACE_SRC) || defined(TRACE_AIRE) || defined(TRACE_PAR) || defined(TRACE_PAR_COL) || defined(TRACE_TRANS) || defined(TRACE3D_CEL) || defined(TRACE_FLUX)

void trace_init();
int trace_init_loop();
void trace_end_loop();

#if defined(TRACE_SRC) || defined(TRACE_AIRE)
void trace_test(int ix);
int trace_point(int x, int z);
#endif

#ifdef TRACE_TRANS
//int trace_trans_init();
void trace_trans(int tr, int ix);
void trace_trans_blocked(int tr, int ix);
void trace_trans_dump();
void trace_trans_quit();
#endif

#ifdef TRACE_PAR
void trace_par(int ix);
//void trace_par_fin();
#endif

#ifdef TRACE_PAR_COL
void trace_par_col_init();
void trace_par_col_rotation(float angle);
void trace_par_col(int ix, int tr);
#endif

#ifdef TRACE3D_CEL
void trace3d_cel();
void trace3d_cel_dump();
void trace3d_cel_quit();
#endif

#ifdef TRACE_FLUX
typedef struct flux_map{
  float angle;
  float *flux_ew;
  float *flux_ns;
} FluxMap;

void trace_flux_init();
void trace_flux(int tr, int ix);
void trace_flux_rotate(float angle);
void trace_flux_dump();
#endif

void params_trace();
void trace_dump(char flag_info);

