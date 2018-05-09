/* ReSCAL - Surface processes
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


typedef struct pos_2d{
  short x;
  short y;
} Pos2;

typedef struct vector_2d{
  float x;
  float y;
} Vec2;

/*
typedef struct fifo_int{
  int *array;
  int length;
  int nb;
  int head;
  int tail;
} FifoInt;
*/

typedef struct fifo_pos_2d{
  Pos2 *array;
  int length;
  int nb;
  int head;
  int tail;
} FifoPos2;

enum ALTI_MODES {ALTI_MODE_BAS, ALTI_MODE_HAUT};

enum AVA_MODES {AVA_NONE, AVA_SYNC, AVA_TRANS, AVA_PROPAG};

#define Alti(_x,_y) alti[_x + (_y)*LEO]

void params_surface();
void modif_alti_cel(int ix, unsigned char typ);
void calcule_alti(unsigned char typ, char alti_mode);
int calcule_alti_max(unsigned char typ);
void calcule_normales();
int check_ava(int ix, void *data);
#ifdef LGCA
void init_cgv();
void reset_grad_vel();
void calcule_grad_vel();
float prob_cgv(float gv, float gvmin, float gvmax);
int check_grad_vel(int ix, void *data);
int check_grad_vel_color(int ix, void *data);
void compute_coef_cgv(/*int idb*/);
#endif
void avalanches(unsigned char typ, short h_lim, short nb_cel_max, char mode);
void avalanches_norm(unsigned char typ, short nb_cel_max, char alti_mode);
void ava_propag(int i, int j);

void dump_surface(char* name, int cpt, int unit);
void dump_sigma_alti();
void dump_autocorrel();
#ifdef LGCA
void dump_grad_vel(int cpt, int unit);
void dump_cgv();
void dump_cgv_coef();
#endif

