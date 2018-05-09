/* ReSCAL - 3D space processes
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

typedef struct vector_3d{
  float x;
  float y;
  float z;
} Vec3;

enum ROT_MODES {ROT_MODE_DISK, ROT_MODE_OVERLAP, ROT_MODE_FULL};
//enum ROT_FLAGS {ROT_FLAG_NONE, ROT_FLAG_REORIENT_TEMP, ROT_FLAG_REORIENT_UNDO};
//rotation flags
//#define ROT_DISK 1
//#define ROT_CYCL 2
//#define ROT_CUBE 4
#define ROT_REORIENT_TEMP 1
#define ROT_REORIENT_UNDO 2

enum BOUNDARY_CONDITIONS {BC_PERIODIC, BC_OPEN, BC_OUT, BC_CLOSE, BC_REINJECTION};
enum SUR_MODES {SUR_MODE_UNIFORM, SUR_MODE_CONE, SUR_MODE_TECTO};

#define CellPtr(_i,_j,_k) (TE + _i + (_j)*L + (_k)*HL)
#define CellType(_i,_j,_k) TE[_i + (_j)*L + (_k)*HL].celltype
#define OutOfSpace(_x,_y) (rot_map[_x + L*(_y)] <= 0)
#define RotMapPos(_x,_y) rot_map_pos[_x + L*(_y)]
#define RotMapPosX(_x,_y) rot_map_pos[_x + L*(_y)].x
#define RotMapPosY(_x,_y) rot_map_pos[_x + L*(_y)].y

void init_terre();
//void wait_csp_ready();
void lock_csp(int log_flag);
void unlock_csp(int log_flag);
void cree_terre();
#ifdef PARALLEL
void cree_tunnels();
#endif
void apocalypse();
//int lecture_terre(char *nom);
void cellule_terre(unsigned char type, int ix);
int hcycle(int ix);
int get_cell_east(int ix);
int get_cell_west(int ix);
int get_cell_south(int ix);
int get_cell_north(int ix);
int get_cell_down(int ix);
int get_cell_up(int ix);
int get_cell_dir(int ix, char dir);
#if CELL_DATA
//void init_cell_data(CellData data, int ix);
//void get_cell_data(CellData data, int ix);
void swap_cell_data(int ix, int ix2);
void copy_cell_data(int ix, int ix2);
#endif
void init_color();
void init_lissage();
void lissage_terre();
void conditions_bords();
void calcule_couloir();
void compute_v_walls();
void compute_h_ceil();
void compute_h_floor();
Vec3 compute_mass_center(int type);

int check_alti(int ix, void *data);
int check_cell_dir(int ix, void *data);
int check_grain_seul(int index, void *data);

#ifdef SURRECTIONS
int surrection();
#endif
void translation(int dx, int dz);
void rotate_map(float angle);
void rotation(float angle, char mode, char flags);
void rotation90(int n);

void dump_terre(char dump_type, int cpt, int unit);
void dump_rugosi();
#ifdef DUMP_SIGNATURE
void dump_signature(int ii);
#endif
