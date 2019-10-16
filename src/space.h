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
 * aint64_t with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

#include <stdint.h>
#include <sys/time.h>

typedef struct vector_3d {
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

double calculate_msec(struct timeval start, struct timeval stop);
void log_time_delta(struct timeval start, char* message, int id, int visible);
void init_terre();
//void wait_csp_ready();
void lock_csp(int32_t log_flag);
void unlock_csp(int32_t log_flag);
void cree_terre();
void apocalypse();
//int32_t lecture_terre(char *nom);
void cellule_terre(uint8_t type, int32_t ix);
int32_t hcycle(int32_t ix);
int32_t get_cell_east(int32_t ix);
int32_t get_cell_west(int32_t ix);
int32_t get_cell_south(int32_t ix);
int32_t get_cell_north(int32_t ix);
int32_t get_cell_down(int32_t ix);
int32_t get_cell_up(int32_t ix);
int32_t get_cell_dir(int32_t ix, char dir);
#ifdef CELL_DATA
void swap_cell_data(int32_t ix, int32_t ix2);
void copy_cell_data(int32_t ix, int32_t ix2);
#endif
void init_color();
void init_lissage();
void lissage_terre();
void conditions_bords();
void calcule_couloir();
void compute_v_walls();
void compute_h_ceil();
void compute_h_floor();
Vec3 compute_mass_center(int32_t type);

int32_t check_alti(int32_t ix, void *data);
int32_t check_cell_dir(int32_t ix, void *data);
int32_t check_grain_seul(int32_t index, void *data);

void translation(int32_t dx, int32_t dz);
void rotate_map(float angle);
void rotation(float angle, char mode, char flags);
void rotation90(int32_t n);

void dump_terre(char dump_type, int32_t cpt, int32_t unit);
void dump_rugosi();
#ifdef DUMP_SIGNATURE
void dump_signature(int32_t ii);
#endif
