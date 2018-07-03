/* ReSCAL - Render and display module
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

#define IMG_W_MAX 5000

#define COUPE MAX_CELL

#ifdef NORM2D
#define VECNORM Vec2
#else
#define VECNORM Vec3
#endif

#define NB_SHADE_MAX 5

//enum IMG_FORMATS {IMG_PNG, IMG_JPEG};
enum VIEW_DIR {VDIR_NONE, VDIR_NORTH, VDIR_WIND};

typedef struct {
  int32_t col0;
  int32_t col1;
  int32_t nb_val;
  int32_t start;
} Shade;

void show_view_options();
int32_t view_init(int, char **);
void view_img_size(int32_t *, int32_t *);
void view_palette(int32_t *);
void rotate_light(float angle);
unsigned char* view();
void reset_zoom();
unsigned char* view_zoom(int32_t w, int32_t h, float coef);
void view_quit();
int32_t update_cv(int32_t x, int32_t y, int32_t flag);

#ifdef USE_LIBPNG
void view_dump_init();
void view_dump_png(char *nom);
#endif

#ifdef USE_GD
void view_dump_init();
//void view_dump_inter(int32_t inter, int32_t format);
void view_dump_png(char *nom);
void view_dump_jpeg(char *nom);
#endif

void dump_image(char *filename, char *format);
void dump_image_inter(int32_t inter, char *format);

