/* ReSCAL - GTK callbacks
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

//#define FLASH_DELAY 50


typedef struct dump_parameters {
  uint8_t type;
  int32_t delay;
} DumpPar;

typedef struct dump_delay {
  uint8_t unit;
  float val;
} DumpDelay;

void callbacks_init();
void lock_display(int32_t log_flag);
int32_t trylock_display(int32_t log_flag);
void unlock_display(int32_t log_flag);
int32_t elapsed(double *sec);
void timer_init();
void *rescal_thread(void *data);
void *lgca_thread(void *data);
void wait_lgca_thread();
void do_thread_sched();
void log_info();

#ifdef __GTK_H__
gboolean gcallback_dump(gpointer  data);
gboolean gcallback_stop(gpointer  data);
gboolean gcallback_quit(gpointer data);
gboolean gcallback_flash(gpointer  data);

void on_drawingarea1_expose_event(GtkObject *object, gpointer user_data);
void on_tool_button_go_toggled(GtkObject *object, gpointer user_data);
void on_action_info_activate(GtkObject *object, gpointer user_data);
gboolean on_drawingarea1_button_press_event(GtkWidget *widget, GdkEventButton *event, gpointer user_data);
void on_window_destroy(GtkObject *object, gpointer user_data);
#endif

void* do_png(void* delay);
void* do_jpeg(void* delay);
void set_ss_timeout(int32_t delay, const char* type);
void* do_stop(void* arg);
void* do_quit(void* arg);


