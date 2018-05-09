/* ReSCAL - GTK interfaces
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


#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <string.h>

#include <gdk/gdkkeysyms.h>
#include <gtk/gtk.h>

#include "defs.h"
#include "macros.h"
#include "callbacks.h"
#include "interface.h"

extern gint area_w, area_h;
extern unsigned char opt_info;
extern unsigned char opt_nv;

GtkWidget *drawingarea;
GtkStatusbar *statusbar;
guint id_status=0;

#ifdef GUI

GtkWidget *window=NULL;
GtkLabel *label;

GtkWidget* create_window(void)
{
  GtkBuilder  *builder;
  //GtkTextView       *textview;
  //GtkTextBuffer     *buffer;
  gchar *ui_file = "rescal-ui.xml";

  gchar* filename = g_build_filename (RESCAL_DATA_DIR, ui_file, NULL);

  LogPrintf("create_window: start\n");

  builder = gtk_builder_new();
  if (!gtk_builder_add_from_file(builder, filename, NULL)){
    //try relative path
    filename = g_build_filename ("../src", ui_file, NULL);
    if (!gtk_builder_add_from_file(builder, filename, NULL)){
      ErrPrintf("ERROR: cannot load : %s\n", filename);
      exit(-1);
    }
  }

  window = GTK_WIDGET(gtk_builder_get_object (builder, "window1"));
  gtk_window_resize(GTK_WINDOW(window), area_w, area_h);

  //textview = (GtkTextView*) gtk_builder_get_object (builder, "textview1");
  //gtk_text_view_set_buffer (textview, "click on Go button ...");

  drawingarea = GTK_WIDGET(gtk_builder_get_object (builder, "drawingarea1"));
  gtk_widget_set_usize (drawingarea, area_w, area_h);

  label = GTK_LABEL(gtk_builder_get_object (builder, "label1"));
  gtk_widget_hide(GTK_WIDGET(label));

  statusbar = GTK_STATUSBAR(gtk_builder_get_object (builder, "statusbar1"));

  gtk_builder_connect_signals(builder, NULL);

  if (opt_info) gtk_toggle_tool_button_set_active (GTK_TOGGLE_TOOL_BUTTON(gtk_builder_get_object (builder, "tool_button_info")), TRUE);

  g_object_unref(G_OBJECT (builder));

  gtk_widget_show(window);

  gtk_label_set_text(label, "click on Go button ...");

  id_status = gtk_statusbar_get_context_id(statusbar, "rescal");
  push_status("ready",id_status);

  LogPrintf("create_window: end\n");

  return window;
}


void update_window()
{
  if (window && id_status){
    gdk_threads_enter ();
    /*while (gtk_events_pending ())*/
    gtk_main_iteration ();
    gdk_threads_leave ();
  }
}

#endif //GUI

void push_status(char* str, int id)
{
  if (opt_nv) return;
  if (!id) id = id_status;
  if (id){
    //gchar* status = g_strdup_printf("%s | %s", MOD_NAME, str);
    gdk_threads_enter ();
    gtk_statusbar_push (statusbar, id, str);
    gdk_threads_leave ();
    //LogPrintf("push_status: id = %d   str = %s\n", id, str);
  }
}

void pop_status(int id)
{
  if (opt_nv) return;
  if (!id) id = id_status;
  if (id){
    gdk_threads_enter ();
    gtk_statusbar_pop (statusbar, id);
    gdk_threads_leave ();
    //LogPrintf("pop_status: id = %d\n", id);
  }
}

