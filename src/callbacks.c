/* ReSCAL - GTK callbacks
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
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */



#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <pthread.h>
#include <sys/time.h>
//#include <gtk/gtk.h>
#include "defs.h"
#include "macros.h"
#include "entry.h"
#include "callbacks.h"
#include "interface.h"
#include "rescal.h"
//#include "monitor.h"
#include "format.h"
#ifdef PARALLEL
#include "synchro.h"
#endif
#include "view.h"
#include "cells.h"
#include "doublets.h"
#include "transitions.h"
#include "space.h"
#include "surface.h"
#include "simul.h"
#include "trace.h"
#ifdef LGCA
#include "lgca.h"
#endif

extern int H, L, D;
extern float frame_delay;
extern unsigned long int iter, md_iter;
extern unsigned char opt_quit;
extern double csp_time;
extern int abs_cv, prof_cv;
extern unsigned char opt_cv;
extern int vdir_mode; //display mode of the current orientation
extern unsigned char reorient_flag;
extern unsigned char csphpp_flag;
#ifdef LGCA
extern int use_lgca;
extern int col_iter;
extern float meanvel, maxvel;
#endif
#ifdef ROTATIONS
extern int rot_mode;
#endif
#ifdef PARALLEL
extern int mode_par;    //mode parallele
extern int proc_id;     //numero de process en mode parallele (0 = serveur)
#endif

#ifdef GUI
extern GtkWidget *drawingarea;
extern GtkWidget *window;
extern GtkLabel *label;
extern GtkStatusbar *statusbar;

gint area_w, area_h;
gint timer=0;
gint end_of_rescal=0;     //semaphore (pour l'arret des callbacks)
gint rescal_paused=0;     //semaphore (pour la mise en pause)
volatile gint dump_wait=0;         //semaphore (pour la sauvegarde)
gint mouse_flag=0;
GdkGC *gc=NULL;
#else
int end_of_rescal=0;
int rescal_paused=0;
volatile int dump_wait=0;
#endif

pthread_mutex_t mutex_pause;
pthread_cond_t cond_pause;
pthread_mutex_t mutex_display;
#ifdef PARALLEL_AUTOMATA
pthread_barrier_t lgca_barrier;
pthread_t thread_id_lgca;
#endif // LGCA

int elapsed(double *sec)
{
  //extern int gettimeofday();
  struct timeval t;
  struct timezone tz;

  int stat;
  stat = gettimeofday(&t, &tz);
  *sec = (double)(t.tv_sec + t.tv_usec/1000000.0);
  return(stat);
}

void *rescal_thread(void *data) {
  sleep(1);
  rescal();
  end_of_rescal = 1;
  pthread_exit(0);
  return 0;
}

#ifdef PARALLEL_AUTOMATA

void *lgca_thread(void *data) {
  while (1) {
    simul_lgca();
    pthread_barrier_wait(&lgca_barrier);
    //LogPrintf("LGCA-barrier\n");
  }
  pthread_exit(0);
  return 0;
}


void wait_lgca_thread()
{
  //  extern pthread_barrier_t lgca_barrier;
  //  extern pthread_t thread_id_lgca;
  static int nRetVal = -1;

  /// wait for the end of lgca cycle
  if (nRetVal == 0) {pthread_barrier_wait(&lgca_barrier); /*LogPrintf("CSP-barrier\n");*/}

  /// start lgca thread
  if (nRetVal == -1){
    pthread_barrier_init(&lgca_barrier,NULL,2);
    nRetVal = pthread_create( &thread_id_lgca, 0, lgca_thread, 0);
    if (nRetVal != 0) {ErrPrintf("ERROR: cannot create lgca_thread\n"); exit(-1);}
  }
}

#endif

void do_thread_sched()
{
  int unlocked=0;

  if (dump_wait){
    //push_status("writing",0); //do not work !?
    //LogPrintf("dump_wait=%d   unlock\n", dump_wait);
    //update_window();
    unlocked=1;
    unlock_csp(0);
    while (dump_wait) sched_yield();
  }
  if (rescal_paused){
//     push_status("paused",0);
    LogPrintf("paused\n");
    if (!unlocked){
      unlocked=1;
      unlock_csp(0);
    }
    pthread_cond_wait(&cond_pause, &mutex_pause);
    LogPrintf("resume\n");
//     pop_status(0);
  }
  if (unlocked){
    lock_csp(0);
    unlocked=0;
    //pop_status(0);
    //LogPrintf("dump_wait=%d   lock\n", dump_wait);
  }
}

void log_info()
{
  //lock_csp(0);
  dump_time();
  #ifdef INFO_CEL
  log_cell();
  #endif
  #ifdef INFO_DBL
  dump_db_info();
  #endif
  #ifdef INFO_TRANS
  dump_trans_info();
  #endif
  #ifdef PARALLEL
  if (mode_par) dump_msg_info();
  #endif
  #if defined(TRACE_TRANS) || defined(TRACE3D_CEL) || defined(TRACE_FLUX)
  trace_dump(1);
  #endif
  #ifdef LGCA
  if (use_lgca){
    dump_densite();
    #ifndef STABILITY_ANALYSIS
    dump_vel();
    #endif
    #ifdef CGV
    dump_cgv_coef();
    #endif
  }
  #endif // LGCA
  //unlock_csp(0);
}

void* do_png(void* delay) {
  while (1) {
    dump_image_inter((intptr_t) delay, "png");
    sleep((intptr_t) delay);
  }
  return 0;
}

void* do_jpeg(void* delay) {
  while (1) {
    dump_image_inter((intptr_t) delay, "jpeg");
    sleep((intptr_t) delay);
  }
  return 0;
}

void set_ss_timeout(int delay, const char* type) {
  pthread_t pth;
  if (type == "png") {
    pthread_create(&pth, 0, do_png,  (void*)(intptr_t) delay);
  } else if (type == "jpeg") {
    pthread_create(&pth, 0, do_jpeg, (void*)(intptr_t) delay);
  }
}

void* do_stop(void* arg) {
  sleep((intptr_t) arg);
  LogPrintf("time to quit !\n");
  //push_status("stopping ...", 0);
  exit(-1);
  return 0;
}

void* do_quit(void* arg) {
  sleep(1);
  if (opt_quit && end_of_rescal) {
    LogPrintf("quit\n");
    exit(-1);
  }
  return 0;
}

#ifdef GUI
void callbacks_init()
{
  view_img_size(&area_w, &area_h);

#if !defined(CSP_MUTEX) && defined(REORIENT_AUTO)
  ErrPrintf("ERROR: CSP_MUTEX option is required when REORIENT_AUTO option is set. Please define CSP_MUTEX option.\n");
  exit(-1);
#endif

  /* Initialize mutex and condition variable objects */
  pthread_mutex_init(&mutex_pause, NULL);
  pthread_cond_init (&cond_pause, NULL);
  pthread_mutex_init(&mutex_display, NULL);
}

void lock_display(int log_flag)
{
  if (log_flag) LogPrintf("lock display\n");
  pthread_mutex_lock(&mutex_display);
  if (log_flag) LogPrintf("display locked\n");
}

//non blocking function
int trylock_display(int log_flag)
{
  int ret;
  if (log_flag>0) LogPrintf("try to lock display\n");
  ret = pthread_mutex_trylock(&mutex_display);
  if (ret && (log_flag)) LogPrintf("trylock_display failed (%d)!\n",ret); //debug
  if (log_flag>0) LogPrintf("display locked\n");
  return ret;
}

void unlock_display(int log_flag)
{
  if (log_flag) LogPrintf("unlock display\n");
  pthread_mutex_unlock(&mutex_display);
  if (log_flag) LogPrintf("display unlocked\n");
}

double get_flash_delay(gpointer  data)
{
  int i;
  double t0, t1;

  elapsed(&t0);
  for (i=1; i<10; i++) gcallback_flash(data);
  elapsed(&t1);
  return (t1-t0)/10.0;
}

void timer_init() {
  //GtkWidget  *drawingarea = lookup_widget(GTK_WIDGET (widget), "drawingarea1");
  if (!frame_delay){
    float flash_delay = get_flash_delay(drawingarea) * 1000.0;
    //LogPrintf("frame_delay = %f\n", frame_delay);
    LogPrintf("flash_delay_avg = %f\n", flash_delay);
    //if (flash_delay < frame_delay) frame_delay -= flash_delay;
    //else frame_delay = flash_delay + frame_delay;
    //if (frame_delay < flash_delay / 2.0) frame_delay = flash_delay / 2.0;
    frame_delay = flash_delay * 5;
    if (frame_delay < 10) frame_delay = 10;
  }
  LogPrintf("frame_delay = %f\n", frame_delay);
  timer = gtk_timeout_add(frame_delay, gcallback_flash, 0);
}

gboolean gcallback_dump (gpointer data){
  static int cpt=0;
  //static guint id = 0;
  DumpPar  *dp;
  //if (!id) id = gtk_statusbar_get_context_id(statusbar, "gcallback_dump");
  dp = (DumpPar*)data;
  dump_wait++;
  lock_csp(0);
  //push_status(status_str[dp->type]);
  //push_status("writing",id);
#ifdef REORIENT_AUTO
  if (dp->type != DUMP_LOG){
    //reorientation_terre(0);
    reorient_flag = 1;
    rotation(0, rot_mode, ROT_REORIENT_TEMP);
    //display wind direction
    vdir_mode = VDIR_WIND;
  }
#endif
  //LogPrintf("gcallback_dump: %d\n", dp->type);
  /*if (!rescal_paused || !dp->delay)*/{
    switch(dp->type){
      case(DUMP_CSP):
      case(DUMP_BIN):
        //dump_terre(dp->type, dp->delay);
        dump_terre(dp->type, cpt, UNIT_COMP);
#ifdef ALTI //MODEL_AVA
        dump_surface("ALTI", cpt, UNIT_COMP);
#endif
#ifdef DUMP_RUGOSI
        dump_rugosi();
#endif
        break;
      case(DUMP_PNG):
        //view_dump_inter(dp->delay, IMG_PNG);
#if !defined(USE_LIBPNG) && !defined(USE_GD)
        gcallback_flash(0); //update image buffer
#endif
        dump_image_inter(dp->delay, "png");
        break;
      case(DUMP_JPEG):
        //view_dump_inter(dp->delay, IMG_JPEG);
#ifndef USE_GD
        gcallback_flash(0); //update image buffer
#endif
        dump_image_inter(dp->delay, "jpeg");
        break;
      case(DUMP_LOG):
        log_info();
        break;
    }
  }
  //sleep(3);
#ifdef REORIENT_AUTO
  if (dp->type != DUMP_LOG){
    //reorientation_terre(1);
    rotation(0, rot_mode, ROT_REORIENT_UNDO);
    vdir_mode = VDIR_NONE;
    reorient_flag = 0;
  }
#endif

  if ((dp->type == DUMP_CSP) || (dp->type == DUMP_BIN)){
#ifdef LGCA
    /// lattice gas and shear stress not reoriented
    if (use_lgca && csphpp_flag){
      dump_mvt(cpt, UNIT_COMP);
#ifdef CGV
      dump_grad_vel(cpt, UNIT_COMP);
#endif
    }
#endif
    cpt += dp->delay;
  }

  //sleep(2);
  //pop_status(id);
  dump_wait--;
  unlock_csp(0);
  //LogPrintf("gcallback_dump: done %d\n", dp->type);
  return !end_of_rescal;
}

gboolean gcallback_stop (gpointer  data)
{
  LogPrintf("time to quit !\n");
  //push_status("stopping ...", 0);
  sleep(1);
  exit(-1);
}

gboolean gcallback_quit (gpointer data)
{
  if (opt_quit && end_of_rescal) {
    LogPrintf("quit\n");
    exit(-1);
  }

  return 1;
}

/// Display text area
gboolean gcallback_text (gpointer  data)
{
  static char start=2;
  static int last_iter=0;
  static int cpt=1;
  static float speed=0;
  static double t0, t1=0, t3=0;//, e0, e1, e2;
  static float avg_delay=1e10;
  static char slabel[1024], slabel2[512];
#ifdef LGCA
  static float speed_col=0;
  static int last_col_iter=0;
#endif

  elapsed(&t0);
  if (!t1) t1 = t0;
  if (!t3) t3 = t0;

  if (t0-t1 >= 1.0) {
    float delay = t0-t1;
    avg_delay = delay/(float)cpt;
    speed = (float)(iter-last_iter)/delay;
    last_iter = iter;
    t1 = t0;
    cpt = 1;
    if (start){
      start--;
      if (!start) LogPrintf("affichage : %.2f fps\n", (1.0/avg_delay));
    }
  }
  else {
    cpt++;
  }

#ifdef LGCA
  if (use_lgca && (t0-t3 >= 3.0)) {
    float delay = t0-t3;
    speed_col = (float)(col_iter-last_col_iter)/delay;
    last_col_iter = col_iter;
    t3 = t0;
  }
#endif

  //sprintf(slabel,"Transitions = %lu%09lu \nTemps=%f \nVitesse = %d tr/min \nAffichage = %.2f img/sec", md_iter, iter, temps, (int)speed, (1.0/avg_delay));
  snprintf(slabel,sizeof(slabel),"Model : %s \nSize = %dx%dx%d \nTime = %.3f", MOD_NAME, H-2, L-2, D-2, csp_time);
#ifdef TIME_SCALE
  extern double real_time;
  /*
  snprintf(slabel2,sizeof(slabel2),"\nReal time = %.3f", real_time);
  strcat(slabel,slabel2);
  */
  long rt = (long)real_time;
  int s = rt % 60; //seconds
  rt /= 60;
  int m = rt % 60; //minuts
  rt /= 60;
  int h = rt % 24; //hours
  rt /= 24;
  int d = rt; //% 365; //days
  //rt /= 365;
  //int y = rt; //years
  snprintf(slabel2,sizeof(slabel2),"\nReal time = %dd%02dh%02dm%02ds", d, h, m, s);
  strcat(slabel,slabel2);
#endif //TIME_SCALE
  snprintf(slabel2,sizeof(slabel2),"\nTransitions = %lu%09lu \nTrans. rate = %d tr/sec \nFrame rate = %.2f fps", md_iter, iter, (int)speed, (1.0/avg_delay));
  strcat(slabel,slabel2);
#ifdef LGCA
  if (use_lgca){
    snprintf(slabel2,sizeof(slabel2),"\nLgca rate = %d cyc/sec \nMean vel. = %.3f \nMax. vel. = %.3f", (int)speed_col, meanvel, maxvel);
    strcat(slabel,slabel2);
  }
#endif
#ifdef REORIENT_AUTO
  extern float csp_angle;
  static int orientation_flag = 0;
  if (csp_angle) orientation_flag = 1;
  if (orientation_flag){
    snprintf(slabel2,sizeof(slabel2),"\nOrientation = %.2f", csp_angle);
    strcat(slabel,slabel2);
  }
#endif
  if (opt_cv){
    sprintf(slabel2,"\n(%d,%d)", abs_cv, prof_cv);
    strcat(slabel,slabel2);
  }
#ifdef PARALLEL
  if (mode_par){
    sprintf(slabel2,"\nPid = %d", proc_id);
    strcat(slabel,slabel2);
  }
#endif
  //snprintf(slabel2,sizeof(slabel2),"\n___________________________________");
  //strcat(slabel,slabel2);
  gtk_label_set_text(label, slabel);
  //gdk_gc_unref (gc);

  return TRUE;
}

/// Display image area
gboolean gcallback_flash (gpointer  data)
{
  static GdkRgbCmap colormap;
  unsigned char *image_zoom;
  GdkDrawable *drawable;
  unsigned char *image;
  int w, h;

  //elapsed(&e0);
  //GtkWidget  *drawingarea = (GtkWidget  *)data;
  if (!drawingarea) return TRUE;

  drawable = drawingarea -> window;
  if (!drawable) return TRUE;

  if (!gc) {
    //LogPrintf("!gc\n");
    gdk_rgb_init();
    gc = gdk_gc_new (drawable);
    view_palette((int*)colormap.colors);
  }

  if (trylock_display(0/*-1*/)) return TRUE;
  image = view();
//#ifdef GTK_OLD_SYNTAX
#if GTK2 && (GTK_MINOR_VERSION<24)
  gdk_drawable_get_size(drawable, &w, &h);
#else
  w = gdk_window_get_width(drawable);
  h = gdk_window_get_height(drawable);
#endif
  if ((w == area_w) && (h == area_h)){
    reset_zoom();
    gdk_draw_indexed_image (drawable, gc, 0, 0, area_w, area_h, GDK_RGB_DITHER_NONE, image, area_w, &colormap);
    //LogPrintf("coef=1.0\n");
  }
  else{
    float wcoef = (float)w/area_w;
    float hcoef = (float)h/area_h;
    float coef = Min(wcoef,hcoef);
    //LogPrintf("coef=%f, w=%d, h=%d\n", coef, w, h);
    /*if (coef>=1*/{
      image_zoom = view_zoom(w, h, coef);
      gdk_draw_indexed_image (drawable, gc, 0, 0, w, h, GDK_RGB_DITHER_NONE, image_zoom, w, &colormap);
    }
    //LogPrintf("coef=%f\n",coef);
  }

  unlock_display(0);

  gcallback_text(data);

  return TRUE;
}

void on_drawingarea1_expose_event(GtkObject *object, gpointer user_data)
{
  if(!timer) {
    timer_init();
  }

  if (area_w) {
    gcallback_flash(0);
  }
}

/*
void on_tool_button_go_toggled (GtkObject *object, gpointer user_data)
{

  rescal_paused = !gtk_toggle_tool_button_get_active (GTK_TOGGLE_TOOL_BUTTON(object));
  LogPrintf(rescal_paused ? "pause\n" : "continue\n");
  if (rescal_paused) push_status("paused"); else pop_status();
  if (!rescal_paused) pthread_cond_signal(&cond_pause);
}
*/
void on_tool_button_go_clicked (GtkObject *object, gpointer user_data)
{
  rescal_paused = 1 - rescal_paused;
  LogPrintf(rescal_paused ? "pause\n" : "continue\n");
  //if (rescal_paused) push_status("paused"); else pop_status();
  if (!rescal_paused) pthread_cond_signal(&cond_pause);
  gtk_tool_button_set_stock_id (GTK_TOOL_BUTTON(object), rescal_paused ? "gtk-media-play" : "gtk-media-pause");
}

void on_action_info_activate (GtkObject *object, gpointer user_data)
{
  static char active = 0;

  active = 1-active;
  LogPrintf("on_action_info_activate : %d\n", active);

  if (active){
    gtk_widget_show(GTK_WIDGET(label));
  }
  else{
    gtk_widget_hide(GTK_WIDGET(label));
    gtk_window_resize(GTK_WINDOW(window), area_w, area_h);
  }
}

void on_action_snapshot_activate (GtkObject *object, gpointer user_data)
{
  DumpPar dp_img;
  dp_img = (DumpPar){DUMP_PNG, 0};
  gcallback_dump(&dp_img);
}

gboolean on_drawingarea1_button_press_event(GtkWidget *widget, GdkEventButton *event, gpointer user_data)
{
  LogPrintf("Event: button pressed\n"); fflush(stdout);
  mouse_flag = update_cv(event->x, event->y, 0);

  return FALSE;
}

gboolean on_drawingarea1_motion_notify_event(GtkWidget *widget, GdkEventButton *event, gpointer user_data)
{
  if (mouse_flag){
    update_cv(event->x, event->y, mouse_flag);
  }
  return FALSE;
}

gboolean on_drawingarea1_button_release_event(GtkWidget *widget, GdkEventButton *event, gpointer user_data)
{
  //LogPrintf("Event: button released\n"); fflush(stdout);
  mouse_flag = 0;
  return FALSE;
}

void on_window_destroy (GtkObject *object, gpointer user_data)
{
  LogPrintf("quit\n");
  if (gc) gdk_gc_unref (gc);
  gtk_main_quit();
#ifdef PARALLEL
  if (mode_par) synchro_quit();
#endif
  view_quit();
  //rescal_quit();
}

#endif //GUI

