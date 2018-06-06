/* ReSCAL - Main entry
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



#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#define _MAIN_
#include <locale.h>
//GTK-JBK
//#include <gtk/gtk.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include <unistd.h>
#include <stdint.h>

#include "defs.h"
#include "macros.h"
#include "entry.h"
#include "param.h"
#include "cells.h"
#include "simul.h"
//GTK interface-JBK
//#include "interface.h"
#include "callbacks.h"
//#include "monitor.h"
#include "view.h"
#include "rescal.h"
#include "space.h"
#include "format.h"
#ifdef PARALLEL
#include "messages.h"
#endif


extern int32_t graine;            // graine pour la generation des nombres aleatoires
extern float dump_delay_png;
extern float dump_delay_csp;
extern uint8_t csphpp_flag;
extern float stop_delay_t0;
extern double stop_time;
extern int32_t mode_pat;
#ifdef PARALLEL
extern int32_t mode_par;    //mode parallele
extern int32_t proc_id;
#endif

int32_t prog=PROG_RESCAL;
int32_t arch=0;
float frame_delay=0;
uint8_t opt_h=0, opt_hm=0, opt_nv=0, opt_info=0,  opt_dcsp=0, opt_dbin=0, opt_stop=0, opt_quit=0;
uint8_t flag_no_input=0;
uint8_t opt_dpng=0, opt_djpeg=0;

//converted from gint32_t to int-JBK
//gint32_t dcsp_delay, dbin_delay, stop_delay;
//gint32_t dpng_delay=0, djpeg_delay=0;
int32_t dcsp_delay, dbin_delay, stop_delay;
int32_t dpng_delay=0, djpeg_delay=0;

int32_t main (int32_t argc, int8_t *argv[])
{
  DumpPar dp_info, dp_csp, dp_bin, dp_png, dp_jpeg;

#ifdef PARALLEL
  if (argc>1) msg_init();
#endif

#ifdef LOG_FILE
  int8_t log_filename[100];
  sprintf(log_filename,"RESCAL.log");
#ifdef PARALLEL
  if (proc_id>=0){
    sprintf(log_filename,"RESCAL%d.log",proc_id);
  }
#endif
  if (argc>1) log_file = fopen(log_filename,"w");
#endif

  //LogPrintf("ReSCAL\n");
  LogPrintf("ReSCAL %s (%s)\n",VER_NUM,VER_DAT);
#ifdef VER_BRA
  LogPrintf("branch %s\n",VER_BRA);
#endif
  LogPrintf("%s model", MOD_NAME);
#ifdef MOD_DESC
  LogPrintf(" - %s", MOD_DESC);
#endif
  LogPrintf("\n");

  if (argc == 1){
    usage_no_input();
    exit(-4);
  }

  flag_no_input = (*argv[1] == '-');


  //GTK to be removed-JBK
  //gtk_set_locale();

#ifdef NUM_MODE
  LogPrintf("format of numerical values: %s \n", NUM_MODE);
  setlocale(LC_NUMERIC, NUM_MODE);
#endif
//#ifdef GTK2
//  LogPrintf("compiled with GTK2.%d\n", GTK_MINOR_VERSION);
//#endif
  if (sizeof(long)==4) arch = W32;
  if (sizeof(long)==8) arch = W64;
  LogPrintf("size of long: %ld bytes (%s)\n", sizeof(long), (arch==W32)?"32b":"64b");
  LogPrintf("size of cell: %lu byte(s)\n", sizeof(Cell));

//this isn't accesible with the GTK free version-JBK
#ifndef GUI
  LogPrintf("graphical interface disabled\n");
#endif


#ifdef USE_GD
  LogPrintf("use of libgd\n");
#endif
#ifdef LGCA
  LogPrintf("lattice gas enabled\n");
#endif
#ifdef CYCLAGE_HOR
  LogPrintf("support for periodic boundary conditions\n");
#endif
#ifdef OPENMP
  LogPrintf("OpenMP enabled\n");
#endif
#ifdef USE_AVX
  LogPrintf("AVX enabled\n");
#endif
#ifdef PARALLEL_AUTOMATA
  //LogPrintf("\"pat\" multithreading supported\n");
#endif
#ifdef PARALLEL
  LogPrintf("parallel execution with MPI\n");
#endif
#ifdef __CYGWIN__
  LogPrintf("compiled under Cygwin\n");
#endif

  general_options(argc, argv);

  if (opt_h || opt_hm){
    rescal_usage();

    if (opt_h){
      show_general_options();
      show_view_options();
    }

    if (opt_hm){
      rescal_params();
      param_usage();
    }

    printf("\nPlease see README file for instructions ");
    printf("(also see http://www.ipgp.fr/rescal).\n");
  }

  if (flag_no_input){
    if (!(opt_h || opt_hm)) usage_no_input();
    exit(-4);
  }


  /// is execution deterministic ?
  int32_t dflag = 1;
#ifdef DETERMINISTIC
#ifdef PARALLEL_AUTOMATA
  if (mode_pat) dflag = 0;
#endif // PARALLEL_AUTOMATA
#ifdef PARALLEL
  if (mode_par) dflag = 0;
#endif // PARALLEL
#else
#if defined(LGCA) && defined(OPENMP)
  dflag = 0;
#endif // LGCA
#endif // DETERMINISTIC
  if (dflag)
    LogPrintf("deterministic execution\n")
  else
    LogPrintf("non-deterministic execution\n");

//Major GTK section to be removed-JBK
  /*/
  if (!opt_nv){
#if (GLIB_MAJOR_VERSION <= 2) && (GLIB_MINOR_VERSION<32)
    g_thread_init(NULL); /// deprecated function since glib 2.32
#endif
    gdk_threads_init();
    gtk_init (&argc, &argv);
  }
   */
  /*else gtk_init_check (&argc, &argv);*/

  rescal_init(argv[1]);
  view_init(argc, argv);
#ifdef GUI
  callbacks_init();
#endif

  //if (opt_nv) rescal(0);

//This is some more GTK code that isn't needed-JBK
/*
#ifdef GUI
  if (!opt_nv) {
    create_window();
    gcallback_flash(0);
  }
#endif
*/

  if (opt_info) {
    log_info();
    dp_info = (DumpPar){DUMP_LOG, 1};
    //GTK output-JBK
    /*
    g_timeout_add_seconds(60, gcallback_dump, &dp_info); //toutes les minutes
    //g_timeout_add_seconds(1, gcallback_dump, &dp_info); //toutes les secondes
  */
   }

//Removed for GTK dependinces-JBK
//TODO Remove GTK Dependinces by channging callback dump, and g_timeout-JBK
/*  if (opt_dcsp) {
    //dump_terre(DUMP_CSP, dcsp_delay);
    dp_csp = (DumpPar){DUMP_CSP, dcsp_delay};
    gcallback_dump(&dp_csp);
    g_timeout_add_seconds(dcsp_delay*60, gcallback_dump, &dp_csp);
  }

  if (opt_dbin) {
    //dump_terre(DUMP_BIN, dbin_delay);
    dp_bin = (DumpPar){DUMP_BIN, dbin_delay};
    gcallback_dump(&dp_bin);
    g_timeout_add_seconds(dbin_delay*60, gcallback_dump, &dp_bin);
  }
*/
#if defined(USE_LIBPNG) || defined(USE_GD)
  view_dump_init();
#endif

  if (opt_dpng) {
    set_ss_timeout(dpng_delay, "png");
    //g_timeout_add_seconds(dpng_delay, gcallback_dump, &dp_png);
  }
  if (opt_djpeg) {
    set_ss_timeout(djpeg_delay, "jpeg");
    //g_timeout_add_seconds(djpeg_delay, gcallback_dump, &dp_jpeg);
  }

  if (opt_stop) {
    pthread_t pth;
    pthread_create(&pth, 0, do_stop, (void*)(intptr_t)(stop_delay*60));
    //TODO remove GTK dependince from timeout_add-JBK
    //g_timeout_add_seconds_full(G_PRIORITY_DEFAULT, stop_delay*60+15, gcallback_stop, 0, NULL);

    //gtk_timeout_add(stop_delay*60000+15000, callback_stop, 0);
  }
  if (opt_quit) {
    pthread_t pth;
    pthread_create(&pth, 0, do_quit, 0);
    //TODO remove GTK dependince from timeout_add-JBK
    //g_timeout_add_seconds_full(G_PRIORITY_DEFAULT, 1, gcallback_quit, 0, NULL);
  }

  int32_t nRetVal;
  pthread_t pth;
  nRetVal = pthread_create( &pth, 0, rescal_thread, 0);
  pthread_join(pth,0); //TODO: THIS IS TEMP CODE-JBK
  //**** This was already commented out? ***-JBK
  /*
  GThread *pth;
  GError *err = NULL;
  g_thread_init(NULL);
  pth = g_thread_create((GThreadFunc)rescal_thread, 0, FALSE, &err);
  */

  //gtk_idle_add(idle, NULL);
  //*************************************

  //More GTK that's not needed-JBK
  //gtk_main ();

#ifdef LOG_FILE
  fclose(log_file);
#endif

  return 0;
}

void read_delay(int8_t *str, DumpDelay *pd)
{
  int32_t n;
  if (strstr(str, "t0")){
    pd->unit = UNIT_T0;
    if (sscanf(str, "%ft0", &pd->val)!=1) {ErrPrintf("ERROR: invalid delay value %s\n",str); exit(-1);}
  }
  else{
    pd->unit = UNIT_COMP;
    if (sscanf(str, "%d", &n)!=1) {ErrPrintf("ERROR: invalid delay value %s\n",str); exit(-1);}
    pd->val = n;
  }
  //LogPrintf("delay: %f (%d)\n", pd->val, pd->unit);
}

void general_options(int32_t argc, int8_t *argv[])
{
  int32_t i;
  DumpDelay dd;

  /// log the command line options
  if (!flag_no_input){
    LogPrintf("args: ");
    for(i=1; i<argc; i++) LogPrintf("%s ", argv[i]);
    LogPrintf("\n");
  }

  /// parse the command line options
  for(i=1; i<argc; i++){
    if (!strcmp(argv[i],"-h") || !strcmp(argv[i],"--help"))
      opt_h = 1;
    else if (!strcmp(argv[i],"-hm"))
      opt_hm = 1;
    else if (!strcmp(argv[i],"-nv"))
      opt_nv = 1;
    else if (!strcmp(argv[i],"-fr")){
      int32_t frame_rate = atoi(argv[++i]);
      frame_delay = (frame_rate<100) ? 1000/frame_rate : 10;
    }
    else if (!strcmp(argv[i],"-info")){
      opt_info = 1;
#ifdef PARALLEL
      if (proc_id > 0) opt_info = 0;
#endif
    }
    else if (!strcmp(argv[i],"-dcsp") || !strcmp(argv[i],"-dcsphpp")){
      if (!strcmp(argv[i],"-dcsphpp")) csphpp_flag = 1;
      //opt_dcsp = 1;
      //dcsp_delay = atoi(argv[++i]);
      read_delay(argv[++i], &dd);
      if (dd.unit==UNIT_COMP){
        opt_dcsp = 1;
        dcsp_delay = dd.val;
      }
      else if (dd.unit==UNIT_T0){
        dump_delay_csp = dd.val;
      }
    }
    else if (!strcmp(argv[i],"-dbin")){
      opt_dbin = 1;
      dbin_delay = atoi(argv[++i]);
    }
    else if (!strcmp(argv[i],"-dpng")){

      //dpng_delay = atoi(argv[++i]);
      read_delay(argv[++i], &dd);
      if (dd.unit==UNIT_COMP){
        opt_dpng = 1;
        dpng_delay = dd.val;
      }
      else if (dd.unit==UNIT_T0){
        dump_delay_png = dd.val;
      }
    }
#ifndef USE_LIBPNG
    else if (!strcmp(argv[i],"-djpeg")){
      opt_djpeg = 1;
      djpeg_delay = atoi(argv[++i]);
      if (!djpeg_delay) djpeg_delay = 1;
    }
#endif
    else if (!strcmp(argv[i],"-stop")){
      read_delay(argv[++i], &dd);
      if (dd.unit==UNIT_COMP){
        opt_stop = 1;
        stop_delay = dd.val; //atoi(argv[++i]);
      }
      else if (dd.unit==UNIT_T0){
        stop_delay_t0 = dd.val + 0.5;
      }
      LogPrintf("stop delay: %s\n", argv[i]);
    }
    else if (!strcmp(argv[i],"-stoptime")){
      read_delay(argv[++i], &dd);
      if (dd.unit==UNIT_T0){
        stop_time = dd.val + 0.5;
      }
      else if (dd.unit==UNIT_COMP){
        ErrPrintf("ERROR: bad unit in -stoptime option\n");
        exit(-1);
      }
      LogPrintf("stop time: %s\n", argv[i]);
    }
    else if (!strcmp(argv[i],"-q")){
      opt_quit = 1;
    }
    else if (!strcmp(argv[i],"-s") || !strcmp(argv[i],"-g")){
      graine = atoi(argv[++i]);
    }
#ifdef PARALLEL_AUTOMATA
    else if (!strcmp(argv[i],"-pat")){
      mode_pat = 1;
      LogPrintf("parallel automata with threads\n");
    }
#endif
#ifdef PARALLEL
    else if (!strcmp(argv[i],"-par")){
      mode_par = 1;
    }
#endif
    /*else if (*argv[i]=='-'){
      ErrPrintf("ERROR: unknown option %s\n", argv[i]);
      exit(-5);
    }*/
  }
  if (opt_dcsp && opt_dbin){
    ErrPrintf("ERROR: -dcsp and -dbin options are mutually exclusive\n");
    exit(-1);
  }
#ifndef GUI
  opt_nv = 1;
#endif
}

void show_general_options()
{
  printf("\nGENERAL OPTIONS\n");
  printf("  -h \t\t display usage information\n");
  printf("  -hm \t\t display parameters of the compiled model\n");
#ifdef GUI
  printf("  -nv \t\t no window mode\n");
#endif
  printf("  -info \t monitoring data every minutes\n");
  printf("  -dcsp <n> \t generation of a CSP file every <n> minutes\n");
  printf("  -dcsp <f>t0 \t generation of a CSP file with delay in t0 unit (float value)\n");
  printf("  -dcsphpp <n> \t generation of CSP and HPP files every <n> minutes\n");
  printf("  -dcsphpp <f>t0 \t generation of CSP and HPP files with delay in t0 unit (float value)\n");
  printf("  -dbin <n> \t generation of a binary file every <n> minutes\n");
  printf("  -dpng <n> \t generation of a PNG image every <n> seconds\n");
  printf("  -dpng <f>t0\t generation of PNG images with delay in t0 unit (float value)\n");
#ifndef USE_LIBPNG
  printf("  -djpeg <n> \t generation of a Jpeg image every <n> seconds\n");
#endif
  printf("  -stop <n> \t quit after <n> minutes\n");
  printf("  -stop <f>t0 \t quit after running <f> t0 in simulated time\n");
  printf("  -stoptime <f>t0 \t quit when simulated time reaches <f> t0\n");
  printf("  -q \t\t quit at the end of simulation\n");
#ifdef GUI
  printf("  -fr <n> \t frame rate (fps)\n");
#endif
  printf("  -s <n> \t random seed\n");
#ifdef PARALLEL_AUTOMATA
  printf("  -pat \t\t parallel computation of lattice-gas and stochastic automata with threading\n");
#endif
//#ifdef MODEL_DUN
//  printf("  -ro <n> : rotation d'angle aleatoire toutes les <n> minutes\n");
//#endif
#ifdef PARALLEL
  printf("  -par \t\t parallel mode enabled by domain decomposition using MPI\n");
  //printf("  -pars : mode parallele, process serveur\n");
  //printf("  -parc <i> : mode parallele, process client numero <i>\n");
#endif
  fflush(stdout);
}

void usage_no_input()
{
  printf("no input file\n");
  printf("Try '-h' option for usage\n");
}


//TODO convert this to non GTK?-JBK
/*
//#ifdef GTK_OLD_SYNTAX
#if GTK2 && (GTK_MINOR_VERSION<14)
void g_timeout_add_seconds(int32_t nsec, GtkFunction gcallback_dump, gpointer data)
{
  gtk_timeout_add(nsec*1000, gcallback_dump, data);
}

void g_timeout_add_seconds_full(int32_t prior, int32_t nsec, GtkFunction gcallback_dump, gpointer data, void *notify)
{
  gtk_timeout_add(nsec*1000, gcallback_dump, data);
}
#endif
*/

