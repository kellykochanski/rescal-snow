/* ReSCAL - Stochastic simulation engine
 *
 * Copyright (C) 2011-2013
 *
 * Author: Olivier Rozier <rozier@ipgp.fr>
 *
 * Code based on dissol program,
 * by Eduardo Sepulveda <edo@espci.fr>
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


#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <memory.h>
#include <math.h>
#include <stdint.h>

// For output
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "defs.h"
#include "macros.h"
#include "param.h"
#include "format.h"
#include "cells.h"
#include "space.h"
#include "surface.h"
#include "doublets.h"
#include "transitions.h"
#include "simul.h"
#include "trace.h"
#include "lgca.h"
#include "callbacks.h"
#include "view.h"
extern uint8_t opt_info, opt_nv;
extern int32_t H, L, D, HL, HLD;       // les dimensions de la terre
extern int32_t LN, LS, LNS, HLN;    //couloir est-ouest (limite nord, limite sud, largeur nord-sud, ...)
extern Cell *TE;            // la 'terre'
extern int32_t *db_pos[];                   // les tableaux contenant la position des doublets actifs
extern Doublet t_doub[MAX_DB];               // la description des doublets
extern int32_t Ncel[];                       // nombre de cellules par type
extern int32_t Ndb[];                        // nombre doublets actifs par type
extern int32_t nb_trans_db, nb_trans_cel, nb_trans; // nombre de transitions de 2 cellules, d'1 cellule et au total
extern TransitionDb t_trans[];            // la description des transitions de doublets
extern TransitionCel t_trans_cel[];           // la description des transitions de cellules
extern LienTransDb   t_lien[];        // la description des liens de correlation entre transitions de 2 cellules
extern int32_t nb_liens;             // nombre de liens entre transitions
extern const char *classes_db[];        // noms des classes de doublet
extern const char *etats[];             // les noms des types de cellules
#ifdef INFO_TRANS
extern int32_t cpt_lien[];      // compteur des correlations effectives
#endif //INFO_TRANS
extern int32_t direction[];             // conversion classe de transition -> direction
extern char *rot_map;        // periodic mapping of the rotating space
extern Pos2 *rot_map_pos;        // periodic mapping of the rotating space
extern int32_t use_lgca;
#ifdef LGCA
extern int32_t lgca_reset;
extern int32_t lgca_speedup;
extern int32_t col_iter;    //nombre de cycles de collisions
extern int32_t CLEO;
extern char *nom_fic_mvt;
extern float maxvel, meanvel;
#ifdef CGV
extern float grdvc_max; //seuil max pour le controle de l'erosion en fonction du gradient de vitesse
extern float grdvc_min; //seuil min pour le controle de l'erosion en fonction du gradient de vitesse
#endif //CGV
#endif //LGCA
extern float lambda_A;
extern float lambda_A_unstable;
#ifdef REFDB_PTR
extern RefDoublets *RefDB;       // references des cellules de la terre vers les doublets actifs
#else
extern RefDoublets_Type *RefDB_Type;       // references des cellules de la terre vers les doublets actifs
extern RefDoublets_Ind *RefDB_Ind;       // references des cellules de la terre vers les doublets actifs
#endif //REFDB_PTR
#ifdef CGV
extern float cgv_coef;
#endif 
extern int32_t vdir_mode; //display mode of the current orientation
extern uint8_t reorient_flag;
uint64_t iter = 0, md_iter = 0; // nombre d'iterations (unites, milliards d'unites)
double init_time = 0.0;               // initial time
double csp_time = 0.0;                // time variable for the cellular space
double stop_time = 0.0;               // stopping time
// temps reel simule
#ifdef TIME_SCALE
double real_time = 0.0;     // real time with time scale
#endif
char end_of_simul = 0;
int32_t cpt_trans_blocked[MAX_TRANSITIONS_DB];  //nombre de transitions blockees
double PII, PII_time;     //probabilités totales
double PI_ij[MAX_TRANSITIONS];                  // les PI_ij
double cum[MAX_TRANSITIONS];     //probabilites cumulees
char *boundary_str = NULL;
#ifdef CYCLAGE_HOR
int32_t boundary = BC_PERIODIC;  //boundary conditions
int32_t pbc_mode = 1; //periodic boundary conditions
#else
int32_t boundary = BC_REINJECTION;  //boundary conditions
int32_t pbc_mode = 0; //periodic boundary conditions
#endif
int32_t ava_mode = AVA_SYNC;  //mode d'avalanches
char *ava_mode_str = NULL;
int32_t ava_trans = 0;        //transitions d'avalanches
int32_t ava_norm = 0;       //avalanches avec normales
int32_t ava_h_lim = 1;  // hauteur limite avant avalanche
int32_t ava_nb_cel_max = 1;  // nb de cellules qui tombent simultanement
float ava_delay = 0.0; //delay between avalanches
float ava_duration = 1.0; //default value for DUN model
float ava_angle = 0.0; //angle of avalanches (degrees)
float ava_angle_col = 0.0; //angle of avalanches (degrees) for colored grains
float ava_angle_stable = 0.0;
float ava_angle_unstable = 0.0;
float ava_angle_stable_col = 0.0;
float ava_angle_unstable_col = 0.0;
float lambda_A_stable = 0.0;
int32_t ava_upwind = 1;
uint8_t simul_dump_flag = 0;
uint8_t csphpp_flag = 0;
uint8_t alti_only_flag; // flag that causes heightmap (ALTI*) files to be written to a file, but not cellspace files (*.csp)
float dump_delay_png = 0.0;
float dump_delay_csp = 0.0;
float stop_delay_t0 = 0.0;
int32_t mode_pat;
int32_t rot_mode = ROT_MODE_OVERLAP;  //default rotation mode
#ifdef LGCA
int32_t rot_mvt_nb_cycles = -1;  //nombre de cycles de stabilisation apres rotation
#endif
#ifdef ROTATIONS
float rot_angle = 0.0;
float rot_delay = 0.0;
char *rot_mode_str = NULL;
#endif
#ifdef CENTERING_AUTO
float centering_delay = 0.0;
#endif
#ifdef LGCA
int32_t init_mvt_nb_cycles = -1;  //nombre de cycles de stabilisation du flux au depart
double lgca_delay = 0.0; //delai entre deux cycles de gaz sur reseau
int32_t lgca_ready = 0;
int32_t lgca_stabilize_flag = 0;
#endif
char *wind_data_filename = NULL;
#ifdef CGV
char *qsat_data_filename = NULL;
#endif
#ifdef TIME_SCALE
char *phys_prop_filename = NULL;
float time_scale = 1.0; //t_0
#endif

char *output_directory = "../out"; // directory for log and output files -KK
void output_headers(); // write headers for output files into output_directory
void output_write(char *output_filename, char *output_content); // write a line of output
void output_path(char *filename); // append directory name and .log to filename
void output_path_png(char *pngname); //append directory name and .png

void compute_prob_dist();
int32_t simul_trans();
int32_t simul_check(int32_t, int32_t, int32_t);
void simul_lien_trans(int32_t, int32_t, int32_t);
#ifdef LGCA
void flow_stabilization(int32_t);
#endif
void simul_stop();
void simul_dump();
void params_simul() {
  //csp_time
  parameter("Time", "initial physical time - optional", &init_time, PARAM_DOUBLE, "GENERAL");
#ifdef CYCLAGE_HOR
#if defined(MODEL_DUN) || defined(MODEL_SNO)
  parameter("Boundary", "boundary conditions (PERIODIC|OPEN|OUT|CLOSE|REINJECTION), PERIODIC by default", &boundary_str, PARAM_STRING, "GENERAL");
#else
  parameter("Boundary", "boundary conditions (PERIODIC|OPEN|CLOSE), PERIODIC by default", &boundary_str, PARAM_STRING, "GENERAL");
#endif
#elif defined(MODEL_DUN)
  parameter("Boundary", "boundary conditions (OPEN|REINJECTION), REINJECTION by default", &boundary_str, PARAM_STRING, "GENERAL");
#else
  parameter("Boundary", "boundary conditions (OPEN|CLOSE), CLOSE by default", &boundary_str, PARAM_STRING, "GENERAL");
#endif //CYCLAGE_HOR
#ifdef AVALANCHES
  param_family("AVALANCHES", "Avalanches parameters");
  parameter("Ava_mode", "mode of avalanches (NONE|SYNC|TRANS|PROPAG, SYNC by default)", &ava_mode_str, PARAM_STRING, "AVALANCHES");
  parameter("Lambda_A", "rate of avalanches in TRANS mode", &lambda_A, PARAM_FLOAT, "AVALANCHES");
  parameter("Lambda_A_unstable", "unstable avalanches rate in TRANS mode", &lambda_A_unstable, PARAM_FLOAT, "AVALANCHES");
  parameter("Ava_h_lim", "height limit for avalanches (1 by default)", &ava_h_lim, PARAM_INT, "AVALANCHES");
  parameter("Ava_nb_cel", "number of cells falling together in SYNC mode (1 by default) - optional", &ava_nb_cel_max, PARAM_INT, "AVALANCHES");
  parameter("Ava_delay", "delay between avalanches - optional", &ava_delay, PARAM_FLOAT, "AVALANCHES");
  parameter("Ava_angle", "angle of repose (degrees) - optional", &ava_angle, PARAM_FLOAT, "AVALANCHES");
  parameter("Ava_duration", "duration of avalanches - optional", &ava_duration, PARAM_FLOAT, "AVALANCHES");
  parameter("Ava_angle_unstable", "unstable angle of repose (degrees) - optional", &ava_angle_unstable, PARAM_FLOAT, "AVALANCHES");
  parameter("Ava_upwind", "enable upwind avalanches in TRANS mode (YES|NO, YES by default) - optional", &ava_upwind, PARAM_BOOLEAN, "AVALANCHES");
#endif //AVALANCHES
#ifdef ROTATIONS
  param_family("ROTATION", "Rotations parameters");
  parameter("Rot_mode", "rotation mode (FULL|DISK|OVERLAP) - optional", &rot_mode_str, PARAM_STRING, "ROTATION");
  parameter("Rot_angle", "rotation angle - optional", &rot_angle, PARAM_FLOAT, "ROTATION");
  parameter("Rot_delay", "rotation delay - optional", &rot_delay, PARAM_FLOAT, "ROTATION");
#ifdef LGCA
  parameter("Rot_ncycl", "number of flow cycles without transitions after rotation - optional", &rot_mvt_nb_cycles, PARAM_INT, "LGCA");
#endif
#endif
#ifdef CENTERING_AUTO
  parameter("Centering_delay", "automatic centering delay - optional", &centering_delay, PARAM_FLOAT, "");
#endif
#ifdef LGCA
  parameter("Init_ncycl", "initial number of flow cycles without transitions - optional", &init_mvt_nb_cycles, PARAM_INT, "LGCA");
  parameter("Lgca_delay", "physical time delay between flow cycles - optional", &lgca_delay, PARAM_DOUBLE, "LGCA");
#endif
  param_family("REAL", "Real-space parameters");
#ifdef CGV
  parameter("Tau_min", "min shear stress value in the linear relationship between erosion rate and shear stress - optional", &grdvc_min, PARAM_FLOAT, "MODEL");
  parameter("Tau_max", "max shear stress value in the linear relationship between erosion rate and shear stress - optional", &grdvc_max, PARAM_FLOAT, "MODEL");
  parameter_hidden("Grdvc_min", "Grdvc_min parameter is obsolete - Please use Tau_min parameter instead !", &grdvc_min, PARAM_FLOAT, "MODEL");
  parameter_hidden("Grdvc_max", "Grdvc_max parameter is obsolete - Please use Tau_max parameter instead !", &grdvc_max, PARAM_FLOAT, "MODEL");
  parameter("Qsat_file", "real PDF data between Qsat and tau1 - optional", &qsat_data_filename, PARAM_STRING, "REAL");
#endif //CGV
#ifdef WIND_DATA
  parameter("Wind_file", "filename for wind data - optional", &wind_data_filename, PARAM_STRING, "REAL");
#endif
#ifdef TIME_SCALE
  parameter("Phys_prop_file", "filename for physical parameters - optional", &phys_prop_filename, PARAM_STRING, "REAL");
#endif
  parameter("Output_directory", "name of directory for output - default value OUT", &output_directory, PARAM_STRING, "GENERAL");
}
void simul_parse() {
  if (boundary_str) {
    if (!strcmp(boundary_str, "PERIODIC")) {
#ifdef CYCLAGE_HOR
      boundary = BC_PERIODIC;
#else
      ErrPrintf("Incorrect value for Boundary : %s, please define CYCLAGE_HOR in defs.h and recompile\n", boundary_str);
      exit(-2);
#endif
    } else if (!strcmp(boundary_str, "OPEN")) {
      boundary = BC_OPEN;
    } else if (!strcmp(boundary_str, "OUT")) {
      boundary = BC_OUT;
    } else if (!strcmp(boundary_str, "CLOSE")) {
      boundary = BC_CLOSE;
    }
    else if (!strcmp(boundary_str, "REINJECTION")){
#if defined(MODEL_DUN) || defined(MODEL_SNO)
      boundary = BC_REINJECTION;
#else
      ErrPrintf("Incorrect value for Boundary : %s, not compatible with model %s\n (as defined in simul.c)", boundary_str, MOD_NAME);
      exit(-2);
#endif
    } else {
      ErrPrintf("Incorrect value for Boundary : %s\n", boundary_str);
      exit(-2);
    }
    pbc_mode = (boundary == BC_PERIODIC);
  }
  if (ava_mode_str) {
    if (!strcmp(ava_mode_str, "NONE")) {
      ava_mode = AVA_NONE;
    } else if (!strcmp(ava_mode_str, "SYNC")) {
      ava_mode = AVA_SYNC;
    } else if (!strcmp(ava_mode_str, "TRANS")) {
      ava_mode = AVA_TRANS;
    } else if (!strcmp(ava_mode_str, "PROPAG")) {
      ava_mode = AVA_PROPAG;
    } else {
      ErrPrintf("Incorrect value for Ava_mode : %s\n", ava_mode_str);
      exit(-2);
    }
  }
#ifdef ROTATIONS
  if (rot_mode_str) {
    if (!strcmp(rot_mode_str, "DISK") || !strcmp(rot_mode_str, "TABLE")) {
      rot_mode = ROT_MODE_DISK;
    } else if (!strcmp(rot_mode_str, "OVERLAP") || !strcmp(rot_mode_str, "CYCLE")) {
      rot_mode = ROT_MODE_OVERLAP;
    } else if (!strcmp(rot_mode_str, "FULL") || !strcmp(rot_mode_str, "SQUARE")) {
      rot_mode = ROT_MODE_FULL;
    } else {
      ErrPrintf("Incorrect value for Rot_mode : %s\n", rot_mode_str);
      exit(-2);
    }
  }
  if ((rot_mode == ROT_MODE_FULL) && (boundary == BC_REINJECTION)) {
    ErrPrintf("Incorrect value for Rot_mode: FULL, not compatible with REINJECTION boundary conditions\n");
    exit(-2);
  }
#endif
}
void init_simul() {

  simul_parse();
  ava_trans = (ava_mode == AVA_TRANS) || (ava_mode == AVA_PROPAG);
  ava_norm = (ava_angle != 0);
  simul_dump_flag = (dump_delay_png>0) || (dump_delay_csp>0);
  output_headers(); //write headers for all output files
}

// computation of the probability distribution
inline void compute_prob_dist() {
  int32_t i;
  double PI_ij_time, ccc;
  TransitionDb *ptr;
  TransitionCel *ptr_cel;


  // doublet transitions
  ptr = t_trans;
  for (i = 0; i < nb_trans_db; i++, ptr++) {
    //regulation of a transition
    if (ptr->regul) {
      ptr->regul(&ptr->intensite);
    }
    PI_ij[i] = (double) Ndb[ptr->depart] * ptr->intensite;
#ifdef CGV
#endif
    //debug
  }
  // cell transitions
  ptr_cel = t_trans_cel;
  for (; i < nb_trans; i++, ptr_cel++) {
    if (ptr_cel->regul) {
      ptr_cel->regul(&ptr_cel->intensite);  //regulation de la transition
    }
    PI_ij[i] = (double) Ncel[ptr_cel->depart] * ptr_cel->intensite;
  }
  // calcul PII : overall probability for any transition to occur in a unit of time
  // calcul PII_time : the same only for the transitions with time evolution (and time correction if transitions are blocked)
  for (i = 0, PII = 0.0, ccc = 0.0, PII_time = 0.0; i < nb_trans; i++) {
    PII  += PI_ij[i];
    if (i < nb_trans_db) {
      if (t_trans[i].time_mode == TIME_EVOL) {
        PI_ij_time = PI_ij[i];
      } else if (t_trans[i].time_mode == TIME_CORR) {
#ifdef CGV
        PI_ij_time = PI_ij[i] * cgv_coef;
#else
        ErrPrintf("calcule_prob: no coefficient for time correction\n");
        exit(-1);
#endif
      } else {
        PI_ij_time = 0;
      }
    } else {
      PI_ij_time = PI_ij[i];
    }
    cum[i] = PI_ij[i] + ccc;
    ccc    = cum[i];
    PII_time += PI_ij_time;
  }
}

void simul_time() {
  double aleat, t;
  //computation of csp and real time evolution
  //using Poisson formula
  if (PII_time) {
    aleat = drand48();
    if (!aleat) {
      ErrPrintf("ERROR: simul_time - aleat = 0 !\n");
      exit(-1);
    }
    t = -log(aleat) / PII_time;
    csp_time += t;
#ifdef TIME_SCALE
    real_time += t * time_scale;
#endif
  } else {
    LogPrintf("simul_time: PII_time = 0\n");
  }
  simul_stop();
}

inline int32_t simul_trans() {
  int32_t i, elt, tr;
  double ccc, aleat;
  int32_t dir;
  static int32_t cpt_blk = 0;
  int32_t max_blk = HLD;
  char flag_evol_time = 0;
  char flag_trans_done = 0;
  while ((!flag_trans_done) && (cpt_blk < max_blk)) {
    //we choose randomly which transition will occur
    tr = 0;
    aleat = drand48();
    ccc = PII * aleat;
    while (ccc >= cum[tr]) {
      tr++;
    }
    //debug
    if (tr > nb_trans) {
      ErrPrintf("ERROR: array overflow in t_trans (%d)\n", tr);
      exit(-1);
    }
    if (tr < nb_trans_db) {
      /// DOUBLET TRANSITION
      //we choose randomly which doublet will undergo the transition
      int32_t db_depart = t_trans[tr].depart;
      elt = (int) floor(drand48() * Ndb[db_depart]);
      int32_t ix = db_pos[db_depart][elt];
      //direction (class) of the transition
      if (t_trans[tr].classe == HORIZONTAL) {
#ifdef REFDB_PTR
        dir = (RefDB[ix][EST] == &db_pos[db_depart][elt]) ? EST : SUD;
#else
        dir = ((RefDB_Type[ix][EST] == db_depart) && (RefDB_Ind[ix][EST] == elt)) ? EST : SUD;
#endif
      } else {
        dir = direction[t_trans[tr].classe];
      }
      //control of the transition
      if (t_trans[tr].nb_chk && !simul_check(tr, ix, dir)) {
        //blocked transition
        cpt_blk++;
        continue;
      }
      do_trans_db(tr, ix, dir);
      flag_trans_done = 1;
      iter++;
      if (iter > 1000000000L) {
        iter -= 1000000000L;
        md_iter++;
      }
#ifdef LIENS_TRANSITIONS
      //linked transitions
      if (t_trans[tr].lien) {
        simul_lien_trans(tr, ix, dir);
      }
#endif
      flag_evol_time = (t_trans[tr].time_mode != TIME_NO_EVOL);
      if (flag_evol_time) {
        cpt_blk = 0;  /// reset counter when last transition makes time evolve
      }
    } else {
      /// CELL TRANSITION
      static int32_t imax = 0;
      int32_t ix;

      tr -= nb_trans_db;
      i = 0;
      //we look for a cell of the desired type
      do {
        ix = (int) floor(drand48() * HLD);
        i++;
      } while (TE[ix].celltype != t_trans_cel[tr].depart);
      if (i > imax) {
        imax = i;
        if (imax > 100) {
          LogPrintf("WARNING: cell transition was hard to achieve (%d trials)\n", imax);
        }
      }
      do_trans_cel(tr, ix);
      iter++;
      flag_trans_done = 1;
      flag_evol_time = 1;
      cpt_blk = 0;
    }
  }


  if (!flag_trans_done) {
    LogPrintf("simul_trans: too many blocked transitions (%d)\n", cpt_blk);
    cpt_blk = 0;
  }
  //time evolution
  if (flag_evol_time) {
    simul_time();
  }
  return flag_trans_done;
}

#ifdef TIME_SCALE
extern float lambda_C;
extern float lambda_T;
//#ifdef CGV
void simul_time_scale(float wind_value) {
  FILE *fp;
  static int32_t nb_rows;
  static float *tau1_values; //=grdvc_min
  static float *erosion_rates;
  static float *qsat_ratios;
  static float kappa; //Karman constant
  static float z_0; //roughness length
  static float z; //height of observation (velovity gradient)
  static float rho_air; //density of air
  static float rho_sand; //density of sand
  static float g; //gravity constant
  static float ut; //threshold velocity (erosion)
  static float d; //diameter of grain
  static float l_0; //length scale
  static char start = 1;
  int32_t i;
  if (!phys_prop_filename) {
    return;
  }
  if (start) {
    //read physical parameters file
    LogPrintf("simul_time_scale: read physical parameters file %s\n", phys_prop_filename);
    if ((fp = fopen(phys_prop_filename, "r")) == NULL) {
      ErrPrintf("ERROR: cannot open file %s\n", phys_prop_filename);
      exit(-4);
    }
    if (fscanf(fp, "kappa=%f\n", &kappa) != 1) {
      ErrPrintf("ERROR: cannot read parameter kappa in %s\n", phys_prop_filename);
      exit(-4);
    }
    LogPrintf("kappa=%f\n", kappa);
    if (fscanf(fp, "z_0=%f\n", &z_0) != 1) {
      ErrPrintf("ERROR: cannot read parameter z_0 in %s\n", phys_prop_filename);
      exit(-4);
    }
    LogPrintf("z_0=%f\n", z_0);
    if (fscanf(fp, "z=%f\n", &z) != 1) {
      ErrPrintf("ERROR: cannot read parameter z in %s\n", phys_prop_filename);
      exit(-4);
    }
    LogPrintf("z=%f\n", z);
    if (fscanf(fp, "rhoair=%f\n", &rho_air) != 1) {
      ErrPrintf("ERROR: cannot read parameter rho_air in %s\n", phys_prop_filename);
      exit(-4);
    }
    LogPrintf("rhoair=%f\n", rho_air);
    if (fscanf(fp, "rhosand=%f\n", &rho_sand) != 1) {
      ErrPrintf("ERROR: cannot read parameter rho_sand in %s\n", phys_prop_filename);
      exit(-4);
    }
    LogPrintf("rhosand=%f\n", rho_sand);
    if (fscanf(fp, "g=%f\n", &g) != 1) {
      ErrPrintf("ERROR: cannot read parameter g in %s\n", phys_prop_filename);
      exit(-4);
    }
    LogPrintf("g=%f\n", g);

    if (fscanf(fp, "d=%f\n", &d) != 1) {
      ErrPrintf("ERROR: cannot read parameter d in %s\n", phys_prop_filename);
      exit(-4);
    }
    LogPrintf("d=%f\n", d);

    fclose(fp);
    ut = 0.1 * sqrt((rho_sand / rho_air) * g * d);
    LogPrintf("ut=%f\n", ut);
    l_0 = 5 * (rho_sand / rho_air) * d / 4.;
    LogPrintf("l_0=%f\n", l_0);
    //read qsat data file
    LogPrintf("simul_time_scale: read qsat data file %s\n", qsat_data_filename);
    if ((fp = fopen(qsat_data_filename, "r")) == NULL) {
      ErrPrintf("ERROR: cannot open file %s\n", qsat_data_filename);
      exit(-4);
    }
    //number of points
    if (fscanf(fp, "%d", &nb_rows) != 1) {
      ErrPrintf("ERROR: cannot read the number of points in %s\n", qsat_data_filename);
      exit(-4);
    };
    LogPrintf("number of rows: %d\n", nb_rows);
    //read qsat data
    AllocMemory(tau1_values, float, nb_rows);
    AllocMemory(erosion_rates, float, nb_rows);
    AllocMemory(qsat_ratios, float, nb_rows);
    for (i = 0; i < nb_rows; i++) {
      if (fscanf(fp, "%f %f %f", &tau1_values[i], &erosion_rates[i], &qsat_ratios[i]) != 3) {
        ErrPrintf("ERROR: cannot read qsat data in %s (row %d)\n", qsat_data_filename, i);
        exit(-4);
      };
      //debug
    }
    fclose(fp);
  }
  //computation of the time_scale and tau1
  float ustar; //friction velovity
  float Qsat_ratio; //ratio of sand fluxes
  //threshold of the shear stress
  float qs; //sand flux in the nature
  float erosion_rate; //erosion rate
  float transport_rate; //transport rate
  float deposition_rate; //deposition rate
  float qsat; //sand flux in the model
  //time scale
  if (wind_value >= 0) {
    ustar = wind_value * kappa / log(z / z_0);
    Qsat_ratio = 1 - (ut / ustar) * (ut / ustar);
    for (i = 0; (i < nb_rows) && (Qsat_ratio < qsat_ratios[i]); i++);
    if (i < nb_rows) {
      //here we use PDF.data to get tau1 according to the value Qsat_ratio[i]
      //linear interpolation
      float t = (Qsat_ratio - qsat_ratios[i]) / (qsat_ratios[i - 1] - qsat_ratios[i]);
      grdvc_min = tau1_values[i - 1] * t + tau1_values[i] * (1 - t); //here we use PDF.data to get tau1 according to the value Qsat_ratio[i]
      erosion_rate = erosion_rates[i];
    } else {
      grdvc_min = tau1_values[nb_rows - 1];
      erosion_rate = erosion_rates[nb_rows - 1];
    }
    if (grdvc_min == 0) {
      grdvc_min = 1.0;  //value 0 is not allowed
    }
    grdvc_max = grdvc_min + 1000.0;
    LogPrintf("grdvc_min = %f\n", grdvc_min);
  } else {
    //if no wind data, then we compute ustar and erosion rate from tau1
    if (grdvc_min == 0) {
      grdvc_min = 1.0;  //value 0 is not allowed
    }
    for (i = 0; (i < nb_rows) && (grdvc_min > tau1_values[i]); i++);
    Qsat_ratio = qsat_ratios[i];
    erosion_rate = erosion_rates[i];
    if (Qsat_ratio < 1) {
      ustar = ut / sqrtf(1 - Qsat_ratio);
    } else {
      ustar = 0; //no erosion
    }
  }
  LogPrintf("ustar=%f\n", ustar);
  if (ustar > ut) {
    qs = 22.0 * (rho_air / rho_sand) * sqrtf(d / g) * (ustar * ustar - ut * ut); //qsat in the natural
    transport_rate = lambda_T;
    deposition_rate = lambda_C;
    qsat = transport_rate * erosion_rate / deposition_rate; //in model
    time_scale = qsat * (l_0 * l_0) / qs; //t_0
  } else {
    //no erosion in this case

    time_scale = 1e6; //not possible to set the time scale accurately
    static char first = 1;
    if (first) {
      ErrPrintf("Warning: no time scale for wind value %f\n", wind_value);
      first = 0;
    }
  }
  LogPrintf("time_scale = %g\n", time_scale);
  start = 0;
}
#endif
void simul_correl(int32_t tr1, int32_t ix, int32_t ix2, int32_t dir) {
  static double cum2[MAX_LIENS];
  static int32_t ln[MAX_LIENS], ind[MAX_LIENS];
  double ccc, aleat;
  int32_t i, db, db2, n;
  int32_t trans2_db_depart;
  //determination du doublet de depart
  db = type_doublet(ix, dir);
  db2 = type_doublet(ix2, dir);
  //construction du tableau des transitions correlees
  for (i = 0, ccc = 0.0, n = 0; i < nb_liens; i++) {
    if (t_lien[i].trans1 == tr1) {
      trans2_db_depart = t_trans[t_lien[i].trans2].depart;
      if (t_lien[i].cel == 1) {
        if (trans2_db_depart == db) {
          ccc += t_lien[i].intensite;
          cum2[n] = ccc;
          ln[n] = i;
          ind[n] = ix;
          n++;
        }
      } else {
        if (trans2_db_depart == db2) {
          ccc += t_lien[i].intensite;
          cum2[n] = ccc;
          ln[n] = i;
          ind[n] = ix2;
          n++;
        }
      }
    }
  }
  if (n > 0) { //test tableau non-vide ?
    //choix d'une transition selon les probabilites cumulees
    if (ccc > 1.0) {
      ErrPrintf("ERROR : cumulated probabilites = %f > 1.0\n", ccc);
      exit(-1);
    }
    aleat = drand48();
    i = 0;
    while ((i < n) && (aleat > cum2[i])) {
      i++;
    }
    if (i < n) {
      int32_t tr2 = t_lien[ln[i]].trans2;
      do_trans_db(tr2, ind[i], dir);
      iter++;
#ifdef INFO_TRANS
      cpt_lien[ln[i]]++;
#endif
      if (t_trans[tr2].lien) {
        simul_lien_trans(tr2, ind[i], dir);
      }
    }
  }
}

void simul_lien_trans(int32_t tr, int32_t ix, int32_t dir) {
  int32_t ix2;
  int32_t ix_west, ix2_west;
  int32_t ix_up, ix2_up;
  int32_t ix_north, ix2_north;
  if (t_trans[tr].classe == VERTICAL) {
    // PREMIERE TRANSITION VERTICALE
    ix2 = get_cell_down(ix);
    // liens horizontaux, doublet est
    simul_correl(tr, ix, ix2, EST);
    // liens horizontaux, doublet ouest
    ix_west = get_cell_west(ix);
    ix2_west = get_cell_down(ix_west);
    simul_correl(tr, ix_west, ix2_west, EST);
    // liens horizontaux, doublet sud
    simul_correl(tr, ix, ix2, SUD);
    // liens horizontaux, doublet nord
    ix_north = get_cell_north(ix);
    ix2_north = get_cell_down(ix_north);
    simul_correl(tr, ix_north, ix2_north, SUD);
  } else {
    // PREMIERE TRANSITION HORIZONTALE
    ix2 = get_cell_dir(ix, dir);
    int32_t flag_eo = (dir == EST/*EST_OUEST*/);
    // liens verticaux, doublet bas
    simul_correl(tr, ix, ix2, BAS);
    // liens verticaux, doublet haut
    ix_up = get_cell_up(ix);
    ix2_up = get_cell_up(ix2);
    simul_correl(tr, ix_up, ix2_up, BAS);
    if (!flag_eo) {
      // liens horizontaux, doublet est
      simul_correl(tr, ix, ix2, EST);
      // liens horizontaux, doublet ouest
      ix_west = get_cell_west(ix);
      ix2_west = get_cell_west(ix2);
      simul_correl(tr, ix_west, ix2_west, EST);
    } else {
      // liens horizontaux, doublet sud
      simul_correl(tr, ix, ix2, SUD);
      // liens horizontaux, doublet nord
      ix_north = get_cell_north(ix);
      ix2_north = get_cell_north(ix2);
      simul_correl(tr, ix_north, ix2_north, SUD);
    }
  }
}

int32_t simul_check(int32_t tr, int32_t ix, int32_t dir) {
  int32_t i;
  int32_t chk = 1;
  DataCheck *pdc;
  pdc = t_trans[tr].checks;
  for (i = 0; (i < t_trans[tr].nb_chk) && chk; i++, pdc++) {
    pdc->dir = dir;
    /// callback function on first cell
    if (pdc->cel == 1) { // test sur la premiere cellule
      chk = pdc->func(ix, (void*)pdc);
    } else { /// callback function on second cell
      int32_t ix2 = get_cell_dir(ix, dir);
      chk = pdc->func(ix2, (void*)pdc);
    }
    /// inversion flag
    if (pdc->inv) {
      chk = !chk;
    }
  }
  if (!chk) {
    cpt_trans_blocked[tr]++;
  }
  return chk;
}

#ifdef ROTATIONS
extern int32_t Ncel[];           // nombre de cellules par type
extern float coef_injection; //regulation de l'injection de grain
void simul_rot() {
  static char start = 1;
  static float time_threshold = 0.0;
  static int32_t cpt_rot = 0;
  static int16_t lg_seq = 5; //2
  static int16_t seq_rot[] = {0, 1, 2, 3, 4};
  static int16_t seq_del[] = {1, 1, 1, 1, 1};
  static int32_t total_seq_rot = 0;
  static int32_t seq_flag = 0; //flag for using sequences of angles/delays/...
  int32_t nseq;
  if (!rot_delay) {
    return;
  }
  if (start) {
    if (!seq_flag) {
      time_threshold = rot_delay * ceil(csp_time / rot_delay);
      if (time_threshold == csp_time) {
        time_threshold += rot_delay;
      }
    }
    LogPrintf("rotation mode = %s\n", (rot_mode == ROT_MODE_DISK) ? "rotating table" : (rot_mode == ROT_MODE_FULL) ? "full space" : "overlap");
    LogPrintf("rotation delay = %f\n", rot_delay);
    if (seq_flag) {
      LogPrintf("period od the sequence: %d\n", lg_seq);
    }
#ifdef LGCA
    if (use_lgca) {
      LogPrintf("number of lattice-gas stabilization cycles after rotation = %d\n", rot_mvt_nb_cycles);
      if (rot_mvt_nb_cycles < VSTEP_TIME) {
        ErrPrintf("WARNING: stabilization after rotation is too int16_t ( Rot_ncycl = %d < VSTEP_TIME = %d )\n", rot_mvt_nb_cycles, VSTEP_TIME);
      }
    }
#endif
    start = 0;
  }
  if (csp_time >= time_threshold) {
    //translation
    if (seq_flag) {
      //sequences de rotations
      //random direction in sequence
      nseq = cpt_rot % lg_seq;
      LogPrintf("nseq=%d   seq_rot=%d   seq_del=%d\n", nseq, seq_rot[nseq], seq_del[nseq]);
      rotation(seq_rot[nseq]*rot_angle, rot_mode, 0);
      total_seq_rot += seq_rot[nseq];
    } else {
      //rotation d'angle rot_angle


      if (1) { //(rot_angle)
        rotation(rot_angle, rot_mode, 0);
      } else {
        rotation(360 * drand48(), rot_mode, 0);
      }
    }
    cpt_rot++;
    if (seq_flag) {
      time_threshold += rot_delay * seq_del[nseq];
    } else {
      time_threshold += rot_delay;
    }
    LogPrintf("cpt_rot = %d\n", cpt_rot);
#if defined(GR) || defined(GRJ)
    LogPrintf("nombre de grains : %d\n", Ncel[GR] + Ncel[GRJ]);
#ifdef AVALANCHES
    if (ava_mode == AVA_SYNC) {
      avalanches(GR, ava_h_lim, ava_nb_cel_max, ALTI_MODE_BAS);
    }
#endif
#endif
#ifdef LGCA
    lgca_stabilize_flag = 1;
#endif
  }
}
#endif
#ifdef CENTERING_AUTO
void simul_centering() {
  static char start = 1;
  static float time_threshold = 0.0;
  static int32_t cx0, cz0;
  if (!centering_delay) {
    return;
  }
  if (start) {
    time_threshold = centering_delay * ceil(csp_time / centering_delay);
    cx0 = roundf(L / 2.0 - 0.5);
    cz0 = roundf(D / 2.0 - 0.5);
    start = 0;
  }
  if (csp_time >= time_threshold) {
    Vec3 center = compute_mass_center(ALTI);
    int32_t trx = -roundf(center.x - cx0);
    int32_t trz = -roundf(center.z - cz0);
    if (trx || trz) {
      translation(trx, trz);
#ifdef LGCA
      lgca_stabilize_flag = 1; /// force flow stabilization
#endif
    }
    time_threshold += centering_delay;
  }
}
#endif
#ifdef AVALANCHES
//synchronized avalanches (SYNC mode)
void simul_ava_sync() {
  static float time_threshold = 0.0;
  static uint64_t inter_iter_ava = 20000;  //nb min de transitions entre 2 avalanches
  static float ava_delay_sync = 0.0;
  static char start = 1;
  if (start) {
#ifdef DUMP_SIGMA
    inter_iter_ava = 1;
#endif
    ava_delay_sync = ava_delay;
    time_threshold = ava_delay_sync * ceil(csp_time / ava_delay_sync);
    LogPrintf("seuil temps avalanches = %f\n", time_threshold);
    LogPrintf("delai avalanches = %f\n", ava_delay_sync);
    LogPrintf("nb min iterations avant avalanches = %" PRIu64 "\n", inter_iter_ava);
    LogPrintf("ava_h_lim = %d\n", ava_h_lim);
    LogPrintf("ava_nb_cel_max = %d\n", ava_nb_cel_max);
    start = 0;
  }
  if ((!csp_time) || (csp_time >= time_threshold)) {
    avalanches(ALTI, ava_h_lim, ava_nb_cel_max, ALTI_MODE_BAS);
    if (ava_angle) {
      avalanches_norm(ALTI, ava_nb_cel_max, ALTI_MODE_BAS);
    }
#ifdef DUMP_SIGMA
    dump_sigma_alti();
#endif
#ifdef DUMP_AUTOCORREL
    dump_autocorrel();
#endif
    time_threshold = csp_time + ava_delay_sync;
  }
}

// dynamical regulation of the rate and angles of avalanches
void simul_ava_dynamics() {
  static float time_threshold = 0.0;
  static char start = 1;
  static char mode_stable = 1;
  if (ava_mode == AVA_PROPAG) {
    if (start) {
      time_threshold = 1 + ava_delay * ceil(csp_time / ava_delay);
      ava_angle_stable = ava_angle;
      ava_angle_stable_col = ava_angle_col;
      lambda_A_stable = lambda_A;
      //we set the unstable values for the transitions, if not null
      if (ava_angle_unstable) {
        ava_angle = ava_angle_unstable;
        mode_stable = 0;
      }
      if (lambda_A_unstable) {
        lambda_A = lambda_A_unstable;
        mode_stable = 0;
      }
      start = 0;
    }
    return;
  }
  if (start) {
    time_threshold = ava_duration * ceil(csp_time / ava_duration);
    ava_angle_stable = ava_angle;
    lambda_A_stable = lambda_A;
    ava_angle_stable_col = ava_angle_col;
    LogPrintf("simul_ava_dynamics:\n");
    LogPrintf("duration of avalanches = %f\n", ava_duration);
    LogPrintf("delay between avalanches = %f\n", ava_delay);
    LogPrintf("ava_angle = %f (stable)| %f (unstable)\n", ava_angle_stable, ava_angle_unstable);
    start = 0;
  }
  if ((!csp_time) || (csp_time >= time_threshold)) {
    if (mode_stable) {
      mode_stable = 0;
      ava_angle = ava_angle_unstable;
      lambda_A = lambda_A_unstable;
      time_threshold = csp_time + ava_delay;
    } else {
      mode_stable = 1;
      ava_angle = ava_angle_stable;
      lambda_A = lambda_A_stable;
      time_threshold = csp_time + ava_duration;
    }
  }
}
#endif //AVALANCHES
#ifdef CGV
//erosion rate
extern float grdvc_max; //seuil max pour le controle de l'erosion en fonction du gradient de vitesse
extern float grdvc_min; //seuil min pour le controle de l'erosion en fonction du gradient de vitesse
void simul_wind_variability() {
  static float delai_vent_max = 200.0;
  static float time_threshold = 0.0;

  static char start = 1;
  if (start) {
    time_threshold = delai_vent_max * ceil(csp_time / delai_vent_max);
    LogPrintf("parametres d'intensite du vent au debut : grdvc_min = %f   grdvc_max = %f\n", grdvc_min, grdvc_max);
    LogPrintf("delai max de variation du vent = %f\n", delai_vent_max);
    start = 0;
  }
  if (csp_time >= time_threshold) {
    grdvc_min = drand48() * 500.0;
    grdvc_max = grdvc_min + 1000.0;
    float delai_vent_alea = delai_vent_max * drand48();
    LogPrintf("variation d'intensite du vent : grdvc_min = %f   grdvc_max = %f   csp_time = %f\n", grdvc_min, grdvc_max, csp_time);
    LogPrintf("delai de variation du vent = %f\n", delai_vent_alea);
    time_threshold += delai_vent_alea;
  }
}
#endif //CGV
#ifdef WIND_DATA
// load wind data from an ascii file
// format:
// nb_wind
// delay angle value
// ...
void simul_wind_data() {
  static float *wind_delays = NULL;
  static float *wind_angles = NULL;
  static float *wind_values = NULL;
  static int32_t nb_wind = 0;
  static int32_t time_unit = 0;
  static int32_t cpt_wind = 0;
  static float time_threshold = 0.0;
  static float angle_dir = 0.0;
  float angle = 0.0;
  float seq_duration = 0.0;
  char time_unit_str[128];
  static char start = 1;
  if (start) {
    LogPrintf("simul_wind_data: read file %s\n", wind_data_filename);
    LogPrintf("rotation mode = %s\n", (rot_mode == ROT_MODE_DISK) ? "rotating table" : (rot_mode == ROT_MODE_FULL) ? "full space" : "overlap");
    //open wind data file
    FILE *fp;
    if ((fp = fopen(wind_data_filename, "r")) == NULL) {
      ErrPrintf("ERROR: cannot open file %s\n", wind_data_filename);
      exit(-4);
    }
    //number of wind regimes
    if (fscanf(fp, "%d,%s", &nb_wind, time_unit_str) != 2) {
      ErrPrintf("ERROR: cannot read the number of wind regimes in %s\n", wind_data_filename);
      exit(-4);
    };
    LogPrintf("number of winds: %d\n", nb_wind);
    if (!strcmp(time_unit_str, "t0")) {
      time_unit = UNIT_T0;
      LogPrintf("time unit is t0\n");
    } else if (!strcmp(time_unit_str, "rt")) {
      time_unit = UNIT_RT;
      LogPrintf("time unit is physical\n");
#ifndef TIME_SCALE
      ErrPrintf("ERROR: no time scale, not possible to compute wind data from file %s\n", wind_data_filename);
      exit(-4);
#endif
    } else {
      ErrPrintf("ERROR: cannot read time unit in %s\n", wind_data_filename);
      exit(-4);
    }
    //read wind data
    int32_t i;
    AllocMemory(wind_delays, float, nb_wind);
    AllocMemory(wind_angles, float, nb_wind);
    AllocMemory(wind_values, float, nb_wind);
    for (i = 0; i < nb_wind; i++) {
      if (fscanf(fp, "%f %f %f", &wind_delays[i], &wind_angles[i], &wind_values[i]) != 3) {
        ErrPrintf("ERROR: cannot read wind data in %s\n", wind_data_filename);
        exit(-4);
      };
      //duration of the wind sequence in t0 unit
#ifdef TIME_SCALE
      if (time_unit == UNIT_RT) {
        simul_time_scale(wind_values[i]);
        seq_duration += wind_delays[i] / time_scale;
      }
#endif
      if (time_unit == UNIT_T0) {
        seq_duration += wind_delays[i];
      }
      if (i < 10) {
        LogPrintf("%d\t%f\t%f\t%f\n", i, wind_delays[i], wind_angles[i], wind_values[i]);
      }
    }
    if (i > 10) {
      LogPrintf("...\n");
    }
    LogPrintf("duration of the wind sequence = %f t0\n", seq_duration);
    fclose(fp);
    //process wind data
#ifdef TIME_SCALE
    if (time_unit == UNIT_RT) {
      //real time
      //initial wind and time_threshold
      time_threshold = seq_duration * floor(csp_time / seq_duration);
      while (time_threshold <= csp_time) {
        simul_time_scale(wind_values[cpt_wind]);
        time_threshold += wind_delays[cpt_wind++] / time_scale;
      }
      cpt_wind--;
      LogPrintf("cpt_wind = %d   time_threshold = %f\n", cpt_wind, time_threshold);
      //first wind orientation (north = 0 degree)
      angle = wind_angles[cpt_wind] + 90;
    }
#endif
    if (time_unit == UNIT_T0) {
      //virtual time unit t0
      //initial wind and time_threshold
      time_threshold = seq_duration * floor(csp_time / seq_duration);
      while (time_threshold <= csp_time) {
        time_threshold += wind_delays[cpt_wind++];
      }
      cpt_wind--;
      LogPrintf("cpt_wind = %d   time_threshold = %f\n", cpt_wind, time_threshold);
      //first wind orientation (west = 0 degree)
      angle = wind_angles[cpt_wind];
#ifdef CGV
      //first value of min shear stress
      if (!param_is_set("Tau_min")) {
        grdvc_min = wind_values[cpt_wind];
      }
#ifdef TIME_SCALE
      //time scale
      simul_time_scale(-1);
#endif
      LogPrintf("grdvc_min = %f\n", grdvc_min);
#endif
    }
    LogPrintf("simul_wind_data: direction = %f, rotation = %f, time = %e\n", wind_angles[cpt_wind], angle, csp_time);
    rotation(-angle, rot_mode, 0);
    angle_dir = wind_angles[cpt_wind]; //first wind direction
    start = 0;
  }
  if (csp_time >= time_threshold) {
    //we change the wind counter
    cpt_wind++;
    if (cpt_wind == nb_wind) {
      cpt_wind = 0;
    }
    //wind rotation
    angle = wind_angles[cpt_wind] - angle_dir;
    LogPrintf("simul_wind_data: direction = %f, rotation = %f, time = %e\n", wind_angles[cpt_wind], angle, csp_time);
    LogPrintf("cpt_wind = %d\n", cpt_wind);
    rotation(-angle, rot_mode, 0);
    //new wind direction
    angle_dir = wind_angles[cpt_wind];
#ifdef TIME_SCALE
    if (time_unit == UNIT_RT) {
      //time scale
      simul_time_scale(wind_values[cpt_wind]);
      //real time
      time_threshold += wind_delays[cpt_wind] / time_scale;
      LogPrintf("real_time = %e\n", real_time);
    }
#endif
    if (time_unit == UNIT_T0) {
      //virtual time unit t0
      time_threshold += wind_delays[cpt_wind];
#ifdef CGV
      //next value of min shear stress
      if (!param_is_set("Tau_min")) {
        grdvc_min = wind_values[cpt_wind];
      }
#ifdef TIME_SCALE
      //time scale
      simul_time_scale(-1);
      LogPrintf("real_time = %e\n", real_time);
#endif
      LogPrintf("grdvc_min = %f\n", grdvc_min);
#endif
    }
  }
#ifdef LGCA
  if (angle) {
    lgca_stabilize_flag = 1;
  }
#endif
}
#endif //WIND_DATA

#ifdef LGCA
extern uint8_t opt_vss;
extern float *grdv;
void simul_profile_lgca(int32_t lgca_flag) {
  static double lgca_elapsed_time = 0.0;
  static double other_elapsed_time = 0.0;
  static double t0 = 0.0;
  static double t1 = 0.0;
  static int32_t cpt = 0;
  if (cpt > 1000) {
    return;
  }
  if (lgca_flag == 0) {
    elapsed(&t0);
    if (t1) {
      other_elapsed_time += (t0 - t1);
    }
  } else if (lgca_flag == 1) {
    elapsed(&t1);
    lgca_elapsed_time += (t1 - t0);
    if (cpt && (cpt % 100 == 0)) {
      LogPrintf("csp_time = %f\n", csp_time);
      if (mode_pat) {
        LogPrintf("time spent in lattice gas : %.2f sec (%.1f %%)\n", lgca_elapsed_time, 100 * lgca_elapsed_time / (other_elapsed_time + lgca_elapsed_time));
        LogPrintf("time spent waiting for stochastic engine : %.2f sec (%.1f %%)\n", other_elapsed_time, 100 * other_elapsed_time / (other_elapsed_time + lgca_elapsed_time));
      } else {
        LogPrintf("time spent in stochastic engine : %.2f sec (%.1f %%)\n", other_elapsed_time, 100 * other_elapsed_time / (other_elapsed_time + lgca_elapsed_time));
        LogPrintf("time spent in lattice gas : %.2f sec (%.1f %%)\n", lgca_elapsed_time, 100 * lgca_elapsed_time / (other_elapsed_time + lgca_elapsed_time));
      }
    }
    cpt++;
  }
}
void flow_stabilization(int32_t nb_cyc) {
  int32_t i;
  char flag_interp = 0;
#ifdef CGV
  calcule_alti(ALTI, ALTI_MODE_BAS);
  calcule_normales();
#endif
  LogPrintf("flow stabilization : %d cycles\n", nb_cyc);
  lgca_ready = 0;
  for (i = 0; i < nb_cyc; i++) {
    //pour ralentir ...
    do_collisions();
    do_propagations();
    /// interpolation only during the last cycle (faster)
    flag_interp = (i == nb_cyc - 1);
    compute_vel(flag_interp);
    do_thread_sched();
  }
  lgca_stabilize_flag = 0;
  lgca_ready = 1;
  LogPrintf("end of flow stabilization : meanvel = %f   maxvel = %f\n", meanvel, maxvel);
}
void simul_lgca() {
  /// start profiling
  simul_profile_lgca(0);
  /// lgca stabilization
  if (lgca_stabilize_flag) {
    if (lgca_reset) {
      init_mvt();
    }
    flow_stabilization(rot_mvt_nb_cycles);
  }
  /// lgca cycle
  if (rot_map) {
    out_of_space_mvt(0);
  }
  do_collisions();
  do_propagations();
  compute_vel(1);
  /// stop profiling
  simul_profile_lgca(1);
}
#endif
#ifdef ALTI
extern uint8_t opt_ls;
void simul_norm_lum() {
  static float seuil_calcule_normales = 0.0;
  static float delai_calcule_normales = 1.0;
  if ((opt_ls) && (csp_time >= seuil_calcule_normales)) {
    //recalcul des normales pour le rendu avec eclairage
    calcule_normales();
    seuil_calcule_normales += delai_calcule_normales;
  }
}
#endif
void simul_stop() {
  if (stop_time && (csp_time >= stop_time)) {
    if (simul_dump_flag) {
      simul_dump();
    }
    LogPrintf("csp_time=%g\n", csp_time);
    LogPrintf("time to quit\n");
    sleep(1);
    exit(0);
  }
}

int32_t simul_csp() {
  int32_t ii;
  int32_t res;
  /// start time
  if (param_is_set("Time")) {
    csp_time = init_time;
  }
  LogPrintf("csp_time = %f\n", csp_time);
  /// stop time
  if (stop_delay_t0) {
    stop_time = csp_time + stop_delay_t0;
    LogPrintf("stop_time = %f\n", stop_time);
  }
#ifdef ALTI
  calcule_alti(ALTI, ALTI_MODE_BAS);
  simul_norm_lum();
#endif
  if (simul_dump_flag) {
    lock_csp(0); //a rotation can occur here
    simul_dump();
    unlock_csp(0);
  }
#ifdef WIND_DATA
  if (wind_data_filename) {
    lock_csp(0); //a rotation can occur here
    simul_wind_data();
    unlock_csp(0);
  }
#endif
#ifdef TIME_SCALE
  /// time scale computation
  if (!wind_data_filename) {
    simul_time_scale(-1);
  }
#endif
#ifdef LGCA
  double lgca_time_threshold = csp_time;
  if (use_lgca && !lgca_delay) {
#if defined(MODEL_DUN) || defined(MODEL_SNO)
    lgca_delay = 1.0 / NB_MVT_EO;
#else
    /// default lgca delay derived from transition rates
    double lambda_max = 0;
    for (i = 0; i < nb_trans_db; i++) {
      if (t_trans[i].intensite > lambda_max) {
        lambda_max = t_trans[i].intensite;
      }
    }
    lgca_delay = 1.0 / lambda_max;
#endif
  }
  LogPrintf("lgca_delay = %f\n", lgca_delay);
#endif //LGCA
  /// CSP ENGINE START
#ifdef DUMP_SIGNATURE
  if (opt_info) {
    dump_signature(0);
  }
#endif
#ifdef DUMP_SIGMA
  dump_sigma_alti();
#endif
#ifdef DUMP_AUTOCORREL
  dump_autocorrel();
#endif
  lock_csp(0);
#ifdef CENTERING_AUTO
  simul_centering();
#endif
#ifdef AVALANCHES
  if (ava_mode == AVA_SYNC) {
    simul_ava_sync();
  }
  if (ava_delay) {
    simul_ava_dynamics();
  }
#endif //AVALANCHES
#ifdef LGCA
  if (use_lgca) {
    /// stabilization of lattice gas
    if (!param_is_set("Init_ncycl")) {
      init_mvt_nb_cycles = 2 * CLEO;  //default value
    } else if (init_mvt_nb_cycles < VSTEP_TIME) {
      init_mvt_nb_cycles = VSTEP_TIME;  //min value
    }
    if (nom_fic_mvt) {
      init_mvt_nb_cycles = VSTEP_TIME;
    } else if (lgca_stabilize_flag) {
      init_mvt(); /// in case of an initial non-zero orientation
    }
    flow_stabilization(init_mvt_nb_cycles);
    if (!param_is_set("Rot_ncycl")) {
      rot_mvt_nb_cycles = lgca_reset && (!lgca_speedup) ? 2 * CLEO : 3 * VSTEP_TIME;
    }
#ifdef CGV
    if (opt_info) {
      dump_cgv();
    }
#endif
  }
#endif //LGCA
  ii = 0;
  ResetMemory(cpt_trans_blocked, int, MAX_TRANSITIONS_DB);
  /// MAIN LOOP
  while (!end_of_simul) {
    do_thread_sched();
    //pour ralentir ...
#ifdef LGCA
    /// LATTICE GAS
    if (use_lgca) {
      if ((csp_time >= lgca_time_threshold) || lgca_stabilize_flag) {
#ifdef PARALLEL_AUTOMATA
        if (mode_pat) {
          wait_lgca_thread();
        } else {
          simul_lgca();
        }
#else
        simul_lgca();
#endif
#ifdef ALTI
        calcule_normales();
#endif
#ifdef CGV
        calcule_grad_vel();
        compute_coef_cgv();
#endif
        lgca_time_threshold += lgca_delay;
        while (lgca_time_threshold < csp_time) {
          lgca_time_threshold += lgca_delay;
        }
      }
    } else {
      simul_norm_lum();
    }
#elif defined(ALTI)
    simul_norm_lum();
#endif
    /// STOCHASTIC ENGINE
    /// compute the probability distribution
    compute_prob_dist();
    if (PII == 0) {
      ///no active doublets
      if (!use_lgca) {
        /// stop simulation
        LogPrintf("no active doublets\n");
        end_of_simul = 1;
        break;
      }
#ifdef LGCA
      else {
        /// force time evolution
        LogPrintf("PII=0, csp_time=%g\n", csp_time);
        csp_time = lgca_time_threshold;
        simul_stop();
      }
#endif
    }

    if (PII > 0) {
      /// try to perform a stochastic transition
      res = simul_trans();
      if (!res) {
        /// no transition performed
#ifdef LGCA
        if (use_lgca) {
          LogPrintf("simul_csp: time evolution forced (%f)\n", csp_time);
          if (csp_time < lgca_time_threshold) {
            csp_time = lgca_time_threshold;
            simul_stop();
          } else {
            while (lgca_time_threshold < csp_time) {
              lgca_time_threshold += lgca_delay;
            }
          }
        }
#endif
      }
    }
    if (simul_dump_flag) {
      simul_dump();
    }
#ifdef CENTERING_AUTO
    simul_centering();
#endif
#ifdef ROTATIONS
    simul_rot();
#endif
#ifdef AVALANCHES
    if (ava_mode == AVA_SYNC) {
      simul_ava_sync();
    }
    if (ava_delay) {
      simul_ava_dynamics();
    }
#endif //AVALANCHES
#ifdef CGV
#endif
#ifdef WIND_DATA
    if (wind_data_filename) {
      simul_wind_data();
    }
#endif
    ii++;
#ifdef DUMP_SIGNATURE
    if (opt_info &&  !(ii & 0x000fffff)) { //toutes les 1048576 transitions
      /// save fingerprint32_t of the cellular space (hash function)
      dump_signature(ii);
    }
#endif
#ifdef TRACE3D_CEL
    //toutes les 8192 transitions
    //toutes les 256 transitions
    trace3d_cel();  //a chaque transition
#endif
#ifdef TRACE_FLUX
    extern float trace_flux_delay;
    if (trace_flux_delay) {
      trace_dump(0);
    }
#endif
#ifdef DUMP_SIGMA
#define SIGMA_DELAY 1.0 //delay in t0 unit
    //toutes les 4096 transitions
    //toutes les 256 transitions
    static double sigma_delay = SIGMA_DELAY;
    static double sigma_time_threshold = SIGMA_DELAY;
    if (csp_time >= sigma_time_threshold) {
      dump_sigma_alti();
      sigma_time_threshold += sigma_delay;
    }
#endif
#ifdef DUMP_AUTOCORREL
#define AUTOCOR_DELAY 10.0 //delay in t0 unit
    static double autocor_delay = AUTOCOR_DELAY;
    static double autocor_time_threshold = AUTOCOR_DELAY;
    if (csp_time >= autocor_time_threshold) {
      dump_autocorrel();
      autocor_time_threshold += autocor_delay;
    }
#endif
  }
  //END OF MAIN LOOP
  unlock_csp(0);
  LogPrintf("end of simulation\n");

  return (end_of_simul);
}

void simul_dump() {
  static float dump_delay = 0.0;
  static double time_threshold = 0.0;
  static double time_threshold_dpng = 0.0;
  static double time_threshold_dcsp = 0.0;
  static int32_t cpt_dump = 0;
  static char name[512];
  static char str[100];
  static double t0 = 0, t1 = 0;
  static char start = 1;
  static double epsilon = 0.1;
  if (start) {
    *str = 0;
    //set dump_delay = minimal non-zero delay
    if (dump_delay_png > 0) {
      time_threshold_dpng = dump_delay_png * ceil((csp_time - epsilon) / dump_delay_png); //time_threshold_dpng is mulitple of dump_delay_png
      dump_delay = dump_delay_png;
      LogPrintf("dump_delay_png = %f (t0)\n", dump_delay_png);
    }
    if (dump_delay_csp > 0) {
      time_threshold_dcsp = dump_delay_csp * ceil((csp_time - epsilon) / dump_delay_csp); //time_threshold_dcsp is mulitple of dump_delay_csp
      if (!dump_delay || (dump_delay_csp < dump_delay)) {
        dump_delay = dump_delay_csp;
      }
      LogPrintf("dump_delay_csp = %f (t0)\n", dump_delay_csp);
    }
    if ((dump_delay_png > 0) && (dump_delay_csp > 0)) {
      time_threshold = Min(time_threshold_dpng, time_threshold_dcsp);
    } else if (dump_delay_png > 0) {
      time_threshold = time_threshold_dpng;
    } else {
      time_threshold = time_threshold_dcsp;
    }
    LogPrintf("time_threshold = %f\n", time_threshold);
    LogPrintf("dump_delay = %f (t0)\n", dump_delay);
    start = 0;
  }
  if ((csp_time >= time_threshold) || end_of_simul) {
    //no more than 10 files per sec., for security
    elapsed(&t1);
    if (!t0 || (t1 - t0 >= 0.1)) {
#ifdef REORIENT_AUTO
      reorient_flag = 1;
      rotation(0, rot_mode, ROT_REORIENT_TEMP);
      //display wind direction
      vdir_mode = VDIR_WIND;
#endif
      if (dump_delay_png && (time_threshold == time_threshold_dpng)) {
        //dump PNG file
        sprintf(name, "%s%05d_t0%s.png", MOD_NAME, cpt_dump, str);
        dump_image(name, "png");
        while (time_threshold_dpng <= csp_time) {
          time_threshold_dpng += dump_delay_png;
        }
      }
      if (dump_delay_csp && (time_threshold == time_threshold_dcsp)) {
        // edited to allow dumping of ALTI files but not .csp files using the altionly flag
        if (!alti_only_flag) {
          dump_terre(DUMP_CSP, cpt_dump, UNIT_T0);
        }
#ifdef ALTI
        dump_surface("ALTI", cpt_dump, UNIT_T0);
#endif
#ifdef LGCA
        if (use_lgca && csphpp_flag) {
          dump_mvt(cpt_dump, UNIT_T0);
#ifdef CGV
          dump_grad_vel(cpt_dump, UNIT_T0);
#endif
        }
#endif
        while (time_threshold_dcsp <= csp_time) {
          time_threshold_dcsp += dump_delay_csp;
        }
      }
#ifdef REORIENT_AUTO
      rotation(0, rot_mode, ROT_REORIENT_UNDO);
      vdir_mode = VDIR_NONE;
      reorient_flag = 0;
#endif
      t0 = t1;
      //update cpt
      cpt_dump++;
      cpt_dump = cpt_dump % 100000; //no more than 100000 files, for security
    } else {
      if (time_threshold == time_threshold_dpng) {
        time_threshold_dpng += dump_delay_png;
      }
      if (time_threshold == time_threshold_dcsp) {
        time_threshold_dcsp += dump_delay_csp;
      }
    }
    //next time threshold
    while (time_threshold <= csp_time) {
      time_threshold += dump_delay;
    }
    if (dump_delay_png && (time_threshold_dpng < time_threshold)) {
      time_threshold = time_threshold_dpng;
    }
    if (dump_delay_csp && (time_threshold_dcsp < time_threshold)) {
      time_threshold = time_threshold_dcsp;
    }
  }
}

void dump_time()
{
  char current_output[256];
  static int32_t cpt = 0;
  static uint64_t md_iter_0 = 0, iter_0 = 0;
  static double csp_time_0 = 0;
  static double real_time_0 = 0;
  static int32_t col_iter_0 = 0;
  int64_t delta_iter, delta_md_iter;

  // First time this function is called
  if (!cpt){ 
    csp_time_0  = csp_time;
    strcat(current_output, "      nb trans.        delta trans.    time            delta time       ");
    if (use_lgca) strcat(current_output, "lgca cyc.    ");
#ifdef TIME_SCALE
    strcat(current_output, "real time     delta real time");
    real_time_0 = real_time;
#endif
    output_write("TIME", current_output);
  }

  // Every call
  delta_iter = iter - iter_0;
  delta_md_iter = md_iter - md_iter_0;
  if ((delta_md_iter > 0) && (delta_iter < 0)) {
    delta_iter += 1000000000;
    delta_md_iter--;
  }

  if (delta_md_iter){
    sprintf(current_output,"\n%04d: %04" PRIu64 "%09" PRId64" %03" PRId64 "%09" PRId64 "       %e    %e     ", cpt++, md_iter, iter, delta_md_iter, delta_iter, csp_time, csp_time - csp_time_0);
  }
  else {
    sprintf(current_output,"\n%04d: %04" PRIu64 "%09" PRId64"    %09" PRId64 "       %e    %e     ", cpt++, md_iter, iter, delta_iter, csp_time, csp_time - csp_time_0);
  }
  output_write("TIME", current_output);
#ifdef LGCA
  if (use_lgca){
    sprintf(current_output,"%09d    ", col_iter - col_iter_0);
    col_iter_0 = col_iter;
  }
#endif
#ifdef TIME_SCALE
  sprintf(current_output,"%e   %e", real_time, real_time - real_time_0);
  output_write("TIME", current_output);
#endif

  md_iter_0 = md_iter;
  iter_0 = iter;
  csp_time_0 = csp_time;
#ifdef TIME_SCALE
  real_time_0 = real_time;
#endif
}

void output_path(char *filename){
 // Replaces string "filename" with string "{output_directory}/filename.log"
 // KK
  
  char name[strlen(filename)+2];
  strcpy(name, filename);

  strcpy(filename, output_directory);
  strcat(filename, "/");
  strcat(filename, name);
  strcat(filename, ".log");
}

void output_path_noext(char *filename){
  // Replaces string "filename" with string "{output_directory}/filename"
  // without adding any extension.
  // KK
  char name[strlen(filename) + 2];
  strcpy(name, filename);

  strcpy(filename, output_directory);
  strcat(filename, "/");
  strcat(filename, name);
}

void output_write(char *output_filename, char *output_content){
  // Write output string output_content
  // to file output_directory/output_name.log
  // (This opens and closes the file every time it is called.
  //  Inefficient if called often for small blocks of text.
  //  Currently this happens in numerous functions.)
  FILE *fp;
  char path[256] = {'\0'};

  if (strlen(output_filename) + strlen(output_content) > 250){
    ErrPrintf("Unexpectedly long file/directory name");}

  strcat(path, output_filename);
  output_path(path);

  fp = fopen(path, "a");
  if (! fp ) {
	  ErrPrintf("ERROR: cannot open file: %s \n", path);
	  exit(-1);
  }
  fprintf(fp, "%s", output_content);
  fclose(fp);
}

void output_headers(){
  // List expected output files, and provide headers and metadata for each of them.
  // This file should be called once, before any output is written.
  // Output goes to a variety of files in directory output_directory/*.log
  // KK
  // TODO add *useful* and *complete* metadata to headers on output
  // Current headers are copied from previous output functions
  // TODO output_directory and output_filename will currently break in the face of typos
  
  void output_write(char *output_filename, char *output_content);

  // Check if output directory exists. If not, create it.
  // Uses #include <sys/types.h>, <sys/stat.h>, <unistd.h>
  struct stat st = {0};
  if (stat(output_directory, &st) == -1) mkdir(output_directory, 0777);
  LogPrintf("\n Trying to write output to directory : %s", output_directory);

  // Cell states (adds additional output in dump_cell)
  output_write("CELL",     "\n# CELL STATES\n");
  // CGV coefficients
  output_write("CGV_COEF", "        cgv_coef");
  // Doublets (additional output in dump_doublets and dump_dbl_info)
  if (opt_info) dump_doublets();
  // DENSITE.log
  output_write("DENSITE", "\tmvt cells \tfluid nodes \tdensity \ttrapped cells\n");
  // genesis
  // lgca to LGCA.log
  if (opt_info)  dump_collisions(); // information about fluid-particle collisions allowed in simulation
  // mvt_io
  // prob_cgv
  // rescal
  // time
  // trans
  if (opt_info) dump_transitions(); // information about transitions allowed in simulation
#ifdef INFO_TRANS
  dump_trans_info_header();
#endif
  // Velocities
  output_write("VEL", "      \t maxvel \t meanvel\n"); 
} 

