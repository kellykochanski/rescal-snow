/* ReSCAL - Predefined models
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


/// USER'S GUIDE
/// ------------
///
/// Basically a model is a set of doublet or cell transitions.
/// Here we describe all the functions that help define those transitions.
///
/// trans( direction, cel1_start, cel2_start, cel1_end, cel2_end, rate ):
///   defines a transition of oriented doublet [cel1_start, cel2_start] -> [cel1_end, cel2_end] occuring at the specified rate.
///   Possible values for directions are EST_OUEST, NORD_SUD, VERTICAL, HORIZONTAL, ISOTROPE.
///
/// trans_ref( ref,  directions,  cel1_start,  cel2_start,  cel1_end,  cel2_end,  rate ):
///   defines a transition of oriented doublet [cel1_start, cel2_start] -> [cel1_end, cel2_end] occuring at the specified rate,
///   with a unique reference number 'ref' (arbitrary positive integer) which identifies the transition.
///
/// trans_cel( cel_start,  cel_end,  rate ):
///   defines a cell transition [cel_start] -> [cel_end] occuring at the specified rate.
///
/// trans_cel( ref,  cel_start,  cel_end,  rate ):
///   defines a cell transition [cel_start] -> [cel_end] occuring at the specified rate, with a unique reference number 'ref'
///
/// trans_link( ref1,  ref2,  ncel, prob ):
///   creates a probability chain rule from doublet transition 'ref1' to doublet transition 'ref2'.
///   The occurrence of transition 'ref1' generates the occurence of transition 'ref2' with probability 'prob' in its surrounding area.
///   The value 'ncel' (1 or 2) specifies the position of the cell which is linking both transitions. Its state must be the same in the ending and starting doublets of the first and second transitions, respectively.
///   Chain of aligned doublets is not supported yet.
///
/// trans_regul( ref,  callback_reg) :
///   regulates the global rate of the doublet transition 'ref' by a callback function 'callback_reg'.
///
/// trans_cel_regul( ref,  callback_reg ):
///   regulates the global rate of the cell transition 'ref' by a callback function 'callback_reg'.
///
/// trans_check( ref,  callback_chk,  ncel ):
///   locally controls the doublet transition 'ref' by a callback function 'callback_chk'. The transition is blocked if a condition is not met, so that the condition has a positive effect on the transition.
///   The value 'ncel' (1 or 2) specifies the position of the cell to be tested in the starting doublet of the transition.
///
/// trans_check_inv( ref,  callback_chk,  ncel ):
///   locally controls the doublet transition 'ref' by a callback function 'callback_chk'. The transition is blocked if a condition is met, so that the condition has a negative effect on the transition.
///   The value 'ncel' (1 or 2) specifies the position of the cell to be tested in the starting doublet of the transition.
///
/// trans_check_color( ref,  ncel, lambda_col ):
///   declares that the transition 'ref' may have a different rate depending on the color of a cell.
///   The value 'ncel' (1 or 2) specifies the position of the cell to be tested in the starting doublet of the transition.
///   The argument 'lambda_col' is the rate of the transition for the colored cells.
///
/// trans_check_cell( ref, ncel, dir, cel_type ):
///   declares that the transition 'ref' will be blocked if a third cell located at 'dir' has type other than 'cel_type'
///
/// trans_check_no_cell( ref, ncel, dir, cel_type ):
///   declares that the transition 'ref' will be blocked if a third cell located at 'dir' has the type 'cel_type'
///
/// trans_type( ref, type ):
///   forces the type of the doublet transition 'ref'.
///   The only possible value for 'type' is TR_TRANSPORT. Generally, the transport transitions are automatically detected, unless the state of the moving cell changes.
///
/// trans_time( ref, time_mode ):
///   declares how time evolution is handled .
///   Possible values for 'time_mode' are TIME_EVOL (time evolution), TIME_CORR (time correction) and TIME_NO_EVOL (no time evolution).
///   By default, the value of time_mode is as follows:
///      - TIME_EVOL when the transition is not controlled,
///      - TIME_CORR when the transition is controlled by check_grad_vel() callback,
///      - TIME_NO_EVOL in all other cases (time correction not implemented).


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "defs.h"
#include "macros.h"
#include "param.h"
#include "transitions.h"
#include "space.h"
#include "lgca.h"

extern int Ncel[];    //number of cells in each state
extern int LNS;       //north-south width
extern int ava_trans;   //flag for avalanche transitions
extern float ava_delay;
extern int boundary;  //boundary conditions

/// parameters and rates for each model

char *model_name = NULL; //name of the model
#ifdef LGCA
int use_lgca = 1;
#else
int use_lgca = 0;
#endif // LGCA

#ifdef MODEL_DUN
float lambda_E, lambda_C, lambda_D, lambda_G, lambda_J;
float lambda_D_mob;
float lambda_T = -1;
float lambda_I = 1.0;
float coef_a = 0.0;
float coef_b = 10.0;
float coef_c = 500.0;
float prob_link_ET = 0.5;
float prob_link_TT = 0.5;
int flag_hm = 0;
#ifdef CELL_COLOR
float lambda_E_col = -1;
float lambda_C_col = -1;
float lambda_T_col = -1;
float lambda_D_col = -1;
float lambda_D_mob_col = -1;
#endif
#ifdef USE_VEGETATION
int use_veg = 1;
float lambda_V_up;
float lambda_V_die;
float lambda_V_die_air;
float lambda_V_die_mob;
float lambda_V_nuc;
float lambda_V_nuc_air;
float lambda_V_swarm;
#endif
#endif

#ifdef MODEL_AVA
float lambda_I = 0.1;
#endif

#ifdef MODEL_RIV
float lambda_I, lambda_D, lambda_L, lambda_C;
float lambda_E_v, lambda_E_h;
int trans_BT = 0;
#endif

#if defined(MODEL_CMB) || defined(MODEL_CRY)
float lambda_D, lambda_L, lambda_C;
#endif

#ifdef MODEL_ICB
float lambda_D, lambda_L, lambda_C, lambda_R;
#endif

#ifdef MODEL_DIF
float lambda_D, lambda_I;
#endif

#ifdef MODEL_D2G
float lambda_L, lambda_I;
#endif

#ifdef MODEL_LIFE
float lambda_B = 1.0;
float lambda_D = 1.0;
extern int B_min, B_max, S_min, S_max;
#endif

#ifdef AVALANCHES
float lambda_A = 1.0;
float lambda_A_unstable = 0;
#endif

float gravity = 0;

#if defined(MODEL_CMB) || defined(MODEL_RIV)
float flux_eo, flux_ns;  // intensite flux est-ouest, flux nord-sud
#else
float flux;                        // intensite flux horizontal
#endif

#if defined(LGCA) && defined(MODEL_CMB)
float prob_asv_D, prob_asv_L, prob_asv_C;  //asservissement des transitions vis-a-vis du flux
#endif


/// callbacks for the regulation of transitions
#ifdef MODEL_DUN
Callback_regul callback_bord_dun; //automatic reinjection of sand grains
Callback_regul callback_injection_coef; //regulation of the injection of sand grains
#endif
#ifdef ALTI
Callback_regul callback_ava_trans; //regulation of the transitions of avalanches
#endif

/// callbacks for the control of transitions
#ifdef ALTI
Callback_check check_alti;
Callback_check check_grain_seul;
Callback_check check_slope;  //checks the slope in the east-west direction (for the erosion)
Callback_check check_ava; //checks the value of the steepest slope in all the directions (for the avalanches)
#endif
#ifdef CELL_COLOR
Callback_check check_color; //checks the color of a cell
#endif
#ifdef LGCA
Callback_check check_mvt_bas;
Callback_check check_mvt_haut;
Callback_check check_mvt_est;
Callback_check check_mvt_ouest;
Callback_check check_mvt_est_et_bas;
Callback_check check_no_mvt_est_et_bas;
Callback_check check_mvt_mean_bas;
Callback_check check_mvt_mean_est;
Callback_check check_no_mvt_mean_bas;
Callback_check check_mvt_EB;
Callback_check check_grad_vel;
#ifdef CELL_COLOR
Callback_check check_grad_vel_color;
#endif
#endif
#ifdef USE_VEGETATION
Callback_check check_vegetation;
#endif
#ifdef MODEL_LIFE
Callback_check check_birth;
Callback_check check_death;
#endif
Callback_check check_cell_top;


void params_modele()
{
  char str[100];

  /// declare families of parameters
  sprintf(str, "%s model parameters", MOD_NAME);
  param_family("MODEL", str);
#ifdef CELL_COLOR
  param_family("COLOR", "Colored cells parameters");
#endif // CELL_COLOR

  /// declare parameters
#ifdef MODEL_DUN
  parameter("Lambda_E", "erosion rate", &lambda_E, PARAM_FLOAT, "MODEL");
  parameter("Lambda_C", "deposition rate", &lambda_C, PARAM_FLOAT, "MODEL");
  parameter("Lambda_T", "transport rate", &lambda_T, PARAM_FLOAT, "MODEL");
  parameter("Lambda_D", "diffusion rate", &lambda_D, PARAM_FLOAT, "MODEL");
  parameter("Lambda_D_mob", "diffusion rate of mobile grains", &lambda_D_mob, PARAM_FLOAT, "MODEL");
  parameter("Lambda_G", "gravitation rate", &lambda_G, PARAM_FLOAT, "MODEL");
  parameter("Lambda_I", "injection rate - optional", &lambda_I, PARAM_FLOAT, "MODEL");
  parameter("Coef_A", "coefficient for the vertical transport of mobile grains (0 by default)", &coef_a, PARAM_FLOAT, "MODEL");
  parameter("Coef_B", "coefficient for the deposition against an obstacle (10 by default)", &coef_b, PARAM_FLOAT, "MODEL");
  parameter("Coef_C", "coefficient for the deposition behind an obstacle (500 by default)", &coef_c, PARAM_FLOAT, "MODEL");
  parameter("Prob_link_ET", "probability of the link between erosion and upward transport (0.5 by default)", &prob_link_ET, PARAM_FLOAT, "MODEL");
  parameter("Prob_link_TT", "probability of the link between upward and horizontal transport (0.5 by default)", &prob_link_TT, PARAM_FLOAT, "MODEL");
  parameter("High_mobility", "additional transition link for higher mobility of grains (0|1) - optional", &flag_hm, PARAM_INT, "MODEL");
#ifdef LGCA
  parameter("Use_lgca", "use lattice gas automaton for the flow (YES|NO) - YES by default", &use_lgca, PARAM_BOOLEAN, "LGCA");
#endif
#ifdef CELL_COLOR
  parameter("Lambda_E_col", "erosion rate for white grains (Lambda_E by default)", &lambda_E_col, PARAM_FLOAT, "COLOR");
  parameter("Lambda_C_col", "deposition rate for white grains (Lambda_C by default)", &lambda_C_col, PARAM_FLOAT, "COLOR");
  parameter("Lambda_T_col", "transport rate for white grains (Lambda_T by default)", &lambda_T_col, PARAM_FLOAT, "COLOR");
  parameter("Lambda_D_col", "diffusion rate for white grains (Lambda_D by default)", &lambda_D_col, PARAM_FLOAT, "COLOR");
  parameter("Lambda_D_mob_col", "diffusion rate for white mobile grains (Lambda_D_mob by default)", &lambda_D_mob_col, PARAM_FLOAT, "COLOR");
#endif
#ifdef USE_VEGETATION
  param_family("VEGETATION", "Vegetation parameters");
  parameter("Use_vegetation", "use vegetated cells and transitions (YES|NO) - YES by default", &use_veg, PARAM_BOOLEAN, "VEGETATION");
  parameter("Lambda_V_up", "vegetation upward diffusion rate", &lambda_V_up, PARAM_FLOAT, "VEGETATION");
  parameter("Lambda_V_die", "vegetation dying rate under the sand", &lambda_V_die, PARAM_FLOAT, "VEGETATION");
  parameter("Lambda_V_die_air", "vegetation dying rate in air", &lambda_V_die_air, PARAM_FLOAT, "VEGETATION");
  parameter("Lambda_V_die_mob", "vegetation dying rate due to mobile grain", &lambda_V_die_mob, PARAM_FLOAT, "VEGETATION");
  parameter("Lambda_V_nuc", "vegetation nucleation rate on the ground", &lambda_V_nuc, PARAM_FLOAT, "VEGETATION");
  parameter("Lambda_V_nuc_air", "vegetation nucleation rate from air", &lambda_V_nuc_air, PARAM_FLOAT, "VEGETATION");
  parameter("Lambda_V_swarm", "vegetation swarming rate", &lambda_V_swarm, PARAM_FLOAT, "VEGETATION");
#endif
  //parameter("Lambda_J", "non-utilisé", &lambda_J, PARAM_FLOAT);

#endif
#ifdef MODEL_AVA
  parameter("Lambda_I", "injection rate - optional", &lambda_I, PARAM_FLOAT, "MODEL");
#endif //MODEL_AVA

#ifdef MODEL_RIV
  parameter("Lambda_D", "diffusion rate", &lambda_D, PARAM_FLOAT, "MODEL");
  parameter("Lambda_L", "dissolution rate", &lambda_L, PARAM_FLOAT, "MODEL");
  parameter("Lambda_E_v", "rate of vertical erosion by transport", &lambda_E_v, PARAM_FLOAT, "MODEL");
  parameter("Lambda_E_h", "rate of horizontal erosion by transport", &lambda_E_h, PARAM_FLOAT, "MODEL");
  parameter("Lambda_C", "deposition rate", &lambda_C, PARAM_FLOAT, "MODEL");
  //parameter("Lambda_G", "gravitation rate", &lambda_G, PARAM_FLOAT, "MODEL");
  parameter("Gravity", "gravity coefficient", &gravity, PARAM_FLOAT, "MODEL");
  parameter("Lambda_I", "injection rate - optional", &lambda_I, PARAM_FLOAT, "MODEL");
  parameter("Trans_BT", "transitions added for the BT state (0|1) - optional", &trans_BT, PARAM_INT, "MODEL");
  parameter("Flux_eo", "flow intensity in the east-west direction - optional", &flux_eo, PARAM_FLOAT, "MODEL");
  parameter("Flux_ns", "flow intensity in the north-south direction - optional", &flux_ns, PARAM_FLOAT, "MODEL");
  //parameter("Lambda_I", "debit a la source", &lambda_I, PARAM_FLOAT, "MODEL");
#endif

#ifdef MODEL_CMB
  //parameter("Lambda_I", "coefficient d'infiltration", &lambda_I, PARAM_FLOAT, "MODEL");
  parameter("Lambda_D", "diffusion rate", &lambda_D, PARAM_FLOAT, "MODEL");
  parameter("Lambda_L", "dissolution rate", &lambda_L, PARAM_FLOAT, "MODEL");
  parameter("Lambda_C", "crystallization rate", &lambda_C, PARAM_FLOAT, "MODEL");
  parameter("Gravity", "gravity coefficient", &gravity, PARAM_FLOAT, "MODEL");
  parameter("Flux_eo", "flow intensity in the east-west direction - optional", &flux_eo, PARAM_FLOAT, "MODEL");
  parameter("Flux_ns", "flow intensity in the north-south direction - optional", &flux_ns, PARAM_FLOAT, "MODEL");
#endif

#ifdef MODEL_CRY
  //parameter("Lambda_I", "coefficient d'infiltration", &lambda_I, PARAM_FLOAT);
  parameter("Lambda_D", "diffusion rate", &lambda_D, PARAM_FLOAT, "MODEL");
  //parameter("Lambda_L", "coefficient de dissolution", &lambda_L, PARAM_FLOAT);
  parameter("Lambda_C", "crystallization rate", &lambda_C, PARAM_FLOAT, "MODEL");
  //parameter("Gravite", "coefficient de gravite", &lambda_G, PARAM_FLOAT);
  parameter("Flux", "flow intensity in the east-west direction - optional", &flux, PARAM_FLOAT, "MODEL");
#endif

#ifdef MODEL_ICB
  //parameter("lambda_P", "coefficient de compaction", &lambda_P, PARAM_FLOAT, "MODEL");
  parameter("Lambda_R", "cooling rate", &lambda_R, PARAM_FLOAT, "MODEL");
  parameter("Lambda_D", "diffusion rate", &lambda_D, PARAM_FLOAT, "MODEL");
  parameter("Lambda_L", "dissolution rate", &lambda_L, PARAM_FLOAT, "MODEL");
  parameter("Lambda_C", "crystallization rate", &lambda_C, PARAM_FLOAT, "MODEL");
  parameter("Gravity", "gravity coefficient", &gravity, PARAM_FLOAT, "MODEL");
  parameter("Flux", "intensity of horizontal flow - optional", &flux, PARAM_FLOAT, "MODEL");
#endif

#ifdef MODEL_DIF
  parameter("Lambda_D", "diffusion rate", &lambda_D, PARAM_FLOAT, "MODEL");
  parameter("Lambda_I", "injection rate - optional", &lambda_I, PARAM_FLOAT, "MODEL");
#endif

#ifdef MODEL_D2G
  parameter("Lambda_L", "diffusion rate", &lambda_L, PARAM_FLOAT, "MODEL");
  parameter("Lambda_I", "injection rate", &lambda_I, PARAM_FLOAT, "MODEL");
#endif

#ifdef MODEL_LIFE
  parameter("Lambda_B", "birth rate", &lambda_B, PARAM_FLOAT, "MODEL");
  parameter("Lambda_D", "death rate", &lambda_D, PARAM_FLOAT, "MODEL");
  parameter("B_min", "min. number of neighbours for birth transitions", &B_min, PARAM_INT, "MODEL");
  parameter("B_max", "max. number of neighbours for birth transitions", &B_max, PARAM_INT, "MODEL");
  parameter("S_min", "min. number of neighbours that prevents death transitions", &S_min, PARAM_INT, "MODEL");
  parameter("S_max", "max. number of neighbours that prevents death transitions", &S_max, PARAM_INT, "MODEL");
#endif
}

/*
void model_parse()
{
}
*/

void init_modele()
{
  //params_modele();

  //model_parse();

  if (model_name && strcmp(model_name, MOD_NAME)){
    ErrPrintf("ERROR: Wrong model name %s in parameter file. Check model in defs.h and rebuild all.\n", model_name);
    exit(-2);
  }

#if defined(MODEL_CMB) || defined(MODEL_RIV)

  double flux_e, flux_o, flux_n, flux_s;

  if (flux_eo >= 0){
    flux_o = 0;
    flux_e = flux_eo;
  }
  else{
    flux_o = -flux_eo;
    flux_e = 0;
  }

  if (flux_ns >= 0){
    flux_n = 0;
    //flux_n = flux_ns;//0; //ortho
    flux_s = flux_ns;
  }
  else{
    flux_n = -flux_ns;
    flux_s = 0;
  }
#else
  double flux1, flux2;

  if (flux >= 0){
    flux1 = 0;
    flux2 = flux;
  }
  else{
    flux1 = -flux;
    flux2 = 0;
  }
#endif

  gravity += 1.0; //= 1.05;

  //printf("ReSCAL - %s model\n",MOD_NAME);


/*****************************************************************************/
/********************************* DUN model *********************************/
/*****************************************************************************/
#ifdef MODEL_DUN

/* Parameters:
  lambda_E   	erosion
  lambda_C    deposition
  lambda_T    transport
  lambda_D	  diffusion (NS)
  lambda_G	  gravity
  lambda_I    injection
  lambda_A    avalanches
  lambda_J    not used
  _______________________________________________

  7 lambdas + 1 not used

*/

  char tmode = TIME_CORR; //time correction (default value for check_grad_vel)

  /***** erosion *****/
  trans_ref(8,  EST_OUEST, EAUC, GR,  EAUC, GRJ, lambda_E );
  trans_ref(7,  VERTICAL,  EAUC, GR,  GRJ,  EAUC, lambda_E/100000.0 );
  trans_type(7, TR_TRANSPORT);
#ifdef CELL_COLOR
  //erosion of colored grains
  if (lambda_E_col >= 0){
    tmode = TIME_NO_EVOL; // time correction not possible here !
    trans_check_color(8, 2, lambda_E_col);
    trans_check_color(7, 2, lambda_E_col/100000.0);
  }
#endif

#ifdef CGV
  //control of the erosion by shear stress
  Callback_check *pchk = check_grad_vel;
  if (use_lgca){
#ifdef CELL_COLOR
    extern float grdvc_max_col;
    if (grdvc_max_col >= 0){
      LogPrintf("check_grad_vel_color() attached to the erosion transitions\n");
      pchk = check_grad_vel_color;
    }
#endif //CELL_COLOR
    trans_check(8, pchk, 2);
    trans_check(7, pchk, 2);
    trans_time(8, tmode);
    trans_time(7, tmode);
  }
#endif //CGV

  /***** gravity *****/
  trans(VERTICAL, GR, EAUC, EAUC,  GR, lambda_G);
  trans(VERTICAL, GRJ, EAUC, EAUC, GRJ, (coef_a>0) ? lambda_T/coef_a : lambda_G );

  /***** transport *****/
  if (lambda_T<0) lambda_T = 3/*6*/*lambda_C; //default value of lambda_T (with correlations of transitions)
  trans_ref(10,  VERTICAL,  EAUC, GRJ,  GRJ,  EAUC, coef_a*lambda_T );
  //trans_check(10, check_slope, 2, TIME_NO_EVOL);
  //trans_ref(11,  EST_OUEST,  GRJ,  EAUC,  EAUC, GRJ, 6/*6*/*lambda_C );
  trans_ref(11,  EST_OUEST,  GRJ,  EAUC,  EAUC, GRJ, lambda_T );
  //trans( NORD_SUD,  GRJ,  EAUC,  EAUC, GRJ, lambda_T/10 ); //test diffusion NS
  //trans( NORD_SUD,  EAUC, GRJ,   GRJ,  EAUC, lambda_T/10 ); //test diffusion NS
  if (flag_hm) trans_ref(50,  EST_OUEST, GRJ, GR, GRJ, GR, lambda_T );
#ifdef CELL_COLOR
  //transport of colored grains
  if (lambda_T_col >= 0){
    trans_check_color(11, 1, lambda_T_col);
  }
#endif

  /***** chains of transitions for saltation *****/
  trans_link(8, 10, 2, prob_link_ET/*0.5*/);
  trans_link(10, 11, 1, prob_link_TT/*0.5*/);
  trans_link(7, 11, 1, 1.0);
  if (flag_hm) trans_link(50, 10, 1, prob_link_TT/*0.5*/);

  /***** deposition *****/
  trans_ref(9,  EST_OUEST, GRJ, EAUC, GR, EAUC, lambda_C );
  trans_ref(16,  EST_OUEST, GRJ, GR, GR, GR, coef_b*lambda_C );
  trans_ref(17,  EST_OUEST, GRJ, DUM, GR, DUM, coef_b*lambda_C );
  trans(  EST_OUEST, GR, GRJ, GR, GR, coef_c*lambda_C/*1000*/ );
#ifdef CELL_COLOR
  //deposition of colored grains
  if (lambda_C_col >= 0){
    trans_check_color(9, 1, lambda_C_col);
    trans_check_color(16, 1, 10*lambda_C_col);
    trans_check_color(17, 1, 10*lambda_C_col);
  }
#endif

  /***** diffusion *****/
  if (LNS > 1){
    trans_ref(12, NORD_SUD,   GR, EAUC, EAUC,  GR, lambda_D );
    trans_ref(13, NORD_SUD,   EAUC,  GR, GR, EAUC, lambda_D );
#ifndef CELL_COLOR
    trans(NORD_SUD,   GRJ, EAUC, EAUC,  GRJ, lambda_D_mob);
    trans(NORD_SUD,   EAUC,  GRJ, GRJ, EAUC, lambda_D_mob);
#else
    if ((lambda_D_mob > 0) || (lambda_D_mob_col > 0)){
      trans_ref(19, NORD_SUD,   GRJ, EAUC, EAUC,  GRJ, lambda_D_mob);
      trans_ref(20, NORD_SUD,   EAUC,  GRJ, GRJ, EAUC, lambda_D_mob);
    }
    //diffusion of colored grains
    if (lambda_D_col >= 0){
      trans_check_color(12, 1, lambda_D_col);
      trans_check_color(13, 2, lambda_D_col);
    }
    if (lambda_D_mob_col >=0){
      trans_check_color(19, 1, lambda_D_mob_col);
      trans_check_color(20, 2, lambda_D_mob_col);
    }
#endif
#ifdef CGV
    if (use_lgca){
      //control of the diffusion by shear stress
      trans_check(12, pchk, 1);
      trans_check(13, pchk, 2);
      //trans_check(19, check_grad_vel, 1);
      //trans_check(20, check_grad_vel, 2);
    }
#endif //CGV
  }

#ifdef AVALANCHES
  /***** avalanches *****/
  if (ava_trans){
    /*if (LNS > 1){
      trans_ref(14, HORIZONTAL, GR, EAUC, EAUC, GR, lambda_A);
      trans_ref(15, HORIZONTAL, EAUC, GR, GR, EAUC, lambda_A);
    }
    else{
      trans_ref(14, EST_OUEST, GR, EAUC, EAUC, GR, lambda_A);
      trans_ref(15, EST_OUEST, EAUC, GR, GR, EAUC, lambda_A);
    }*/
    trans_ref(14, HORIZONTAL, GR, EAUC, EAUC, GR, lambda_A);
    trans_ref(15, HORIZONTAL, EAUC, GR, GR, EAUC, lambda_A);
    trans_check(14, check_ava, 1);
    trans_check(15, check_ava, 2);
    if (ava_delay){
      trans_regul(14, callback_ava_trans);
      trans_regul(15, callback_ava_trans);
    }
  }
#endif

  /***** boundary transitions *****/
  if ((boundary == BC_OUT) || (boundary==BC_REINJECTION)){
    LogPrintf("transitions sur le bord\n");
    //trans(EST_OUEST,  BORD, GR,  BORD, GRJ, 0.01);
    //trans(EST_OUEST,  BORD, EAUC,  BORD, GRJ, 0.0002);
    trans(EST_OUEST,  GRJ,  BORD, EAUC, BORD, 20);
    trans(EST_OUEST,  GR,  BORD, EAUC, BORD, 20);
  }
  else if (boundary == BC_OPEN){
    /// only limited transport
    LogPrintf("transitions lentes sur le bord\n");
    trans(EST_OUEST,  GRJ,  BORD, EAUC, BORD, 1);
  }

  /***** injection *****/
  //injection of sand grains (for the computation of Qsat and Lsat)
  if (Ncel[IN]>0){
    trans_ref(2, VERTICAL, IN,  EAUC, IN, GR, lambda_I/*0.5*/); // value 2 sets the (arbitrary) reference number of the transition
    trans(VERTICAL, EAUC, IN, GRJ, IN, lambda_I/*0.05*/);
  }

  /***** output *****/
  //output of sand
  if (Ncel[OUT]>0){
    trans(VERTICAL, GR, OUT, EAUC, OUT, 1000);
    trans(VERTICAL, GRJ, OUT, EAUC, OUT, 1000);
  }

  /***** reinjection *****/
  if ((boundary==BC_REINJECTION) && (Ncel[IN]>0)){
    //reinjection of sand grains
    LogPrintf("reinjection de grains\n");
    trans_regul(2, callback_bord_dun);

    //temporal regulation of the reinjection rate
    //trans_regul(2, callback_injection_coef);
  }

  /***** vegetation *****/
#ifdef USE_VEGETATION
//#define VEG_VARIANT 1 /// VEG state is considered as vegetation without grain
#define VEG_VARIANT 2 /// VEG state is considered as vegetated grain

  extern int veg_h_max;

  if (use_veg){
    LogPrintf("VEG_VARIANT=%d\n", VEG_VARIANT);

  #if VEG_VARIANT==1
    /// upward diffusion
    trans_ref(25, VERTICAL, GR, VEG, VEG, GR, lambda_V_up);
    if (veg_h_max > 0) trans_check(25, check_vegetation, 2);

    /// vegetation may die in the sand
    trans(VERTICAL, GR, VEG, EAUC, GR, lambda_V_die);

    /// gravity
    trans(VERTICAL, VEG, EAUC, EAUC, VEG, lambda_G);
    //trans(VERTICAL, VEG, GR, GR, VEG, lambda_G);

    /// mobile grains are stabilized by the vegetation
    trans(EST_OUEST, GRJ, VEG, GR, VEG, 10*lambda_C);
  #endif //VEG_VARIANT

  #if VEG_VARIANT==2
    /// nucleation of vegetation
    trans_ref(24, VERTICAL, EAUC, GR, EAUC, GRV, lambda_V_nuc_air);
    trans(VERTICAL, EAUC, DUM, EAUC, VEG, lambda_V_nuc);

    /// upward diffusion
    trans_ref(25, VERTICAL, GR, GRV, GRV, GR, lambda_V_up);
    //trans(VERTICAL, GR, DUM, VEG, DUM, 0.0003);
    //trans(VERTICAL, GR, DUM, GRV, DUM, lambda_V_up);
    //trans(VERTICAL, GR, DV, VEG, DUM, 0.0005);
    /*if (Ncel[VEG]>0)*/ trans_ref(26, VERTICAL, GR, VEG, GRV, DUM, lambda_V_up);
    if (veg_h_max > 0){
      trans_check(24, check_vegetation, 2);
      trans_check(25, check_vegetation, 2);
    }

    /// death of vegetation
    trans(VERTICAL, GR, GRV, GR, GR, lambda_V_die);
    trans(VERTICAL, GR, VEG, GR, DUM, lambda_V_die);
    trans_ref(35,VERTICAL, EAUC, GRV, EAUC, GR, lambda_V_die_air);
    trans_ref(36,VERTICAL, EAUC, VEG, EAUC, DUM, lambda_V_die_air);
  #ifdef CGV
    if (use_lgca){
      trans_check(35, pchk, 2);
      trans_check(36, pchk, 2);
    }
  #endif // CGV
    trans_ref(27, VERTICAL, GRV, GRV, GRV, GR, lambda_V_die);
    trans_ref(28, VERTICAL, GRV, VEG, GRV, DUM, lambda_V_die);
    trans(EST_OUEST, GRJ, GRV, GR, GR, lambda_V_die_mob);

    /// gravity
    trans(VERTICAL, GRV, EAUC, EAUC, GRV, lambda_G);
    //trans(VERTICAL, VEG, GR, GR, VEG, lambda_G);

    /// mobile grains are stabilized by the vegetated grains
    trans(EST_OUEST, GRJ, GRV, GR, GRV, 10*lambda_C);

    /// swarming
    if (lambda_V_swarm > 0){
      //trans(HORIZONTAL, GRV, GR, GRV, GRV, lambda_V_swarm);
      //trans(HORIZONTAL, GR, GRV, GRV, GRV, lambda_V_swarm);
      /// horizontal swarming of GRV
      trans_ref(29, HORIZONTAL, GRV, GR, GRV, GRV, lambda_V_swarm);
      trans_ref(30, HORIZONTAL, GR, GRV, GRV, GRV, lambda_V_swarm);
      //trans_ref(29, EST_OUEST, GRV, GR, GRV, GRV, lambda_V_swarm);
      //trans_ref(30, EST_OUEST, GR, GRV, GRV, GRV, lambda_V_swarm);
  #ifdef CGV
      if (use_lgca){
        trans_check_inv(29, pchk, 2);
        trans_check_inv(30, pchk, 1);
      }
  #endif // CGV
      /// upward swarming
      trans_check_no_cell(29, 2, HAUT, GRV);
      trans_check_no_cell(30, 1, HAUT, GRV);
      trans_link(29, 25, 2, 1.0);
      trans_link(30, 25, 1, 1.0);
      //trans_link(29, 27, 2, 1.0);
      //trans_link(30, 27, 1, 1.0);
      /// downward swarming
      trans_ref(31, HORIZONTAL, GRV, EAUC, GRV, EAUC, lambda_V_swarm);
      trans_ref(32, HORIZONTAL, EAUC, GRV, EAUC, GRV, lambda_V_swarm);
      trans_check_cell(31, 2, BAS, GR);
      trans_check_cell(32, 1, BAS, GR);
  #ifdef CGV
      if (use_lgca){
        trans_check_inv(31, pchk, 2);
        trans_check_inv(32, pchk, 1);
      }
  #endif // CGV
      trans_link(31, 24, 2, 1.0);
      trans_link(32, 24, 1, 1.0);
      /// horizontal swarming of VEG
      trans_ref(33, HORIZONTAL, VEG, DUM, VEG, VEG, lambda_V_swarm);
      trans_ref(34, HORIZONTAL, DUM, VEG, VEG, VEG, lambda_V_swarm);
      /// upward swarming VEG -> GRV
      trans_check_no_cell(33, 2, HAUT, GRV);
      trans_check_no_cell(34, 1, HAUT, GRV);
      trans_link(33, 26, 2, 1.0);
      trans_link(34, 26, 1, 1.0);
      //trans_link(33, 28, 2, 1.0);
      //trans_link(34, 28, 1, 1.0);
    }
  #endif //VEG_VARIANT
  }
#endif //USE_VEGETATION

#endif //MODEL_DUN


/*****************************************************************************/
/********************************* AVA model *********************************/
/*****************************************************************************/
#ifdef MODEL_AVA
  /***** gravity *****/
  trans(VERTICAL, GR, AIR, AIR,  GR, 1000);

  /***** injection *****/
  trans(VERTICAL, IN, AIR, IN, GR, lambda_I);

  /***** output *****/
  //trans(VERTICAL, GR, OUT, AIR, OUT, 1000);
  trans(VERTICAL, GR, OUT, AIR, OUT, 0.1);

  /***** avalanches *****/
  if (ava_trans){
#if 1 //1 pour conserver le couplage des transitions
    trans_ref(0, HORIZONTAL, GR, AIR, AIR, GR, lambda_A);
    trans_ref(1, HORIZONTAL, AIR, GR, GR, AIR,  lambda_A);
    trans_check(0, check_ava, 1);
    trans_check(1, check_ava, 2);
    if (ava_delay){
      trans_regul(0, callback_ava_trans);
      trans_regul(1, callback_ava_trans);
    }
#else //decouplage
    trans_ref(0, EST_OUEST, GR, AIR, AIR, GR, lambda_A);
    trans_ref(1, EST_OUEST, AIR, GR, GR, AIR,  lambda_A);
    trans_ref(2, NORD_SUD, GR, AIR, AIR, GR, lambda_A);
    trans_ref(3, NORD_SUD, AIR, GR, GR, AIR,  lambda_A);
    trans_check(0, check_ava, 1);
    trans_check(1, check_ava, 2);
    trans_check(2, check_ava, 1);
    trans_check(3, check_ava, 2);
#endif
  }
#endif //MODEL_AVA


/*****************************************************************************/
/********************************* CMB model *********************************/
/*****************************************************************************/
#ifdef MODEL_CMB

  /***** dissolution *****/
  trans(  ISOTROPE,  PLUS, MOINS,  ZERO,  ZERO, lambda_L );
  trans(  ISOTROPE,  PLUS, MOINS,  PLUS,  ZERO, lambda_L );

  /***** crystallization *****/
  trans(  ISOTROPE,  PLUS,  ZERO,  PLUS,  PLUS, lambda_C );
  trans(  ISOTROPE,  PLUS,  ZERO,  PLUS, MOINS, lambda_C*2 );

  /***** transport *****/
  trans(  VERTICAL, MOINS,  ZERO,  ZERO, MOINS, lambda_D*gravity);
  trans(  VERTICAL,  ZERO, MOINS, MOINS,  ZERO, lambda_D/gravity);

  if (flux_eo || flux_ns){
    trans(EST_OUEST, MOINS,  ZERO,  ZERO, MOINS, lambda_D+flux_o );
    trans(NORD_SUD, MOINS,  ZERO,  ZERO, MOINS, lambda_D+flux_n );
    trans(EST_OUEST,  ZERO, MOINS, MOINS,  ZERO, lambda_D+flux_e );
    trans(NORD_SUD,  ZERO, MOINS, MOINS,  ZERO, lambda_D+flux_s );
  }
  else{
    trans(HORIZONTAL, MOINS,  ZERO,  ZERO, MOINS, lambda_D );
    trans(HORIZONTAL,  ZERO, MOINS, MOINS,  ZERO, lambda_D );
  }

/*
  //etat PIERRE
  trans(  ISOTROPE,PIERRE,  ZERO,PIERRE,  PLUS, lambda_C );
  trans(  ISOTROPE,PIERRE,  ZERO,PIERRE, MOINS, lambda_C*2 );

  //PIERRE = source de ZERO
  trans(  ISOTROPE,PIERRE, MOINS,PIERRE,  ZERO, 0.001 );
*/
  //porosite PLUS/ZERO
  //trans(  ISOTROPE,  ZERO,  PLUS,  PLUS,  ZERO, 0.001 );

  /***** boundary transitions *****/
  if (boundary == BC_OPEN){
    LogPrintf("transitions sur le bord\n");
    //trans(  VERTICAL,  BORD,  ZERO,  BORD, MOINS, lambda_D*gravite );

    //modelisation fracture
    //trans(  VERTICAL,  BORD,  MOINS,  BORD, ZERO, 0.01 );
    //trans(  VERTICAL,  ZERO,  BORD,  MOINS, BORD, lambda_D/gravite );
    //trans(  VERTICAL,  BORD,  MOINS,  BORD, ZERO, 0.1/*lambda_D/gravite*/ );
    //trans(  VERTICAL,  ZERO,  BORD,  MOINS, BORD, 0.01 );

    trans(  ISOTROPE,  BORD,  MOINS, BORD,  ZERO, 0.1/*lambda_D/gravite*/ );

    //trans(  VERTICAL,  BORD,  MOINS, BORD,  ZERO, 0.1/*lambda_D/gravite*/ );
    //trans(  VERTICAL, MOINS,  BORD,  ZERO,  BORD, 0.1 );
    //trans(HORIZONTAL,  BORD,  MOINS, BORD,  ZERO, 0.1 );
    //trans(HORIZONTAL, MOINS,  BORD,  ZERO,  BORD, 0.1 );

    //porosite PLUS/ZERO aux bords
    //trans(  ISOTROPE,  BORD,  PLUS,  BORD,  ZERO, 0.0001 );
  }

#endif //MODEL_CMB


/*****************************************************************************/
/********************************* ICB model *********************************/
/*****************************************************************************/
#ifdef MODEL_ICB

  /***** dissolution *****/
  trans(  VERTICAL,  PLUS, MOINS,  ZERO,  ZERO, lambda_L );
  trans(  VERTICAL,  PLUS, MOINS,  PLUS,  ZERO, lambda_L*2 );
  trans(  VERTICAL, MOINS,  PLUS,  ZERO,  ZERO, lambda_L );
  trans(  VERTICAL, MOINS,  PLUS,  ZERO,  PLUS, lambda_L*2 );
  trans(HORIZONTAL,  PLUS, MOINS,  ZERO,  ZERO, lambda_L );
  trans(HORIZONTAL,  PLUS, MOINS,  PLUS,  ZERO, lambda_L*2 );
  trans(HORIZONTAL, MOINS,  PLUS,  ZERO,  ZERO, lambda_L );
  trans(HORIZONTAL, MOINS,  PLUS,  ZERO,  PLUS, lambda_L*2 );

  /***** crystallization *****/
  trans(  VERTICAL,  PLUS,  ZERO,  PLUS,  PLUS, lambda_C );
  trans(  VERTICAL,  PLUS,  ZERO,  PLUS, MOINS, lambda_C*3 );
  trans(  VERTICAL,  ZERO,  PLUS,  PLUS,  PLUS, lambda_C );
  trans(  VERTICAL,  ZERO,  PLUS, MOINS,  PLUS, lambda_C*3 );
  trans(HORIZONTAL,  PLUS,  ZERO,  PLUS,  PLUS, lambda_C );
  trans(HORIZONTAL,  PLUS,  ZERO,  PLUS, MOINS, lambda_C*3 );
  trans(HORIZONTAL,  ZERO,  PLUS,  PLUS,  PLUS, lambda_C );
  trans(HORIZONTAL,  ZERO,  PLUS, MOINS,  PLUS, lambda_C*3 );

  /***** transport *****/
  trans(  VERTICAL, MOINS,  ZERO,  ZERO, MOINS, lambda_D/gravity );
  trans(  VERTICAL,  ZERO, MOINS, MOINS,  ZERO, lambda_D*gravity );
  trans(HORIZONTAL, MOINS,  ZERO,  ZERO, MOINS, lambda_D+flux1 );
  trans(HORIZONTAL,  ZERO, MOINS, MOINS,  ZERO, lambda_D+flux2 );

  /***** cooling *****/
  //LogPrintf("version avec refroidissement\n");
  //LogPrintf("lambda_R = %f\n",lambda_R);
  trans_cel(MOINS, ZERO, lambda_R);

#endif

/*****************************************************************************/
/********************************* RIV model *********************************/
/*****************************************************************************/
#ifdef MODEL_RIV
  double coef_v_L = 1;//0.01;
  double coef_h_L = 1;//100; //dans le sens du flux seulement
  double coef_h_C = 1;//0.01;
  //double coef_h_C2 = 1;//0.1;
  double coef_C = 1;//0.01

  /// parameters for BT transitions
  double coef_C_BT = 1;//10;//1
  double coef_D_BT = 10.0;

  double prob_v_BOUE_BT = 1.0;  //probabilite pour que BOUE devienne BT lors d'un transport vertical
  double prob_v_BOUE_BOUE = (1.0-prob_v_BOUE_BT);
  double prob_h_BT_BOUE = 0.1;  //probabilite pour que BT devienne BOUE lors d'un transport horizontal
  double prob_h_BT_BT = (1.0-prob_h_BT_BOUE);

  LogPrintf("coef_v_L = %f\n",coef_v_L);
  LogPrintf("coef_h_L = %f\n",coef_h_L);
  LogPrintf("coef_h_C = %f\n",coef_h_C);

  if (trans_BT){
    LogPrintf("Parameters for BT transitions:\n");
    LogPrintf("coef_C_BT = %f\n", coef_C_BT);
    LogPrintf("coef_D_BT = %f\n", coef_D_BT);
    LogPrintf("prob_v_BOUE_BT = %f\n", prob_v_BOUE_BT);
    LogPrintf("prob_v_BOUE_BOUE = %f\n", prob_v_BOUE_BOUE);
    LogPrintf("prob_h_BT_BOUE = %f\n", prob_h_BT_BOUE);
    LogPrintf("prob_h_BT_BT = %f\n", prob_h_BT_BT);
  }

  /***** transport *****/
#ifndef LISSAGE
  trans(  VERTICAL, EAU,  BOUE,  BOUE, EAU, lambda_D/gravity );
#endif
  if (!trans_BT){
    trans_ref(0,  VERTICAL,  BOUE, EAU, EAU,  BOUE, lambda_D*gravity );
  }
  else{
    trans_ref(0,  VERTICAL,  BOUE, EAU, EAU,  BOUE, lambda_D*gravity*prob_v_BOUE_BOUE );
    trans_ref(0,  VERTICAL,  BOUE, EAU, EAU,  BT, lambda_D*gravity*prob_v_BOUE_BT );
  }

  //trans_ref(1, HORIZONTAL, EAU,  BOUE,  BOUE, EAU, lambda_D+flux1 );
  //trans_ref(2, HORIZONTAL,  BOUE, EAU, EAU,  BOUE, lambda_D+flux2 );
  if (flux_eo == flux_ns){
    trans_ref(1, HORIZONTAL, EAU,  BOUE,  BOUE, EAU, lambda_D+flux_o );
    trans_ref(2, HORIZONTAL,  BOUE, EAU, EAU,  BOUE, lambda_D+flux_e );
  }
  else{
    trans_ref(1, EST_OUEST, EAU,  BOUE,  BOUE, EAU, lambda_D+flux_o );
    trans_ref(2, EST_OUEST,  BOUE, EAU, EAU,  BOUE, lambda_D+flux_e );
    trans(NORD_SUD, EAU,  BOUE,  BOUE, EAU, lambda_D+flux_n );
    trans(NORD_SUD,  BOUE, EAU, EAU,  BOUE, lambda_D+flux_s );
  }

  /***** dissolution *****/
  trans_ref(3,  VERTICAL,  TERRE, EAU,  BOUE,  BOUE, lambda_L*coef_v_L );
  trans_ref(4,  VERTICAL,  TERRE, EAU,  TERRE,  BOUE, lambda_L*coef_v_L );
  trans_ref(5,  VERTICAL, EAU,  TERRE,  BOUE,  BOUE, lambda_L*coef_v_L );
  trans_ref(6,  VERTICAL, EAU,  TERRE,  BOUE,  TERRE, lambda_L*coef_v_L );
  trans_ref(7, HORIZONTAL,  TERRE, EAU,  BOUE,  BOUE, lambda_L );
  trans_ref(8, HORIZONTAL,  TERRE, EAU,  TERRE,  BOUE, lambda_L );
  trans_ref(9, HORIZONTAL, EAU,  TERRE,  BOUE,  BOUE, lambda_L*coef_h_L );
  trans_ref(10, HORIZONTAL, EAU,  TERRE,  BOUE,  TERRE, lambda_L*coef_h_L );

  /// dilution BOUE/EAU
  //trans(  VERTICAL,  EAU,  BOUE, EAU,  EAU, 0.01 );

  /***** deposition *****/
  trans(  VERTICAL,  TERRE,  BOUE,  TERRE,  TERRE, lambda_C*coef_C);
  trans(  VERTICAL,  BOUE,  TERRE,  TERRE,  TERRE, lambda_C*coef_C);
#ifndef LISSAGE
  trans(  VERTICAL,  TERRE,  BOUE,  TERRE, EAU, lambda_C*2*coef_C);
  trans(  VERTICAL,  BOUE,  TERRE, EAU,  TERRE, lambda_C*2*coef_C);
#endif
  trans(HORIZONTAL,  TERRE,  BOUE,  TERRE,  TERRE, lambda_C*coef_h_C*coef_C);
  trans(HORIZONTAL,  BOUE,  TERRE,  TERRE,  TERRE, lambda_C*coef_h_C*coef_C);
#ifndef LISSAGE
  trans(HORIZONTAL,  TERRE,  BOUE,  TERRE, EAU, lambda_C*2*coef_h_C*coef_C);
  trans(HORIZONTAL,  BOUE,  TERRE, EAU,  TERRE, lambda_C*2*coef_h_C*coef_C);
#endif
/*
  trans(  VERTICAL,  BOUE,  BOUE,  BOUE,  TERRE, lambda_C ); //test deposition BOUE+BOUE
  trans(  VERTICAL,  BOUE,  BOUE,  EAU,  BOUE, lambda_C*2 ); //test deposition BOUE+BOUE
  trans(HORIZONTAL,  BOUE,  BOUE,  TERRE,  BOUE, lambda_C*coef_h_C2 ); //test deposition BOUE+BOUE
  trans(HORIZONTAL,  BOUE,  BOUE,  BOUE,  TERRE, lambda_C*coef_h_C2 ); //test deposition BOUE+BOUE
  trans(HORIZONTAL,  BOUE,  BOUE,  BOUE,  EAU, lambda_C*2*coef_h_C2 ); //test deposition BOUE+BOUE
  trans(HORIZONTAL,  BOUE,  BOUE,  EAU,  BOUE, lambda_C*2*coef_h_C2 ); //test deposition BOUE+BOUE
*/

  /***** BT transitions *****/
  if (trans_BT){
    /// transport
    trans(  VERTICAL, EAU,  BT,  BT, EAU, lambda_D/gravity );
    trans_ref(20,  VERTICAL,  BT, EAU, EAU,  BT, lambda_D*gravity );
    //trans_ref(1, HORIZONTAL, EAU,  BOUE,  BOUE, EAU, lambda_D+flux1 );
    //trans_ref(2, HORIZONTAL,  BOUE, EAU, EAU,  BOUE, lambda_D+flux2 );
    if (flux_eo == flux_ns){
      trans_ref(21, HORIZONTAL, EAU,  BT,  BT, EAU, (lambda_D+flux_o)*coef_D_BT*prob_h_BT_BT );
      trans_ref(22, HORIZONTAL,  BT, EAU, EAU,  BT, (lambda_D+flux_e)*coef_D_BT*prob_h_BT_BT );
      trans_ref(23, HORIZONTAL, EAU,  BT,  BOUE, EAU, (lambda_D+flux_o)*coef_D_BT*prob_h_BT_BOUE );
      trans_ref(24, HORIZONTAL,  BT, EAU, EAU,  BOUE, (lambda_D+flux_e)*coef_D_BT*prob_h_BT_BOUE );
    }
    else{
      trans_ref(21, EST_OUEST, EAU,  BT,  BT, EAU, (lambda_D+flux_o)*coef_D_BT*prob_h_BT_BT );
      trans_ref(22, EST_OUEST,  BT, EAU, EAU,  BT, (lambda_D+flux_e)*coef_D_BT*prob_h_BT_BT );
      trans_ref(23, EST_OUEST, EAU,  BT,  BOUE, EAU, (lambda_D+flux_o)*coef_D_BT*prob_h_BT_BOUE );
      trans_ref(24, EST_OUEST,  BT, EAU, EAU,  BOUE, (lambda_D+flux_e)*coef_D_BT*prob_h_BT_BOUE );
      trans(NORD_SUD, EAU,  BT,  BT, EAU, (lambda_D+flux_n)*coef_D_BT*prob_h_BT_BT );
      trans(NORD_SUD,  BT, EAU, EAU,  BT, (lambda_D+flux_s)*coef_D_BT*prob_h_BT_BT );
      trans(NORD_SUD, EAU,  BT,  BOUE, EAU, (lambda_D+flux_n)*coef_D_BT*prob_h_BT_BOUE );
      trans(NORD_SUD,  BT, EAU, EAU,  BOUE, (lambda_D+flux_s)*coef_D_BT*prob_h_BT_BOUE );
    }

    /// deposition
    trans(  VERTICAL,  TERRE,  BT,  TERRE,  TERRE, lambda_C*coef_C_BT);
    trans(  VERTICAL,  BT,  TERRE,  TERRE,  TERRE, lambda_C*coef_C_BT);
    trans(  VERTICAL,  TERRE,  BT,  TERRE, EAU, lambda_C*2*coef_C_BT);
    trans(  VERTICAL,  BT,  TERRE, EAU,  TERRE, lambda_C*2*coef_C_BT);
    trans(  VERTICAL,  BT,  BT,  BT,  TERRE, lambda_C*coef_C_BT*10 ); //test deposition BT+BT
    trans(  VERTICAL,  BT,  BT,  EAU,  BT, lambda_C*2*coef_C_BT*10 ); //test deposition BT+BT

    /// erosion-deposition (test ...)
    //trans(VERTICAL,  BT, TERRE,  BOUE,  BOUE, lambda_D*lambda_E_v*2 );
    //trans(VERTICAL,  BT, TERRE,  TERRE,  TERRE, lambda_D*lambda_E_v );
  }


  /***** boundary transitions *****/
  if (boundary == BC_OPEN){
    //trans(HORIZONTAL,  BORD,  BOUE, BORD,  EAU, (!flux2)?lambda_D:0 );
    //trans(HORIZONTAL,  BOUE,  BORD, EAU,  BORD, (!flux2)?lambda_D:flux2 /*/2*/ );
    if (flux_eo == flux_ns){
      trans(HORIZONTAL,  BORD,  BOUE, BORD,  EAU, (!flux_o)?lambda_D*gravity:flux_o );
      trans(HORIZONTAL,  BOUE,  BORD, EAU,  BORD, (!flux_e)?lambda_D*gravity:flux_e /*/2*/ );
    }
    else{
      trans(EST_OUEST,  BORD,  BOUE, BORD,  EAU, (!flux_o)?lambda_D:flux_o );
      trans(EST_OUEST,  BOUE,  BORD, EAU,  BORD, (!flux_e)?lambda_D:flux_e /*/2*/ );
      trans(NORD_SUD,  BORD,  BOUE, BORD,  EAU, (!flux_n)?lambda_D:flux_n );
      trans(NORD_SUD,  BOUE,  BORD, EAU,  BORD, (!flux_s)?lambda_D:flux_s /*/2*/ );
    }
    //trans(VERTICAL,  BOUE,  BORD, EAU,  BORD, lambda_D*gravite );
    //trans(VERTICAL,  BOUE,  BORD, DUM,  BORD, lambda_D*gravite );
    //trans(VERTICAL,  BOUE,  BORD, TERRE,  BORD, lambda_C );

    //trans(HORIZONTAL, TERRE, BORD, EAU, BORD, 1); //avalanches aux bords
    //trans(HORIZONTAL, BORD, TERRE, BORD, EAU, 1);
  }

  /***** injection *****/
  trans(  HORIZONTAL,  IN,  EAU,  IN, BOUE, lambda_I );
  //trans(  HORIZONTAL,  IN,  TERRE,  IN, BOUE, lambda_I );  //pour ne pas boucher la source

  /***** gravity *****/
  trans(  VERTICAL, TERRE, EAU, EAU, TERRE, lambda_D*gravity);
  trans(  VERTICAL, TERRE, BOUE, BOUE, TERRE, lambda_D*gravity);
  if (trans_BT){
    trans(  VERTICAL, TERRE, BT, BT, TERRE, lambda_D*gravity);
  }
  if (Ncel[PIERRE] > 0){
    trans(  VERTICAL, PIERRE, EAU, EAU, PIERRE, 10);
    trans(  VERTICAL, PIERRE, BOUE, BOUE, PIERRE, 1);
  }

#ifdef AVALANCHES
  /***** avalanches *****/
  if (ava_trans){
    trans_ref(14, HORIZONTAL, TERRE, EAU, EAU, TERRE, lambda_A);
    trans_ref(15, HORIZONTAL, EAU, TERRE, TERRE, EAU, lambda_A);
    trans_check(14, check_ava, 1);
    trans_check(15, check_ava, 2);
    if (ava_delay){
      trans_regul(14, callback_ava_trans);
      trans_regul(15, callback_ava_trans);
    }
  }
#endif

  /***** chains of transitions for erosion *****/

  /// erosion par transport horizontal, liens verticaux
  if (flux_e == 0){
    trans_link(1, 5, 2, lambda_E_v);
    trans_link(1, 6, 2, lambda_E_v);
  }
  if (flux_o == 0){
    trans_link(2, 5, 1, lambda_E_v);
    trans_link(2, 6, 1, lambda_E_v);
  }

  /// erosion par transport horizontal, liens horizontaux
  if (flux_e == 0){
    trans_link(1, 7, 2, lambda_E_h);
    trans_link(1, 8, 2, lambda_E_h);
    trans_link(1, 9, 2, lambda_E_h);
    trans_link(1, 10, 2, lambda_E_h);
  }
  if (flux_o == 0){
    trans_link(2, 7, 1, lambda_E_h);
    trans_link(2, 8, 1, lambda_E_h);
    trans_link(2, 9, 1, lambda_E_h);
    trans_link(2, 10, 1, lambda_E_h);
  }

  /// erosion par transport vertical, liens horizontaux
  /*trans_link(0, 7, lambda_E_h);
  trans_link(0, 8, lambda_E_h);
  trans_link(0, 9, lambda_E_h);
  trans_link(0, 10, lambda_E_h);*/

  if (trans_BT){
    /// erosion par transport horizontal de cellules BT, liens verticaux
    //trans_link(21, 5, lambda_E_v*5);
    //trans_link(21, 6, lambda_E_v*5);
    //trans_link(22, 5, lambda_E_v*5);
    //trans_link(22, 6, lambda_E_v*5);
    //trans_link(23, 5, lambda_E_v*5);
    //trans_link(23, 6, lambda_E_v*5);
    //trans_link(24, 5, lambda_E_v*5);
    //trans_link(24, 6, lambda_E_v*5);

    //trans_link(23, 5, lambda_E_v);
    //trans_link(23, 6, lambda_E_v);
    //trans_link(24, 5, lambda_E_v);
    //trans_link(24, 6, lambda_E_v);

    trans_link(23, 5, 2, 0.5);
    trans_link(23, 6, 2, 0.5);
    trans_link(24, 5, 1, 0.5);
    trans_link(24, 6, 1, 0.5);

    /// erosion par transport horizontal de cellules BT, liens horizontaux
    /*trans_link(23, 7, 0.05);
    trans_link(23, 8, 0.05);
    trans_link(23, 9, 0.05);
    trans_link(23, 10, 0.05);
    trans_link(24, 7, 0.05);
    trans_link(24, 8, 0.05);
    trans_link(24, 9, 0.05);
    trans_link(24, 10, 0.05);*/

    /// erosion par transport vertical, liens horizontaux
    /*trans_link(20, 7, lambda_E_h);
    trans_link(20, 8, lambda_E_h);
    trans_link(20, 9, lambda_E_h);
    trans_link(20, 10, lambda_E_h);*/ //test
  }
#endif


/*****************************************************************************/
/********************************* CRY model *********************************/
/*****************************************************************************/
#ifdef MODEL_CRY

  /***** crystallization *****/
  trans(VERTICAL, PLUS, ZERO, PLUS, PLUS, lambda_C/100000);
  trans(VERTICAL, ZERO, PLUS, PLUS, PLUS, lambda_C/100000);
  trans(HORIZONTAL, PLUS, ZERO, PLUS, PLUS, lambda_C);
  trans(HORIZONTAL, ZERO, PLUS, PLUS, PLUS, lambda_C);

  /***** diffusion *****/
  trans(VERTICAL, ZERO, AIR, AIR, ZERO, lambda_D);
  trans(VERTICAL, AIR, ZERO, ZERO, AIR, lambda_D);
  trans(HORIZONTAL, ZERO, AIR, AIR, ZERO, lambda_D+flux2);
  trans(HORIZONTAL, AIR, ZERO, ZERO, AIR, lambda_D);

  /***** injection *****/
  trans(VERTICAL, IN, AIR, IN, ZERO, lambda_D/100);

  /***** boundary transitions *****/
  if (boundary == BC_OPEN){
    trans(HORIZONTAL, ZERO, BORD, AIR, BORD, lambda_D+flux2);
  }

#endif


/*****************************************************************************/
/********************************* DIF model *********************************/
/*****************************************************************************/
#ifdef MODEL_DIF

  /***** diffusion *****/
  trans(ISOTROPE, ZERO, ONE, ONE, ZERO, lambda_D);

  /***** injection *****/
  trans(VERTICAL, IN, ZERO, IN, ONE, lambda_I);

  /***** boundary transitions *****/
  if (boundary == BC_OPEN){
    trans(HORIZONTAL, ONE, BORD, ZERO, BORD, lambda_D);
    trans(HORIZONTAL, BORD, ONE, BORD, ZERO, lambda_D);
  }

#endif

/*****************************************************************************/
/********************************* D2G model *********************************/
/*****************************************************************************/
#ifdef MODEL_D2G

  /***** diffusion *****/
  trans(VERTICAL, EAU, AIR, AIR, EAU, 1000);
  trans(HORIZONTAL, EAU, AIR, AIR, EAU, 1);
  trans(HORIZONTAL, AIR, EAU, EAU, AIR, 1);

  /***** dissolution *****/
  //trans(VERTICAL, EAU, TERRE, EAU, BOUE, 1000);
  trans(VERTICAL, EAU, TERRE, AIR, BOUE, lambda_L);
  trans(VERTICAL, EAU, BOUE, AIR, BOUE, lambda_L);

  /***** injection *****/
  trans(VERTICAL, IN, AIR, IN, EAU, lambda_I);

  /***** boundary transitions *****/
  if (boundary == BC_OPEN){
    trans(HORIZONTAL, EAU, BORD, AIR, BORD, 1);
    trans(HORIZONTAL, BORD, EAU, BORD, AIR, 1);
  }

#endif

  //LogPrintf("gravite = %f\n",gravite);

/*****************************************************************************/
/********************************* LIFE model ********************************/
/*****************************************************************************/
#ifdef MODEL_LIFE

  /***** birth *****/
  //trans_ref(0,  ISOTROPE,  ALIVE, DEAD,  ALIVE,  ALIVE, lambda_B );
  trans_ref(0,  HORIZONTAL,  ALIVE, DEAD,  ALIVE,  ALIVE, lambda_B );
  trans_ref(1,  HORIZONTAL,  DEAD,  ALIVE, ALIVE,  ALIVE, lambda_B );
  trans_check(0, check_birth, 2);
  trans_check(1, check_birth, 1);
  trans_time(0, TIME_EVOL);
  trans_time(1, TIME_EVOL);

  /***** death *****/
  //trans_ref(1,  ISOTROPE,  ALIVE, DEAD,  DEAD,  DEAD, lambda_D );
  //trans_ref(2,  ISOTROPE,  ALIVE, ALIVE,  ALIVE,  DEAD, lambda_D );
  trans_ref(2,  HORIZONTAL,  ALIVE, DEAD,  DEAD,  DEAD, lambda_D );
  trans_ref(3,  HORIZONTAL,  DEAD,  ALIVE, DEAD,  DEAD, lambda_D );
  trans_ref(4,  HORIZONTAL,  ALIVE, ALIVE, ALIVE, DEAD, lambda_D );
  trans_ref(5,  HORIZONTAL,  ALIVE, ALIVE, DEAD,  ALIVE,lambda_D );
  trans_check(2, check_death, 1);
  trans_check(3, check_death, 2);
  trans_check(4, check_death, 2);
  trans_check(5, check_death, 1);
  trans_time(2, TIME_EVOL);
  trans_time(3, TIME_EVOL);
  trans_time(4, TIME_EVOL);
  trans_time(5, TIME_EVOL);
#endif
}


#ifdef MODEL_DUN

float coef_injection = 1.0; //regulation de l'injection de grain

int callback_bord_dun(void *data)
{
  static int Ncel_GR_debut = 0;
  int Ncel_GR;
  double *ref_intensite = (double*)data;

  Ncel_GR = Ncel[GR] + Ncel[GRJ];

  if (!Ncel_GR_debut){
    Ncel_GR_debut = Ncel_GR;
    coef_injection = *ref_intensite;
    *ref_intensite = 0.0;
  }
  else{
    //si la quantite de grain diminue, on reinjecte de la matiere
    //*ref_intensite = (Ncel_GR < Ncel_GR_debut) ? 1.0 : 0.0;
    *ref_intensite = (Ncel_GR < Ncel_GR_debut) ? coef_injection : 0.0;
    //if (*ref_intensite) LogPrintf("intensite injection = %f\n",coef_injection);
  }

  return 0;
}

int callback_injection_coef(void *data)
{
  static char start = 1;
  static double lambda_injection_debut;

  double *ref_intensite = (double*)data;

  if (start){
    lambda_injection_debut = *ref_intensite;
    start = 0;
  }
  *ref_intensite = lambda_injection_debut*coef_injection;
  //LogPrintf("intensite injection = %f\n", *ref_intensite);

  return 0;
}

#endif  //MODEL_DUN

#ifdef AVALANCHES

//modulation intensite des transitions d'avalanches
int callback_ava_trans(void *data)
{
  double *ref_intensite = (double*)data;

  *ref_intensite = lambda_A;

  return 0;
}

#endif

