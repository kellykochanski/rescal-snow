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
 * aint64_t with this program; if not, write to the Free Software
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
#include <stdint.h>

#include "defs.h"
#include "macros.h"
#include "param.h"
#include "transitions.h"
#include "space.h"
#include "lgca.h"

extern int32_t Ncel[];    //number of cells in each state
extern int32_t LNS;       //north-south width
extern int32_t ava_trans;   //flag for avalanche transitions
extern float ava_delay;
extern int32_t boundary;  //boundary conditions

/// parameters and rates for each model

char *model_name = NULL; //name of the model
#ifdef LGCA
int32_t use_lgca = 1;
#else
int32_t use_lgca = 0;
#endif // LGCA

#if defined(MODEL_SNO) || defined(MODEL_DUN)
float lambda_E, lambda_C, lambda_D, lambda_G, lambda_J;
float lambda_D_mob;
float lambda_T = -1;
float lambda_I = 1.0;
float coef_a = 0.0;
float coef_b = 10.0;
float coef_c = 500.0;
float prob_link_ET = 0.5;
float prob_link_TT = 0.5;
int32_t flag_hm = 0;
#endif // DUN or SNO


#ifdef MODEL_SNO // options for sno but not dun model
// how much harder is it to erode GRV than GR grains
float lambda_F = 0.1;
// how fast do grains sinter into place
float lambda_S;
#endif // SNO

#ifdef MODEL_AVA
float lambda_I = 0.1;
#endif

#ifdef AVALANCHES
float lambda_A = 1.0;
float lambda_A_unstable = 0;
#endif

float gravity = 0;


/// callbacks for the regulation of transitions
#if defined(MODEL_DUN) || defined(MODEL_SNO)
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
#ifdef LGCA
// Callback_check check_mvt_bas;
// Callback_check check_mvt_haut;
// Callback_check check_mvt_est;
// Callback_check check_mvt_ouest;
// Callback_check check_mvt_est_et_bas;
// Callback_check check_no_mvt_est_et_bas;
// Callback_check check_mvt_mean_bas;
// Callback_check check_mvt_mean_est;
// Callback_check check_no_mvt_mean_bas;
// Callback_check check_mvt_EB;
Callback_check check_grad_vel;
#endif
Callback_check check_cell_top;


void params_modele() {
  char str[100];

  /// declare families of parameters
  sprintf(str, "%s model parameters", MOD_NAME);
  param_family("MODEL", str);

  /// declare parameters
#if defined(MODEL_DUN) || defined(MODEL_SNO)
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
#endif // sno or dun

#ifdef MODEL_SNO
  parameter("Lambda_F", "Factor decrease in erosion rate if particle is sintered (0.1 by default)", &lambda_F, PARAM_FLOAT, "MODEL");
  parameter("Lambda_S", "Sintering rate (turns loose grains into cohesive)", &lambda_S, PARAM_FLOAT, "MODEL");
#endif //sno


}

void init_modele()
{
  if (model_name && strcmp(model_name, MOD_NAME)){
    ErrPrintf("ERROR: Wrong model name %s in parameter file. Check model in defs.h and rebuild all.\n", model_name);
    exit(-2);
  }

  gravity += 1.0; //= 1.05;



/*****************************************************************************/
/******************************** SNO model **********************************/
/*****************************************************************************/
// KK 08/May/2018
// Snow also inherits many transitions from DUN model (below)
#ifdef MODEL_SNO

  /************cohesion******************/
/* Parameters:
 * lambda_F   ratio of erosion of loose to sintered snow (<=1)
 *              lambda_F = 0 -> sintered grains are not erodible at all
 * lambda_S   rate of sintering (turns loose snow into sintered snow)
*/
  double vert_erosion = 100000.0; // ratio of horizontal:vertical erosion
  				 // value taken from sand dun

  /****** erosion of cohesive (sintered) grains ******/
  // grains turn from GRV (cohesive) into GRJ (moving)
  // when they land, they will be loose (GR), their cohesive bonds broken
  if ( lambda_F > 0){
  	trans_ref(101, EST_OUEST, EAUC, GRV, EAUC, GRJ, lambda_E*lambda_F );
  	trans_ref(102, VERTICAL,  EAUC, GRV, GRJ, EAUC, lambda_E*lambda_F/vert_erosion );
  }

  // KK -- loose grains (GR) to turn cohesive (GRV)
  // This requires a grain to transition in place, does not appear to be implemented yet
  // Loose grains (GR) sitting on solids (DUM, BORD, GRV, sinter into place at rate lambda_S:
  trans_ref(103, VERTICAL,  GR, DUM,  GRV, DUM,  lambda_S);
  trans_ref(104, VERTICAL,  GR, BORD, GRV, BORD, lambda_S);
  trans_ref(104, VERTICAL,  GR, GRV,  GRV, GRV,  lambda_S);

#endif //sno

/*****************************************************************************/
/********************************* DUN model *********************************/
/*****************************************************************************/
#if defined(MODEL_SNO) || defined(MODEL_DUN)

  /* Parameters:
    lambda_E    erosion
    lambda_C    deposition
    lambda_T    transport
    lambda_D    diffusion (NS)
    lambda_G    gravity
    lambda_I    injection
    lambda_A    avalanches
    lambda_J    not used
    _______________________________________________

    7 lambdas + 1 not used

  */

  char tmode = TIME_CORR; //time correction (default value for check_grad_vel)

  /***** erosion *****/
  trans_ref(8,  EST_OUEST, EAUC, GR,  EAUC, GRJ, lambda_E);
  trans_ref(7,  VERTICAL,  EAUC, GR,  GRJ,  EAUC, lambda_E / 100000.0);
  trans_type(7, TR_TRANSPORT);

#ifdef CGV
  //control of the erosion by shear stress
  Callback_check *pchk = check_grad_vel;
  if (use_lgca) {
    trans_check(8, pchk, 2);
    trans_check(7, pchk, 2);
    trans_time(8, tmode);
    trans_time(7, tmode);
  }
#endif //CGV

  /***** gravity *****/
  trans(VERTICAL, GR, EAUC, EAUC,  GR, lambda_G);
  trans(VERTICAL, GRJ, EAUC, EAUC, GRJ, (coef_a > 0) ? lambda_T / coef_a : lambda_G);

  /***** transport *****/
  if (lambda_T < 0) {
    lambda_T = 3*lambda_C;  //default value of lambda_T (with correlations of transitions)
  }
  trans_ref(10,  VERTICAL,  EAUC, GRJ,  GRJ,  EAUC, coef_a * lambda_T);
  trans_ref(11,  EST_OUEST,  GRJ,  EAUC,  EAUC, GRJ, lambda_T);
  if (flag_hm) {
    trans_ref(50,  EST_OUEST, GRJ, GR, GRJ, GR, lambda_T);
  }

  /***** chains of transitions for saltation *****/
  trans_link(8, 10, 2, prob_link_ET/*0.5*/);
  trans_link(10, 11, 1, prob_link_TT/*0.5*/);
  trans_link(7, 11, 1, 1.0);
  if (flag_hm) {
    trans_link(50, 10, 1, prob_link_TT/*0.5*/);
  }

  /***** deposition *****/
  trans_ref(9,  EST_OUEST, GRJ, EAUC, GR, EAUC, lambda_C);
  trans_ref(16,  EST_OUEST, GRJ, GR, GR, GR, coef_b * lambda_C);
  trans_ref(17,  EST_OUEST, GRJ, DUM, GR, DUM, coef_b * lambda_C);
  trans(EST_OUEST, GR, GRJ, GR, GR, coef_c * lambda_C/*1000*/);

  /***** diffusion *****/
  if (LNS > 1) {
    trans_ref(12, NORD_SUD,   GR, EAUC, EAUC,  GR, lambda_D);
    trans_ref(13, NORD_SUD,   EAUC,  GR, GR, EAUC, lambda_D);
    trans(NORD_SUD,   GRJ, EAUC, EAUC,  GRJ, lambda_D_mob);
    trans(NORD_SUD,   EAUC,  GRJ, GRJ, EAUC, lambda_D_mob);
#ifdef CGV
    if (use_lgca) {
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
    trans_ref(14, HORIZONTAL, GR, EAUC, EAUC, GR, lambda_A);
    trans_ref(15, HORIZONTAL, EAUC, GR, GR, EAUC, lambda_A);
    trans_check(14, check_ava, 1);
    trans_check(15, check_ava, 2);
    if (ava_delay) {
      trans_regul(14, callback_ava_trans);
      trans_regul(15, callback_ava_trans);
    }
  }
#endif

  /***** boundary transitions *****/
  if ((boundary == BC_OUT) || (boundary == BC_REINJECTION)) {
    LogPrintf("transitions sur le bord\n");
    trans(EST_OUEST,  GRJ,  BORD, EAUC, BORD, 20);
    trans(EST_OUEST,  GR,  BORD, EAUC, BORD, 20);
  } else if (boundary == BC_OPEN) {
    /// only limited transport
    LogPrintf("transitions lentes sur le bord\n");
    trans(EST_OUEST,  GRJ,  BORD, EAUC, BORD, 1);
  }

  /***** injection *****/
  //injection of sand grains (for the computation of Qsat and Lsat)
  if (Ncel[IN] > 0) {
    trans_ref(2, VERTICAL, IN,  EAUC, IN, GR, lambda_I/*0.5*/); // value 2 sets the (arbitrary) reference number of the transition
    trans(VERTICAL, EAUC, IN, GRJ, IN, lambda_I/*0.05*/);
  }

  /***** removal *****/
  //output of sand
  if (Ncel[OUT] > 0) {
    trans(VERTICAL, GR, OUT, EAUC, OUT, 1000);
    trans(VERTICAL, GRJ, OUT, EAUC, OUT, 1000);
  }

  /***** reinjection *****/
  if ((boundary == BC_REINJECTION) && (Ncel[IN] > 0)) {
    //reinjection of sand grains
    LogPrintf("reinjection de grains\n");
    trans_regul(2, callback_bord_dun);
  }
#endif //MODEL_DUN
}


#if defined(MODEL_SNO) || defined(MODEL_DUN)
float coef_injection = 1.0; //regulation de l'injection de grain

int32_t callback_bord_dun(void *data) {
  static int32_t Ncel_GR_debut = 0;
  int32_t Ncel_GR;
  double *ref_intensite = (double*)data;

  Ncel_GR = Ncel[GR] + Ncel[GRJ];

  if (!Ncel_GR_debut) {
    Ncel_GR_debut = Ncel_GR;
    coef_injection = *ref_intensite;
    *ref_intensite = 0.0;
  } else {
    //si la quantite de grain diminue, on reinjecte de la matiere
    *ref_intensite = (Ncel_GR < Ncel_GR_debut) ? coef_injection : 0.0;
  }

  return 0;
}

int32_t callback_injection_coef(void *data) {
  static char start = 1;
  static double lambda_injection_debut;

  double *ref_intensite = (double*)data;

  if (start) {
    lambda_injection_debut = *ref_intensite;
    start = 0;
  }
  *ref_intensite = lambda_injection_debut * coef_injection;
  //LogPrintf("intensite injection = %f\n", *ref_intensite);

  return 0;
}

#endif  //MODEL_DUN or MODEL_SNO

#ifdef AVALANCHES

//modulation intensite des transitions d'avalanches
int32_t callback_ava_trans(void *data) {
  double *ref_intensite = (double*)data;

  *ref_intensite = lambda_A;

  return 0;
}

#endif

