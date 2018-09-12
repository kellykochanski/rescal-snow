/* ReSCAL - Global definitions
 *
 * Copyright (C) 2011-2016
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

//#ifdef HAVE_CONFIG_H
//#  include <config.h>
//#endif

/// version number and date
#define VER_NUM "1.6"
#define VER_DAT "2016-10"

/// numeric format
#define NUM_MODE "C"  //C language
//#define NUM_MODE "en_US.UTF-8"  //english
//#define NUM_MODE "fr_FR.UTF-8"  //french

/// ======================================================
/// programs
enum PROGS {PROG_RESCAL, PROG_TOOL};

/// ======================================================

/// choice of a model
//#define MODEL_DUN    //Dune morphodynamics
#define MODEL_SNO    //Snow bedform morphodynamics

/// ======================================================

/// common options
// #define OPENMP   //use OpenMP parallelization (on lattice gas)
// #define GUI    //graphical interface
#define CYCLAGE_HOR     //enable horizontal cycling
#define INFO_CEL        //log the number of cells
#define INFO_DBL        //log the number of active doublets
#define INFO_TRANS      //log the number of effective transitions
#define DETERMINISTIC   //try to keep a deterministic execution: same parameters produce same results (it may be slower when DETERMINISTIC is defined, but not always !)
#define DUMP_SIGNATURE  //footprints for the comparison of two simulations (identical simulations => identical footprints)
//#define LOG_FILE        //ReSCAL logs redirected into a file
#define USE_LIBPNG        //use of libpng
//#define USE_GD          //use of libgd (alternative to USE_LIBPNG, but obsolete)
//#define ROTATE_LIGHT    //automatic rotation of the light source in ROT_CYCLE mode
#define REORIENT_AUTO   //automatic reorientation when saving (images of) the cellular space
#define CSP_MUTEX       //mutual exclusion on the cellular space (no transition while saving data)
//#define GTK_OLD_SYNTAX  //compile with old version of GTK+2.x (obsolete)
// #define PARALLEL_AUTOMATA //enable concurrent execution of stochastic engine and lattice gas  (non-deterministic only, not working on MacOS X)



/// ===================== SNO MODEL =====================
// KK 10 May 2018
// Most cell types are inherited from the DUN MODEL (below)
#ifdef MODEL_SNO
/// name of the model
#define MOD_NAME "SNO"
#define MOD_DESC "simulation of snow surface under air flow"

#endif

/// ===================== DUN MODEL =====================
//  Most properties of dunes - except name and vegetation - are inherited by SNO MODEL
#ifdef MODEL_DUN

/// name of the model
#define MOD_NAME "DUN"

/// description
#define MOD_DESC "simulation of sand dunes under fluid flow"
#endif // MODEL_DUN

#if defined(MODEL_DUN) || defined(MODEL_SNO)
/// cell states
#define GR    0 //grain
#define GRJ   1 //mobile grain
//#define BR    2 //bedrock (no transitions)
#define GRV    2 //vegetated grain
#define EAUC  3 //air (or water)
#define VEG   4 //vegetation
#define BORD  5 //boundary
#define DUM   6 //neutral
#define IN    7 //input of sand
#define OUT   8 //output of sand
#define TUNNEL 9 //tunnel to another cellular space (with PARALLEL option)
#define EAUT  10 //(used only graphically, no transitions)
#define GRC   11 //colored grain (used only graphically, no transitions)

/// names of cell states
//#define ETATS { "GR", "GRJ", "BR", "EAUC", "VEG", "BORD", "DUM", "IN", "OUT", "TUNNEL" }
#define ETATS { "GR", "GRJ", "GRV", "EAUC", "VEG", "BORD", "DUM", "IN", "OUT", "TUNNEL" }

/// phases of cell states
#define PHASES { SOLID, SOLID, SOLID, FLUID, SOLID, SOLID, SOLID, SOLID, SOLID, SOLID }
//#define PHASES { SOLID, FLUID, SOLID, FLUID, SOLID, SOLID, SOLID, SOLID, SOLID, SOLID }

/// specific options of the model
#define ALTI GR //compute elevation map for cells of type GR
#define AVALANCHES
#define ROTATIONS
#define VARIABLE_FLOW // KK - vary the intensity of the flow
#define CENTERING_AUTO //automatic centering
#define LIENS_TRANSITIONS //correlations of transitions
#define LGCA //lattice gas
#define CGV //compute shear stress (with LGCA option)
#define WIND_DATA //option for importing wind data from a file (with ROTATIONS option)
#define TIME_SCALE //option for computing the time scale (using additional physical parameters)
#define CELL_TIME GR //save time for each cell of type GR
//#define NORM2D //CGV with 2d normals
//#define TRACE3D_CEL
//#define TRACE_FLUX //save the mean flux of grains (2d)
//#define STABILITY_ANALYSIS //produce data for the 2D linear stability analysis (find the length scale of the model)
//#define APPEL_MCC //external matlab application

#endif // DUN or SNO


/// verification of compatibility between the options

#ifndef ROTATIONS
#undef REORIENT_AUTO
#endif

#ifndef CGV
#undef TIME_SCALE
#endif

#ifndef CELL_COLOR
#undef TRACE_PAR_COL
#endif

#ifdef STABILITY_ANALYSIS
#define NORM2D
#define DUMP_SIGMA
#endif

#ifndef LGCA
#undef PARALLEL_AUTOMATA
#endif

//#ifdef PARALLEL_AUTOMATA
//#undef DETERMINISTIC
//#endif

/// ======================================================

/// common constants

#define MAX_CELL  12  //max number of cell types
#define MAX_DB    100 //max number of doublets in the model
#define MAX_CHK   3   //max number of check functions for one transition

#define EST       0
#define BAS       1
#define SUD       2

#define OUEST     3
#define HAUT      4
#define NORD      5

#define VERTICAL   0
#define EST_OUEST  1
#define NORD_SUD   2
#define HORIZONTAL 3 // EST_OUEST ou NORD_SUD
#define ISOTROPE   4
#define CLASSES_DB { "VERTICAL", "EAST_WEST", "NORTH_SOUTH", "HORIZONTAL" ,"ISOTROPIC"}

#define MAX_TRANSITIONS_DB 100
#define MAX_TRANSITIONS_CEL 10
#define MAX_LIENS  100

#define MAX_TRANSITIONS MAX_TRANSITIONS_DB+MAX_TRANSITIONS_CEL

#ifdef PARALLEL
#define MAX_NODE  16
#endif

#define PI 3.14159265

#define W32 1
#define W64 2

#define SOLID 0
#define FLUID 1

/// ======================================================

