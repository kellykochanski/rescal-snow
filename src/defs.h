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
//#define MODEL_AVA    //Avalanches
//#define MODEL_CMB    //Core-mantle boundary
//#define MODEL_ICB    //Innercore-core boundary
//#define MODEL_RIV    //River (landscape evolution)
//#define MODEL_CRY    //Crystallization (of calcite)
//#define MODEL_DIF    //Diffusion
//#define MODEL_D2G    //Diffusion 2d with gravity
//#define MODEL_LIFE   //Stochastic game of life (2D)
#define MODEL_SNO    //Snow bedform morphodynamics

/// ======================================================

/// common options
#define OPENMP		//use OpenMP parallelization (on lattice gas)
#define GUI		//graphical interface
#define CYCLAGE_HOR     //enable horizontal cycling
#define INFO_CEL        //log the number of cells
#define INFO_DBL        //log the number of active doublets
#define INFO_TRANS      //log the number of effective transitions
#define DETERMINISTIC	  //try to keep a deterministic execution: same parameters produce same results (it may be slower when DETERMINISTIC is defined, but not always !)
#define DUMP_SIGNATURE  //footprints for the comparison of two simulations (identical simulations => identical footprints)
//#define LOG_FILE        //ReSCAL logs redirected into a file
#define USE_LIBPNG        //use of libpng
//#define USE_GD          //use of libgd (alternative to USE_LIBPNG, but obsolete)
//#define ROTATE_LIGHT    //automatic rotation of the light source in ROT_CYCLE mode
#define REORIENT_AUTO   //automatic reorientation when saving (images of) the cellular space
#define CSP_MUTEX       //mutual exclusion on the cellular space (no transition while saving data)
//#define GTK_OLD_SYNTAX  //compile with old version of GTK+2.x (obsolete)
#define PARALLEL_AUTOMATA //enable concurrent execution of stochastic engine and lattice gas  (non-deterministic only, not working on MacOS X)

#if !defined(_OPENMP) || defined(__CYGWIN32__)
#undef OPENMP
#endif

#if !defined(USE_AVX) && defined(__AVX__)
#define USE_AVX
#endif

/// ===================== CMB MODEL =====================
#ifdef MODEL_CMB

/// name of the model
#define MOD_NAME "CMB"

/// description
#define MOD_DESC "simulation of the core-mantle boundary"

/// cell states
#define PLUS  0 //solide
#define ZERO  1 //liquide sature en elements legers
#define MOINS 2 //liquide non-sature
#define BORD  3 //bord
#define DUM   4 //inerte
#define PIERRE 7 //granit
#define TUNNEL 9 //tunnel vers une terre parallele

/// names of cell states
#define ETATS { "PLUS", "ZERO", "MOINS", "BORD", "DUM", "", "", "GRANIT", "", "TUNNEL" }

/// phases of cell states
#define PHASES { SOLID, FLUID, FLUID, SOLID, SOLID, SOLID, SOLID, SOLID, SOLID, SOLID }

/// specific options of the model
#define DUMP_RUGOSI
//#define ROTATIONS
#endif

/// ===================== SNO MODEL =====================
// KK 10 May 2018
#ifdef MODEL_SNO
/// name of the model
#define MOD_NAME "SNO"
#define MOD_DESC "simulation of snow surface under fluid flow"

/// cell states
#define GR    0 //grain, no cohesion
#define GRJ   1 //grain, mobile
#define GRH   2 //grain, cohesive/hardened
#define EAUC  3 //air
#define VEG   4 //vegetation
#define BORD  5 //boundary
#define DUM   6 //neutral (e.g. bedrock)
#define IN    7 //input of snow
#define OUT   8 // output of snow
#define TUNNEL  9 // tunnel to another cellular space (with PARALLEL option)
#define EAUT    10 // (used graphically, no transitions)
#define GRC   11 // grain, no cohesion (used only graphically)

//names of cell states
#define ETATS {"GR", "GRJ", "GRH", "EAUC", "VEG", "BORD", "DUM", "IN", "OUT", "TUNNEL"}

// phases of cell states
#define PHASES { SOLID, SOLID, SOLID, FLUID, SOLID, SOLID, SOLID, SOLID, SOLID, SOLID }

/// specific options of the model
#define ALTI GR //compute elevation map for cells of type GR
#define AVALANCHES
#define ROTATIONS
#define CENTERING_AUTO //automatic centering
#define LIENS_TRANSITIONS //correlations of transitions
#define LGCA //lattice gas
#define CGV //compute shear stress (with LGCA option)
#define WIND_DATA //option for importing wind data from a file (with ROTATIONS option)
#define TIME_SCALE //option for computing the time scale (using additional physical parameters)
#define CELL_TIME GR //save time for each cell of type GR

#endif

/// ===================== DUN MODEL =====================
#ifdef MODEL_DUN

/// name of the model
#define MOD_NAME "DUN"

/// description
#define MOD_DESC "simulation of sand dunes under fluid flow"

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
#define CENTERING_AUTO //automatic centering
#define LIENS_TRANSITIONS //correlations of transitions
#define LGCA //lattice gas
#define CGV //compute shear stress (with LGCA option)
#define WIND_DATA //option for importing wind data from a file (with ROTATIONS option)
#define TIME_SCALE //option for computing the time scale (using additional physical parameters)
//#define CELL_TIME GR //save time for each cell of type GR
//#define CELL_COLOR //use colored cells
//#define NORM2D //CGV with 2d normals
//#define TRACE_TRANS
//#define TRACE3D_CEL
//#define TRACE_PAR_COL //save the courses of colored cells (one cell per color)
//#define TRACE_FLUX //save the mean flux of grains (2d)
//#define STABILITY_ANALYSIS //produce data for the 2D linear stability analysis (find the length scale of the model)
//#define DUMP_SIGMA //computation of the standard deviation in height (when using WAVES_2D template)
//#define DUMP_AUTOCORREL //computation of the auto-correlation in height (when using WAVES_2D template)
//#define APPEL_MCC //external matlab application
//#define USE_VEGETATION //option for using cells in VEG state
#endif

/// ===================== AVA MODEL =====================
#ifdef MODEL_AVA

/// name of the model
#define MOD_NAME "AVA"

/// description
#define MOD_DESC "simulation of granular avalanches"

/// cell states
#define GR    0 //grain stable
#define AIR   1 //air
#define BORD  2 //bord
#define DUM   3 //inerte
#define IN    4 //source de grains
#define OUT   5 //sortie de grains
#define TUNNEL 9 //tunnel vers une terre parallele

#define TERRE    GR

/// names of cell states
#define ETATS { "GR", "AIR", "BORD", "DUM", "IN", "OUT", "", "", "", "TUNNEL" }

/// phases of cell states
#define PHASES { SOLID, FLUID, SOLID, SOLID, SOLID, SOLID }

/// specific options of the model
#define ALTI GR //calcul de la hauteur des piles de cellules GR
#define AVALANCHES  //avalanches
//#define ROTATIONS //table tournante
//#define CELL_TIME GR //datation des cellules GR
#define CELL_COLOR //cellules colorees
//#define NORM2D //pour l'affichage
//#define TRACE_TRANS
#endif

/// ===================== ICB MODEL =====================
#ifdef MODEL_ICB

/// name of the model
#define MOD_NAME "ICB"

/// description
#define MOD_DESC "simulation of the Earth's inner-core boundary"

/// cell states
#define PLUS  0 //fer solide
#define ZERO  1 //liquide sature en elements legers
#define MOINS 2 //liquide non-sature
#define BORD  3 //bord
#define DUM   4 //inerte
#define TUNNEL 9 //tunnel vers une terre parallele

/// names of cell states
#define ETATS { "PLUS", "ZERO", "MOINS", "BORD" , "DUM"}

/// phases of cell states
#define PHASES { SOLID, FLUID, FLUID, SOLID, SOLID }

/// specific options of the model
#define DUMP_RUGOSI

#endif

/// ===================== RIV MODEL =====================
#ifdef MODEL_RIV

/// name of the model
#define MOD_NAME "RIV"

/// description
#define MOD_DESC "simulation of erosion by runoff"

/// cell states
#define TERRE  0 //terre
#define BOUE   1 //boue
#define EAU    2 //eau
#define BT     3 //boue turbulent
#define PIERRE 4 //pierre
#define BORD   5 //bord
#define DUM    6 //inerte
#define IN     7 //source
#define OUT    8 //sortie
#define TUNNEL 9 //tunnel vers une terre parallele

/// names of cell states
#define ETATS { "TERRE", "BOUE", "EAU", "BT", "PIERRE", "BORD", "DUM", "IN", "OUT", "TUNNEL"}

/// phases of cell states
#define PHASES { SOLID, FLUID, FLUID, FLUID, SOLID, SOLID, SOLID, SOLID, SOLID, SOLID }

/// specific options of the model
#define LIENS_TRANSITIONS
#define ALTI TERRE //calcul de la hauteur des piles de cellules TERRE
#define AVALANCHES
#define SURRECTIONS
//#define TRACE_SRC //calcul bassin versant
//#define TRACE_AIRE  //calcul aires drainees
//#define TRACE_TRANS
//#define TRACE_PAR //calcul debit + pente
//#define APPEL_GMT
//#define LISSAGE  //couche uniforme de BOUE + cristallisation BOUE->TERRE + pas de transport vers le haut => lissage avant mode TRACE
//#define CELL_TIME TERRE //datation des cellules TERRE


#endif


/// ===================== CRYSTAL MODELE =====================
#ifdef MODEL_CRY

/// name of the model
#define MOD_NAME "CRYSTAL"

/// description
#define MOD_DESC "simulation of calcite crystal growth"

/// cell states
#define PLUS  0 //solide
#define ZERO  1 //gaz
#define BORD  3 //bord
#define DUM   4 //inerte
#define IN    5 //source du gaz
#define AIR   6 //air
#define TUNNEL 9 //tunnel vers une terre parallele

/// names of cell states
#define ETATS { "SOLIDE", "GAZ", "", "BORD", "DUM", "SOURCE", "AIR" }

/// phases of cell states
#define PHASES { SOLID, FLUID, SOLID, SOLID, SOLID, SOLID, FLUID }

#endif

/// ===================== DIF MODEL =====================
#ifdef MODEL_DIF

/// name of the model
#define MOD_NAME "DIF"

/// description
#define MOD_DESC "isotropic diffusion"

/// cell states
#define ZERO  0 //gaz 0
#define ONE   1 //gaz 1
#define IN    2 //source
#define BORD  3 //bord
#define DUM   4 //inerte


/// names of cell states
#define ETATS { "ZERO", "ONE", "IN", "BORD", "DUM" }

/// phases of cell states
#define PHASES { FLUID, FLUID, SOLID, SOLID, SOLID }

#endif

/// ===================== D2G MODEL =====================
#ifdef MODEL_D2G

/// name of the model
#define MOD_NAME "D2G"

/// description
#define MOD_DESC "diffusion 2d with gravity"

/// cell states
#define TERRE 0 //terre
#define EAU  1 //eau
#define AIR   2 //air
#define BORD  3 //bord
#define DUM   4 //inerte
#define IN    5 //source
#define BOUE  6 //boue

/// names of cell states
#define ETATS { "TERRE", "EAU", "AIR", "BORD", "DUM", "SOURCE", "BOUE" }

/// phases of cell states
#define PHASES { SOLID, FLUID, FLUID, SOLID, SOLID, SOLID, FLUID }

#endif

/// ===================== LIFE MODEL =====================
#ifdef MODEL_LIFE

/// name of the model
#define MOD_NAME "LIFE"

/// description
#define MOD_DESC "stochastic game of life (2D)"

/// cell states
#define DEAD 0
#define ALIVE 1
#define BORD  2 //bord
#define DUM   3 //inerte

/// names of cell states
#define ETATS { "DEAD", "ALIVE", "BORD", "DUM" }

#endif

/// ======================================================

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

