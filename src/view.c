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
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#include "defs.h"
#include "macros.h"
#include "view.h"
#include "cells.h"
#include "space.h"
#include "surface.h"
#include "trace.h"
#include "transitions.h"
#include "lgca.h"

extern int prog;    //main prog
extern int H, L, D, HL, HLD;       // les dimensions de la terre
extern int L_bounds, D_bounds;
extern int LEO, LN, LS, LNS, HLN;    //largeur est-ouest, limite nord, limite sud, largeur nord-sud
extern Cell  *TE;	           // la 'terre'
extern const unsigned char Phase[MAX_CELL];	//phase (fluide ou solide) des types de cellules
#ifdef ALTI
extern short *alti; // elevations locales du terrain
extern Vec2 *norm2d;
extern Vec3 *norm3d;          // normale 3d a la surface definie par alti
#endif
extern char *psol;                    // indicateur des plans solides verticaux est-ouest
extern char *csol;
extern int h_ceil;  //epaisseur du plafond
extern int h_floor;
#ifdef LGCA
extern int VH, VL, VHL, CLNS ;	// dimensions du tableau des vitesses
//extern int *Velx, *Vely;	//vitesses locales moyennees
extern int *Velx_sum, *Vely_sum;	//vitesses locales additives
extern float *Velx_interp, *Vely_interp;	//vitesses locales interpolees
extern int use_lgca;
extern int lgca_ready;
#endif
#ifdef USE_VEGETATION
extern int use_veg;
#endif
extern char *rot_map;        // periodic mapping of the rotating space
extern int rot_mode;

unsigned char opt_ch=0;
unsigned char opt_cv=0;
unsigned char opt_lc=0;
unsigned char opt_tr=0;
unsigned char opt_al=0;
unsigned char opt_lal=0;
unsigned char opt_vel=0;
unsigned char opt_vss=0;
unsigned char opt_ls=0;
int hh=0;
int khh=0;
int abs_cv=0, prof_cv=0;
unsigned char reorient_flag=0;
unsigned char horizontal_display=0;

unsigned char visible[MAX_CELL];  //visibilite des cellules par type
#ifdef CELL_COLOR
char colored[MAX_CELL];  //affichage de la couleur par type
#endif
unsigned char *image=NULL;
unsigned char *image_end=NULL;
int nb_shades=1;
int col_shade_index=0;
int img_w, img_h;
Shade shades[NB_SHADE_MAX];
int cel_shade_index[MAX_CELL];  //correspondance celltype --> shade

//direction of the incident light source
float light_elevation = 45; // in degrees, from the horizon
float light_azimuth = -90; //in degrees, clockwise from north
Vec3 light;

//display mode of the current orientation
int vdir_mode = VDIR_NONE;

//intervalle pour les pointilles
int ilc = 2;

//zooming variables
int zoom_offset_x = 0;
int zoom_offset_y = 0;
float zoom_coef = 1.0;

void set_light_xyz();

void bad_options(char *opt)
{
  ErrPrintf("erreur lecture option : %s\n", opt);

  exit(-4);
}

void show_view_options()
{
  printf("\nGRAPHICAL OPTIONS\n");
  printf("  -ch <val> \t horizontal cross section of height <val>\n");
  printf("  -cv0 \t\t centered vertical cross sections\n");
  printf("  -cv <a>x<o> \t vertical cross sections at abscissa <a> and ordinate <o>\n");
  //printf("  -diag \t diagonal layers of cells \n");
  printf("  -lc \t\t show cut lines\n");
#ifdef ALTI
  printf("  -ls <angle> \t surface rendering with ligth source elevation = <angle> degrees\n");
#endif
  printf("  -al \t\t height map\n");
  printf("  -lal \t\t low height map\n");
  printf("  -tr <cel> \t transparent cells of type <cel>\n");
//#if defined(MODEL_DUN) || defined(MODEL_RIV) || defined(MODEL_AVA)
#ifdef LGCA
  if (prog == PROG_RESCAL)
    printf("  -vel \t\t view flow mean velocity on a vertical plan\n");
    printf("  -vel2 \t view flow mean velocity and vector field on a vertical plan\n");
#if defined(CGV)
  if (prog == PROG_RESCAL)
    printf("  -vss \t\t view shear stress\n");
    printf("  -vss2 \t view renormalized shear stress between min and max thresholds\n");
#endif
#endif
  fflush(stdout);
}

int view_init(int argc, char *argv[])
{
  int i;
  unsigned char typ;

  hh = H>>1;
  memset(visible, 1, MAX_CELL);

  for(i=1; i<argc; i++){
    if (!strcmp(argv[i],"-cv0")){
      opt_cv = 1;
      abs_cv = L>>1;
      prof_cv = D>>1;
    }
    else if (!strcmp(argv[i],"-cv")){
      char c;
      int n=0;
      opt_cv = 1;
      if (i+1<argc) n = sscanf(argv[++i], "%d%c%d", &abs_cv, &c, &prof_cv);
      //LogPrintf("n=%d\n",n);
      if (n<3) bad_options("-cv <abs>x<prof>");
      abs_cv+=L_bounds;
      prof_cv+=D_bounds;
      //LogPrintf("abs_cv=%d   prof_cv=%d\n", abs_cv, prof_cv);
    }
    else if (!strcmp(argv[i],"-lc"))
      opt_lc = 1;
    else if (!strcmp(argv[i],"-ch")){
      opt_ch = 1;
      hh = atoi(argv[++i]);
      if (hh>H-2) hh=H-2;
    }
    else if (!strcmp(argv[i],"-tr")){
      opt_tr = 1;
      typ = atoi(argv[++i]);
      visible[typ] = 0;
    }
    else if (!strcmp(argv[i],"-al")){
      opt_al = 1;
      if ((i+1<argc) && (*argv[i+1]!='-')) i++; /// for backward compatibility
    }
    else if (!strcmp(argv[i],"-lal")){
      opt_lal = 1;
      if ((i+1<argc) && (*argv[i+1]!='-')) i++; /// for backward compatibility
    }
#ifdef ALTI
    else if (!strcmp(argv[i],"-ls")){
      opt_ls = 1;
      light_elevation = atof(argv[++i]);
      set_light_xyz();
    }
#endif
#ifdef LGCA
    else if ((prog == PROG_RESCAL) && (use_lgca) && (!strcmp(argv[i],"-vel")))
      opt_vel = 1;
    //LogPrintf("opt_vel = %d\n", opt_vel);
    else if ((prog == PROG_RESCAL) && (use_lgca) && (!strcmp(argv[i],"-vel2")))
      opt_vel = 2;
    //LogPrintf("opt_vel = %d\n", opt_vel);
#if defined(CGV)
    else if ((prog == PROG_RESCAL) && (use_lgca) && (!strcmp(argv[i],"-vss")))
      opt_vss = 1;
    else if ((prog == PROG_RESCAL) && (use_lgca) && (!strcmp(argv[i],"-vss2")))
      opt_vss = 2;
#endif
#endif
  }

  if (!visible[BORD]){
    ErrPrintf("ERROR: do not specify %d cells as transparent\n", BORD);
    exit(-1);
  }

  if (LNS == 1) opt_lc = opt_vss = 0;

  if (opt_vel && !opt_cv) prof_cv = D>>1;

  if (opt_cv || opt_vel){
    LogPrintf("abs_cv = %d\n", abs_cv);
    LogPrintf("prof_cv = %d\n", prof_cv);
  }

  img_w = (L<=IMG_W_MAX)? L : IMG_W_MAX;
  if (L>IMG_W_MAX) {LogPrintf("largeur de l'image fixee a IMG_W_MAX = %d\n", IMG_W_MAX); }
  img_h = D;//LNS+2;//L;
  if (opt_cv){
    if (LNS>1) img_w += H;
    img_h += H; //(LNS>1) ? 2*H : H;
  }
  if (opt_vel && (!opt_vss || !opt_cv || (L>D))) img_h += H;
#if defined(CGV)
  if (opt_vss){
    horizontal_display = (L<=D);
    if (horizontal_display)
      img_w += L;
    else
      img_h += LNS+2;
  }
#ifdef REORIENT_AUTO
  LogPrintf("rot_mode = %d\n",rot_mode);
  //if (opt_vss && (rot_mode == ROT_MODE_FULL)){ WarnPrintf("WARNING: display of the shear stress is not fully supported with FULL rotations\n");}
#endif
#endif

  if (!image){
    /*image = (unsigned char *) malloc(img_w*img_h);
  	if (!image){
          ErrPrintf("erreur malloc image (img_w = %d   img_h = %d)\n", img_w, img_h);
		  exit(-4);
	  }*/
    AllocMemoryPrint("image",image, unsigned char, img_w*img_h);
    PrintTotalMemory();
    memset(image, BORD, img_w*img_h);
    image_end = image + img_w*img_h;
  }

  return 0;
}

void view_img_size(int *pL, int *pH)
{
  *pL = img_w;
  *pH = img_h;
}

void init_shading(int index, int nb_values, int col0, int col1)
{
  assert(index < MAX_CELL);
  shades[index].col0 = col0; //dark brown
  shades[index].col1 = col1; //light brown
  shades[index].nb_val = nb_values;
  shades[index].start = (index==0) ? MAX_CELL+1 : shades[index-1].start + shades[index-1].nb_val;
  assert(shades[index].start + nb_values <256);
}

unsigned char view_shading(int n, float val)
{
  return shades[n].start + val*(shades[n].nb_val-1);
}


void view_palette(int *colors)
{
  static int start=1;
  static int palette[256];
  int i;

  if (start==0){
    memcpy(colors, palette, 256*sizeof(int));
    return;
  }
  start=0;

#if defined(MODEL_DUN)
//#define GRC GRV //colored grain

  palette[GR] = 0x00400000; //0x00000000; //grain
  palette[GRJ] = 0x00ff0000; //mobile grain
  //palette[GRV] = 0x00ff00ff; //not used
#ifdef BR
  palette[BR] = 0x007f3f00; //bedrock
#endif
#ifdef GRV
  palette[GRV] = 0x00009f00; //0x003f5f00; //vegetated grain
#endif
  palette[EAUC] = 0x0000ffff; //air (or water)
  palette[VEG] = 0x00009f00; //vegetation
  palette[BORD] = 0x00c0c0c0; //boundary
  palette[DUM] = 0x007f7f7f;  //0x00ffff00; //neutral
  palette[IN] = 0x00ffff00; //0x007f7f00; //input of sand
  palette[OUT] = 0x009f0000; //output of sand
  palette[EAUT] = 0x000000ff; //turbulence
  palette[GRC] = 0x00ffffff; //colored grain

#elif defined(MODEL_SNO)
  palette[EAUC] = 0x00e5d0ff; // light purple fluid
  palette[EAUT] = 0x00bf8bff; // mid-purple turbulence
  palette[GR]   = 0x00ffffff; // white grain
  palette[GRJ]  = 0x00ff0000; // red mobile grain
  palette[GRV]  = 0x00adcbe3; // light-blue hardened grain
  palette[BR]   = 0x00000000; // black bedrock
  palette[DUM]  = 0x00000000; // black DUM/neutral
  palette[IN]   = 0x00ffdc73; // yellow input
  palette[OUT]  = 0x00ffbf00; // dark yellow output
  palette[GRC]  = 0x00ffc0cb; // pink colored grain
  palette[BORD] = 0x00c0c0c0; // grey boundary

#elif defined(MODEL_AVA)

#define GRC 6 //colored grain

  palette[0] = 0x00000000; //grain
  palette[1] = 0x0000ffff; //air
  palette[2] = 0x007f7f7f; //bord
  palette[3] = 0x007f7f7f; //0x00ffff00; //inerte
  palette[4] = 0x00ffff00; //0x007f7f00; //source de grains
  palette[5] = 0x009f0000; //sortie de grains
  palette[GRC] = 0x00ffffff; //0x00009f00; //grain colore

#elif defined(MODEL_RIV)

  palette[TERRE] = 0x00400000; //0x00000000; //terre
  palette[BOUE] = 0x00406080; //boue
  palette[EAU] = 0x0000ffff; //eau
  palette[BT] = 0x0080a0c0; //0x00ff7f00; //boue turbulente
  palette[PIERRE] = 0x0000c8c8; //pierre
  palette[BORD] = 0x00c0c0c0; //0x007f7f7f; //bord
  //palette[6] = 0x000000ff; //inerte
  palette[DUM] = 0x00ffff00; //inerte
  palette[IN] = 0x007f7f00; //source de grains
  palette[OUT] = 0x009f0000; //sortie de grains

#elif defined(MODEL_DIF)

  palette[ZERO] = 0x00ffffff;
  palette[ONE] = 0x00000000;
  palette[IN] = 0x007f7f00;
  palette[DUM] = 0x00ffff00;

#elif defined(MODEL_D2G)

  palette[TERRE] = 0x00000000;
  palette[AIR] = 0x000000ff;
  palette[EAU] = 0x00ffffff;
  palette[BOUE] = 0x00cf9f5f;
  palette[DUM] = 0x00ffff00;
  palette[IN] = 0x007f7f00;

#elif defined(MODEL_LIFE)

  palette[DEAD] = 0x00ffffff;
  palette[ALIVE] = 0x00000000;
  palette[BORD] = 0x007f7f7f;
  palette[DUM] = 0x000000ff;

#else

  palette[0] = 0x00000000; //PLUS
  palette[1] = 0x00ffffff; //ZERO
  palette[2] = 0x00ff0000; //MOINS
  palette[BORD] = 0x007f7f7f;
  palette[DUM] = 0x000000ff;
  palette[5] = 0x007f7f00; //IN
  palette[6] = 0x00afafaf; //BT
  palette[7] = 0x0000c8c8; //PIERRE

#endif

#ifdef PARALLEL
  palette[TUNNEL] = 0x007f007f;
#endif

#ifdef MODEL_CMB
  palette[COUPE] = 0x00007f00;
#else
  palette[COUPE] = 0x00000000;
#endif

#ifdef ALTI
  /// color palette for shadings
  unsigned char red, green, blue;
  int r0, g0, b0, r1, g1, b1;
  float coef;

  //ResetMemory(cel_shade_index, int, MAX_CELL);
  for(i=0; i<MAX_CELL; i++){
    cel_shade_index[i] = -1;
  }
#ifdef MODEL_DUN
  visible[EAUC] = 0; /// EAUC cells not visible
  if (opt_ls | opt_al | opt_lal) visible[GRJ] = 0; /// GRJ cells not visible
#endif // MODEL_DUN

#ifdef MODEL_SNO
  visible[EAUC] = 0; // EAUC cells not visible
  if (opt_ls | opt_al | opt_lal) visible[GRJ] = 0; // jumping GRJ cells not visible
#endif // MODEL_SNO

#ifdef MODEL_AVA
  visible[AIR] = 0; /// AIR cells not visible
#endif // MODEL_DUN

#ifdef MODEL_RIV
  visible[EAU] = 0; /// EAU cells not visible
#endif // MODEL_RIV

  /// number of shadings
  nb_shades = 1;
#ifdef MODEL_DUN
  nb_shades++;
#endif // MODEL_DUN
#ifdef MODEL_SNO
  nb_shades++;
#endif // MODEL_SNO
#ifdef CELL_COLOR
  nb_shades++;
#endif
#ifdef USE_VEGETATION
  if (use_veg) nb_shades++;
#endif
#ifdef MODEL_RIV
  nb_shades+=2;
#endif

  int nb_values = (255 - MAX_CELL - 1)/nb_shades;
  LogPrintf("nb_shades=%d\n", nb_shades);
  LogPrintf("nb_values=%d\n", nb_values);

  /// default shading
  int index=0;
  init_shading(index, nb_values, 0x00400000, 0x00ffc37f); //brown shading
  cel_shade_index[ALTI] = index;

#ifdef MODEL_DUN
   /// DUM shading
  index++;
  init_shading(index, nb_values, 0x00101010, 0x00909090); //dark grey shading
  cel_shade_index[DUM] = index;
#endif // defined

#ifdef MODEL_SNO
  // DUM shading
  index ++;
  init_shading(index, nb_values, 0x0036454f, 0x00708090); //blue-grey shading
  cel_shade_index[DUM] = index;
#endif

#ifdef CELL_COLOR
  /// "colored" cells
  index++;
  init_shading(index, nb_values, 0x00202020, 0x00e0e0e0); //light grey shading
  col_shade_index = index;
#endif // CELL_COLOR

#ifdef USE_VEGETATION
  if (use_veg){
    /// shading on vegetation
    index++;
    init_shading(index, nb_values, 0x00002000, 0x003f9f1f); //green shading
    cel_shade_index[GRV] = index;
    cel_shade_index[VEG] = index;
  }
#endif // USE_VEGETATION

#ifdef MODEL_RIV
  /// shading on mud cells
  index++;
  //init_shading(index, nb_values, 0x00e0a0a0, 0x00ffffff); //light brown shading
  init_shading(index, nb_values, 0x00406080, 0x0080a0c0); //bluish grey shading
  cel_shade_index[BOUE] = index;
  index++;
  //init_shading(index, nb_values, 0x00ff7f00, 0x00ffffaf); //orange shading
  //init_shading(index, nb_values, 0x004080c0, 0x0070b0f0); //light brown shading
  init_shading(index, nb_values, 0x00406080, 0x0080a0c0); //bluish grey shading
  cel_shade_index[BT] = index;
#endif // MODEL_RIV

  /// color palette for all shadings
  for(index=0; index<nb_shades; index++){
    r0 = (shades[index].col0 & 0x00ff0000) >> 16;
    g0 = (shades[index].col0 & 0x0000ff00) >> 8;
    b0 = (shades[index].col0 & 0x000000ff);
    r1 = (shades[index].col1 & 0x00ff0000) >> 16;
    g1 = (shades[index].col1 & 0x0000ff00) >> 8;
    b1 = (shades[index].col1 & 0x000000ff);
    for(i=0; i<shades[index].nb_val; i++) {
      coef = (float)i/(shades[index].nb_val-1);
      red = r1*coef + r0*(1-coef);
      green = g1*coef + g0*(1-coef);
      blue = b1*coef + b0*(1-coef);
      palette[shades[index].start+i] = (red << 16) | (green << 8) | blue;
    }
  }
#endif // ALTI

  memcpy(colors, palette, 256*sizeof(int));
}


// draw line
void view_line(unsigned char *pix0, int dx, int dy, unsigned char color)
{
  int sdx,sdy,px,py,dxabs,dyabs,i;
  float slope;
  unsigned char *pix;

  /*if (dx*dx+dy*dy > 0){
    LogPrintf("view_line() dx = %d, dy = %d\n", dx, dy);
  }*/

  //horizontal clipping
  int ipix = (int)(pix0 - image);
  int xpix = ipix % img_w;
  Check_min_max(dx, -xpix, img_w-1-xpix);

  dxabs = abs(dx);
  dyabs = abs(dy);
  sdx = sgn(dx);
  sdy = sgn(dy);
  i = 0;
  if (dxabs >= dyabs) /* the line is more horizontal than vertical */
  {
    slope = (float)dy / (float)dx;
    while(1)
    {
      py = roundf(slope*i);
      pix = pix0+i+py*img_w;
      if ((pix >= image) && (pix < image_end)) *pix = color;
      if (i==dx) break;
      i+=sdx;
    }
  }
  else /* the line is more vertical than horizontal */
  {
    slope = (float)dx / (float)dy;
    while(1)
    {
      px = roundf(slope*i);
      pix = pix0+px+i*img_w;
      if ((pix >= image) && (pix < image_end)) *pix = color;
      if (i==dy) break;
      i+=sdy;
    }
  }
}

// draw filled rectangle
void view_rect(unsigned char *pix0, int width, int height, unsigned char color)
{
  unsigned char *pix;
  int y;

  assert(pix0 + (height-1)*img_w + width <= image_end);

  if (width <= 0) return;

  pix = pix0;

  for (y=0; y < height; y++, pix+=img_w)
    memset(pix, color, width);
}

int view_alti(int i, int j, int low_flag)
{
  Cell *pt;
  int k;

  pt = TE + i + j*HL;
  if (!low_flag){
    /// find first visible cell starting from the top
    k=1+h_ceil;
    if (opt_ch) k=khh;
    while ((k < H) && (visible[(pt+k*L)->celltype]==0)) k++;
  }
  else{
    /// find highest visible cell starting from the bottom
    int kmin=0;
    if (opt_ch) kmin=khh;
    k=H-1-h_floor;
    while ((k > kmin) && visible[(pt+(k-1)*L)->celltype]) k--;
  }

  return k;
}

unsigned char view_color(Cell *pt)
{
#if defined(CELL_COLOR) && defined(GRC)
  return (pt->color) ? GRC : pt->celltype;
#else
  return pt->celltype;
#endif // CELL_COLOR
}

void set_light_xyz()
{
  float alpha = light_azimuth*PI/180;
  float phi = light_elevation*PI/180;
  light.x = cosf(phi)*sinf(alpha); //1.0;
  light.y = sinf(phi); //0.0;
  light.z = -cosf(phi)*cosf(alpha);
}

void rotate_light(float angle)
{
  light_azimuth += angle;
  //LogPrintf("rotate_light: angle = %f   azimuth = %f\n", angle, light_azimuth);
  set_light_xyz();
}

float view_light(VECNORM *ptn)
{
  float ldiff;

#ifdef NORM2D
  ldiff = (ptn->x)*light.x + (ptn->y)*light.y; //lumiere diffuse (2d)
#else
  ldiff = (ptn->x)*light.x + (ptn->y)*light.y + (ptn->z)*light.z; //lumiere diffuse
#endif
  if (ldiff < 0.0) ldiff = 0.0;

  return ldiff;
}

unsigned char * view()
{
  Cell *pt;
  unsigned char *pix, *cur_image;
  int i,j,k,typ;
  int shade_index;

  pix = cur_image = image;

  view_rect(image, img_w, img_h, BORD);

  /// update the position of the horizontal cross section
  khh = H-1-hh;

#ifdef ALTI
  if (opt_ls){
    /// relief light shading
    if (norm3d){
      VECNORM *ptn;
      assert (csol);
      for(j=1; j < D-1; j++){
#ifdef NORM2D
        ptn = norm2d + (j-LN)*LEO;
#else
        ptn = norm3d + (j-LN)*LEO;
#endif
        //if ((j-LN<0) && !psol[j]) {ErrPrintf("ERROR: j-LN=%d (%d - %d)\n",j-LN,j,LN); exit(-1);}
        pix = cur_image + j*img_w + 1;
        for(i=1; i < L-1; i++, ptn++){
          if (csol[j*L+i])
              *pix++ = csol[j*L+i];
          else{
            k = view_alti(i,j,0);
            pt = TE+j*HL+i+k*L;
            typ = pt->celltype;
            shade_index = cel_shade_index[typ];
#ifdef CELL_COLOR
            if (pt->color) shade_index = col_shade_index;
#endif // CELL_COLOR
            //if (shade_index >= 0)
            if ((shade_index < 0) || (opt_ch && (k==khh))){
              *pix++ = view_color(pt);
            }
            else{
              *pix++ = view_shading(shade_index, view_light(ptn));
            }
          }
        }
      }
    }
  }
  else if (opt_al || opt_lal){
    /// height map
    for(j=1; j < D-1; j++){
      pix = cur_image + j*img_w + 1;
      for(i=1; i < L-1; i++){
//#ifdef TRACE_SRC
#if defined(TRACE_SRC) || defined(TRACE_AIRE)
        if (trace_point(i,j)){
          *pix = COUPE;
        }
        else
#endif
        {
          k = view_alti(i,j,opt_lal);
          pt = TE+j*HL+i+k*L;
          typ = pt->celltype;
          shade_index = cel_shade_index[typ];
#ifdef CELL_COLOR
          if (pt->color) shade_index = col_shade_index;
#endif // CELL_COLOR
          if ((shade_index < 0) || (opt_ch && (k==khh))){
            *pix++ = view_color(pt);
          }
          else{
            *pix++ = view_shading(shade_index, 1.0*(H-k)/H);
          }
        }
      }
    }
  }
  else
#endif //ALTI
  {
    /// by default, display horizontal view of the (predefined) visible cells
    /// or cross section (with opt_ch option)
    /// without shading
    for(j=1; j < D-1; j++){
      pix = cur_image + j*img_w + 1;
      for(i=1; i < L-1; i++){
        k = opt_ch ? khh : view_alti(i,j,0);
        pt = TE+j*HL+i+k*L;
        *pix++ = view_color(pt);
      }
    }
  }

  memset(pix,BORD,img_w);
  //pix += img_w;
  cur_image += D*img_w;

  if (opt_cv){
    pix = cur_image;
    /// vertical cross section, east-west
    for(j=1; j < H-1; j++){
      pt = TE + HL*(prof_cv) + j*L + 1;
      pix = cur_image + j*img_w + 1;
      for(i=1; i < L-1; i++, pt++){
        *pix++ = view_color(pt);
      }
    }
    if (LNS>1){
      /// vertical cross section, north-south
      for(i=1; i < D-1; i++){
        pt = TE + i*HL + abs_cv + L;
        pix = image + i*img_w + L + 1;
        for(j=1; j < H-1; j++, pt+=L) {
          *pix++ = view_color(pt);
        }
      }
      // empty space (antibug)
      view_rect(cur_image + L, H, H, BORD);
    }


#ifdef ROTATIONS
    if ((vdir_mode != VDIR_NONE) && (LNS>1)){
      /// representation of current orientation
      extern float csp_angle;
      int Hd2 = H/2;
      int l1 = Hd2;
      int l2 = l1/2;
      float alpha = 0.0;
      if (vdir_mode == VDIR_NORTH){
        //TODO
      } else if (vdir_mode == VDIR_WIND){
        alpha = -csp_angle*PI/180.0;
      }
      float co = cosf(alpha);
      float si = sinf(alpha);
      pix = cur_image + L + Hd2 + Hd2*img_w;
      pix -= (int)(0.5*l1*co);
      pix -= (int)(0.5*l1*si)*img_w;
      view_line(pix, l1*co, l1*si, COUPE);
      pix += (int)(0.5*l2*si);
      pix -= (int)(0.5*l2*co)*img_w;
      view_line(pix, -l2*si, l2*co, COUPE);
    }
#endif
    cur_image += H*img_w;
  }
  else if (opt_vel && opt_lc && !reorient_flag){
    /// cut line for the vertical velocity plane
    pix = image + img_w*(prof_cv) + 1;
    for (i=1; i<img_w-1; i+=ilc, pix+=ilc)
      if (i!=L-1) *pix = COUPE;
  }

#ifdef LGCA
//#if 0
  unsigned char *img_col = cur_image;//image + img_w*(img_h-H);
  if (use_lgca){
    if (opt_vel) cur_image += H*img_w;
    if (horizontal_display)
      img_col = image + L + D*img_w + (opt_cv? H : 0);
    //else if (opt_vss)
    //  view_rect(img_col + L, img_w - L, H, BORD);

#if 0
    /// basic rendering (vertical vs horizontal) of lgca particles
    if (opt_vel){
    //if (0){
      unsigned char *pix0;
      int index, vcol;
      pix0 = img_col;
      index = HL*prof_cv;//(CLNS-1);//(CLNS>>1);
      for(i=0; (i<HL) ; i++, index++, pix0++){
        //*pix0 = check_mvt(index, NULL) ? COUPE : check_mvt_solid(index, NULL) ? DUM : BORD;
        //*pix0 = check_mvt_ouest(index, NULL) ? GRJ : check_mvt_est(index, NULL) ? EAUT : check_mvt_solid(index, NULL) ? DUM : BORD;
        //*pix0 = check_mvt_haut(index, NULL) ? GRJ : check_mvt_bas(index, NULL) ? EAUT : check_mvt_solid(index, NULL) ? DUM : BORD;
        *pix0 = check_mvt_haut(index, NULL) || check_mvt_bas(index, NULL) ? GRJ : check_mvt_ouest(index, NULL) || check_mvt_est(index, NULL) ? EAUT : check_mvt_solid(index, NULL) ? DUM : BORD;
      }
    }
#elif 1
    /// mean velocity (interpolation) of the flow in a vertical plane
    if (opt_vel && Velx_interp && Vely_interp){
      if ((prof_cv>=LN) && (prof_cv<=LS)){
        unsigned char *pix0;
        int index, sol_flag;
        float vcol;
        float vcolmax = sqrtf(2)*VSTEP_H*VSTEP_L*VSTEP_TIME;
        pix0 = img_col;
        index = HL*prof_cv;
        for(j=0; j<H; j++){
          pix0 = img_col + j*img_w;
          for(i=0; (i<L) ; i++, index++, pix0++){
            vcol = sqrtf(Velx_interp[index]*Velx_interp[index]+Vely_interp[index]*Vely_interp[index]);
            if (vcol > vcolmax) vcol = vcolmax;
            sol_flag = reorient_flag ? check_mvt_solid(index, NULL) : (Phase[TE[index].celltype]==SOLID);
            *pix0 = sol_flag ? BORD : (Velx_interp[index]<0) ? EAUT : view_shading(0, vcol/vcolmax);
          }
        }
      }
      else{
        for(j=0; j<H; j++) memset(img_col + j*img_w, DUM, L);
      }
    }

    /// rendering of the velocity field with arrows (vertical plane)
    if ((opt_vel >= 2) && Velx_interp && Vely_interp){
      int iv, jv;
      float *pvx0, *pvy0, *pvx, *pvy;
      unsigned char *pix0;
      //int fv=VSTEP*VSTEP_TIME*NB_MVT_EO*NB_MVT_VER/4;
      int fv=VSTEP_H*VSTEP_L*VSTEP_TIME/20;
      int planz = (int)(prof_cv-LN+1)/DIST_MVT_NS; //(CLNS>>1); //LNS-1;//(CLNS>>1);
      Check_min_max(planz, 0, CLNS-1);
      //int fy=2*VSTEP_TIME*NB_MVT_EO*NB_MVT_VER;
      //memset(img_col, BORD, HL);
      pvx0 = Velx_interp+planz*HL;
      pvy0 = Vely_interp+planz*HL;
      pix0 = img_col; // + (VSTEP>>1)*img_w + (VSTEP>>1);
      //view_rect(img_col, L, H, BORD);
      for(jv=1; jv<VH-1; jv++){ //hauteur
        pix = pix0 + jv*VSTEP_H*img_w + VSTEP_L;
        pvx = pvx0 + jv*VSTEP_H*L + VSTEP_L;
        pvy = pvy0 + jv*VSTEP_H*L + VSTEP_L;
        for(iv=1; iv<VL-3; iv+=2, pvx+=2*VSTEP_L, pvy+=2*VSTEP_L, pix += 2*(VSTEP_L/NB_MVT_EO)){ //largeur
          if ((jv<VH-1) && (iv<VL-1)){
            //view_line(pix, (*pvx)>>1, (*pvy)>>1, 0);
            //view_line(pix, ((*pvx)/VSTEP_TIME)>>1, ((*pvy)/VSTEP_TIME)>>1, 0);
            view_line(pix, (*pvx)/fv, (*pvy)/fv, GRC);
      /*if ((iv>0) && (jv>0)){
        //calcul vorticite
            int vort = (*pvy) - (*(pvy-1)) - ((*pvx) - (*(pvx-VL)));
            if (vort > 10) *pix = *(pix+1) = *(pix+L) = *(pix+L+1) = MOINS;
            else if (vort < -10) *pix = *(pix+1) = *(pix+L) = *(pix+L+1) = DUM;
            }*/
          }
        }
      }
    }
#else
      /// basic rendering (east vs west) of lgca particles
      unsigned char *pix0;
      int ii0, ii, x, y;
      pix0 = img_col;
      ii0 = (HL*NB_MVT_EO)*(int)((prof_cv-LN)/DIST_MVT_NS);//+(L*NB_MVT_EO*0.15);
      //ii0 = (HL*NB_MVT_EO)*(L>>1)+(L*(NB_MVT_EO-1));
      for(y=0; y<H ; y++, ii0+=L*NB_MVT_EO){
        for(x=0, ii=ii0; x<L; x++, pix0++, ii+=NB_MVT_EO){
          *pix0 = mvt_ouest(ii) ? EAUT : mvt_est(ii) ? GRJ : mvt(ii) ? EAUC : mvt_solid(ii) ? DUM : GR;
          //*pix0 = mvt_haut(ii) ? GRJ : mvt_est(ii) ? EAUT : mvt_solid(ii) ? GR : EAUC;
          //*pix0 = mvt_haut(ii) ? GRJ : mvt_bas(ii) ? EAUT : mvt_solid(ii) ? DUM : BORD;
          //*pix0 = mvt_haut(ii) || mvt_bas(ii) ? GRJ : mvt_ouest(ii) || mvt_est(ii) ? EAUT : mvt_solid(ii) ? DUM : BORD;
        }
      }
#endif

#ifdef CGV
    /// display normals in the vertical plane of the flow
    if (opt_vel /*&& opt_vss*/ && alti && norm3d && Velx_interp && (prof_cv>=LN) && (prof_cv<=LS) && !reorient_flag){
    //if (0){
      float ln=10;
#ifdef CGV
      extern float dist_grdv;
      ln = dist_grdv;
#endif
      VECNORM *pn;
      short *alt, x, y, z, stepx;
      z = prof_cv-LN;
      stepx = 10;//10; //VSTEP*2; //*2;
      alt = alti+z*LEO+(stepx>>1);
#ifdef NORM2D
      pn = &norm2d[z*LEO+(stepx>>1)];
#else
      pn = &norm3d[z*LEO+(stepx>>1)];
#endif
      for(x=(stepx>>1); x<LEO; x+=stepx, alt+=stepx, pn+=stepx){
      //for(x=1; x<L-2; x+=2, alt+=2, pnx+=2, pny+=2){
        if (rot_map && OutOfSpace(x+1, z+1)) continue;
        y = (H-1)-(*alt);
        pix = img_col + y*img_w + x + 1;
        view_line(pix, (int)(ln*(pn->x))/*/NB_MVT_EO*/, (int)(-ln*(pn->y)), GRJ);
        //if (*pnx > 0.5) LogPrintf("nx = %f   ny = %f\n", *pnx, *pny);
      }
      //view_rect(img_col+L,img_w-L,H,BORD);
    }
#endif

#ifdef CGV
    /// display shear stress in dedicated horizontal plane
    extern float *grdv;
    extern float grdvc, grdvc_min, grdvc_max;

    if (horizontal_display){
      // horizontal disposition
      cur_image = image + L;
      if (opt_cv) cur_image += H;
    }

    if (opt_vss /*&& !reorient_flag*/){
      if (grdv && lgca_ready){
        //view_rect(cur_image, L, LN, BORD);
        float *ptgv, prob, vcol;
        float vcolmax = 4*grdvc;
        for(j=0; j < LNS; j++){
          pix = cur_image+(LN+j)*img_w;
          if (!psol[LN+j] || reorient_flag){ //test si plan solide
            ptgv = grdv + j*LEO;
            *pix++ = BORD;
            for(i=0; i < LEO; i++, ptgv++){
              //if (rot_map && OutOfSpace(1+i, LN+j)){*pix++ = BORD; continue;}
              if (opt_vss==1){
                vcol = (*ptgv);
                //if (vcol > vcolmax) vcol = vcolmax;
                *pix++ = (vcol < 0) ? EAUT : (vcol <= vcolmax) ? view_shading(0, vcol/vcolmax) : GRJ;
              }
              else{
                prob = prob_cgv(*ptgv, grdvc_min, grdvc_max);
                *pix++ = view_shading(0, prob);
              }
            }
            *pix++ = BORD;
          }
          else memset(pix, BORD, L);
        }
        //view_rect(cur_image + LS*img_w, L, D-1-LS, BORD);
        /*if (opt_lc && (prof_cv>=LN) && (prof_cv<=LS)){
          //ligne de coupe
          //extern int CHL;
          //int cx, cy, cz;
          //Calcule_cxyz(HL*prof_cv+1, cx, cy, cz);
          //int cz = (int)(prof_cv-LN)/DIST_MVT_NS;
          //memset(cur_image+img_w+img_w*cz*DIST_MVT_NS, COUPE, L);
          //memset(cur_image+img_w*(1+prof_cv-LN), COUPE, L);
          pix = cur_image+img_w*(prof_cv-LN);
          for(i=0; i<L; i+=ilc, pix +=ilc) *pix = COUPE;
        }*/
      }
      else{
        //view_rect(cur_image, LEO, LNS, BORD);
      }
    }
#ifndef NORM2D
    extern float dist_grdv;
    /// display normals on the horizontal plane of the shear stress
    if ((opt_vss==1) && alti && norm3d && lgca_ready && !reorient_flag){
    //if (0){
      float ln=10;
      ln = dist_grdv;
      Vec3 *pn;
      short x, z, stepx, stepz;
      stepx = 10;
      stepz = 10;
      for(z=(stepz>>1); z<LNS; z+=stepz){
        pn = &norm3d[z*LEO+(stepx>>1)];
        for(x=(stepx>>1); x<LEO; x+=stepx, pn+=stepx){
          if (rot_map && OutOfSpace(x+1, z+1)) continue;
          pix = cur_image + (z+1)*img_w + x + 1;
          view_line(pix, (int)(ln*(pn->x))/*/NB_MVT_EO*/, (int)(ln*(pn->z)), GRJ);
        }
      }
    }
#endif //NORM2D
#endif //CGV
  }
#endif //LGCA

  /// display cut lines (opt_lc option)
  if (opt_cv && opt_lc){
    /// east-west and north-south cut lines
    //memset(image+L*(1+prof_cv-LN), COUPE, L);
    pix = image + img_w*(prof_cv) + 1;
    for (i=1; i<img_w-1; i+=ilc, pix+=ilc)
      if (i!=L-1) *pix = COUPE;

    pix = image + img_w + abs_cv;
    for (j=1; j<D-2+H; j+=ilc, pix+=ilc*img_w)
      if (j!=D-1) *pix = COUPE;

    /// horizontal cut line in east-west vertical cross section (with opt_ch option)
    //if (!opt_tr && !opt_al && !opt_lal && !opt_ls){
    if (opt_ch){
      //memset(cur_image+hh*L, COUPE, L);
      pix = image + (D+khh)*img_w + 1;
      for (i=1; i<L-1; i+=ilc, pix+=ilc){
        *pix = COUPE;
      }
    }
  }

  return image;
}

/// Update the position of horizontal/vertical cross sections
int update_cv(int x, int y, int flag)
{
  int ret=0;
  int dist=5; //max distance to catch the line
  static int flagx=0;
  static int flagy=0;
  static int flagz=0;

  if (!opt_cv || !opt_lc) return ret;

  /// compute position in case of a zoomed area
  x = (x - zoom_offset_x)/zoom_coef;
  y = (y - zoom_offset_y)/zoom_coef;

  if (!flag){
    /// caption of a cut line
    if ((x<L+H) && (y<D+H)){
      flagx = (abs(abs_cv-x) <= dist);
      flagy = (abs(prof_cv-y) <= dist);
      flagz = (opt_ch && (abs(D+khh-y) <= dist));
    }
    else{
      flagx = flagy = flagz = 0;
    }
  }

  /// controlled motion of a cut line
  if (flagx) abs_cv = x;
  if (flagy) prof_cv = y;
  if (flagz) hh = D+H-y;

  /// check boundaries
  if (abs_cv < 1) abs_cv = 1;
  if (abs_cv > L-2) abs_cv = L-2;
  if (prof_cv < 1) prof_cv = 1;
  if (prof_cv > D-2) prof_cv = D-2;
  if (hh < 1) hh = 1;
  if (hh > H-2) hh = H-2;

  ret = 1;

  return ret;
}

void reset_zoom()
{
  zoom_coef = 1.0;
  zoom_offset_x = 0;
  zoom_offset_y = 0;
}

unsigned char* view_zoom(int width, int height, float coef)
{
  static unsigned char *image_zoom = NULL;
  static int zoom_size = 0;
  int zoom_w = Min(img_w*coef, width);
  int zoom_h = Min(img_h*coef, height);
  int size = width*height;
  unsigned char *pix;
  int i, j, imin, imax, jmin, jmax, i0, j0;

  if (zoom_size < size){
    LogPrintf("width=%d, height=%d\n", width, height);
    if (!zoom_size) zoom_size = img_w*img_h;
    while (zoom_size < size) zoom_size *= 1.5;
    if (!image_zoom){
      AllocMemoryPrint("image_zoom", image_zoom, unsigned char, zoom_size);
    }
    else{
      ReallocMemoryPrint("image_zoom", image_zoom, unsigned char, zoom_size);
    }
  }
  memset(image_zoom, BORD, size);
  pix = image_zoom;
  jmin = (height-zoom_h) >> 1;
  jmax = jmin + zoom_h;
  imin = (width-zoom_w) >> 1;
  imax = imin + zoom_w;
  zoom_offset_x = imin;
  zoom_offset_y = jmin;
  zoom_coef = coef;

  for(j=jmin; j < jmax; j++){
    j0 = roundf((j-jmin)/coef);
    if (j0>=img_h) j0 = img_h-1;
    pix = image_zoom + j*width + imin;
    for(i=imin; i < imax; i++, pix++){
      i0 = roundf((i-imin)/coef);
      *pix = *(image + j0*img_w + i0);
    }
  }

  return image_zoom;
}

#ifdef USE_LIBPNG

#include <png.h>

png_color* png_palette;

void view_dump_init()
{
  int colors[256];
  int i;

  LogPrintf("view_dump_init\n");

  AllocMemory(png_palette, png_color, 256);
  ResetMemory(png_palette, png_color, 256);

  view_palette(colors);

  for(i=0; i<256; i++){
    png_palette[i].red   = (colors[i] & 0x00ff0000) >> 16;
    png_palette[i].green = (colors[i] & 0x0000ff00) >> 8;
    png_palette[i].blue  = (colors[i] & 0x000000ff);
  }
}

void view_dump_png(char *filename)
{
  unsigned char *pix;
  int j;

  LogPrintf("view_dump_png: writing file %s\n", filename);

  FILE *fp = fopen(filename, "wb");
  if (!fp){
    ErrPrintf("ERROR: cannot open file %s\n", filename);
    exit(-4);
  }

  png_structp png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
  if (!png_ptr){
    ErrPrintf("ERROR: cannot create structure png_struct\n");
    exit(-4);
  }

  png_infop info_ptr = png_create_info_struct(png_ptr);
  if (!info_ptr)
  {
    png_destroy_write_struct(&png_ptr, (png_infopp)NULL);
    ErrPrintf("ERROR: cannot create structure png_info\n");
    exit(-4);
  }

  if (setjmp(png_jmpbuf(png_ptr)))
  {
    png_destroy_write_struct(&png_ptr, &info_ptr);
    fclose(fp);
    ErrPrintf("ERROR: cannot set jump in view_dump_png()\n");
    exit(-4);
  }

  png_init_io(png_ptr, fp);

  png_set_IHDR(png_ptr, info_ptr, img_w, img_h,
   8, PNG_COLOR_TYPE_PALETTE, PNG_INTERLACE_NONE,
   PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);

  png_set_PLTE(png_ptr, info_ptr, png_palette, 256);

  png_write_info(png_ptr, info_ptr);
  //png_write_png(png_ptr, info_ptr, 0, NULL);

  pix = view();

  for(j=0; j<img_h; j++, pix+=img_w)
    png_write_row(png_ptr, (png_byte*)pix);

  //png_write_flush(png_ptr);

  png_write_end(png_ptr, NULL);

  png_destroy_write_struct(&png_ptr, &info_ptr);

  fclose(fp);
}

void dump_image(char *filename, char *format)
{
  if (!strcmp(format,"png"))
    view_dump_png(filename);
  else{
    ErrPrintf("ERROR: jpeg format not supported\n");
    exit(-1);
  }
}

#elif defined(USE_GD)

#include <gd.h>

gdImagePtr gd_img=NULL;
int gd_col[256];

void view_dump_init()
{
  int colors[256];
  int i, red, green, blue;

  gd_img = gdImageCreate(img_w, img_h);

  view_palette(colors);

  for(i=0; i<256; i++){
    red   = (colors[i] & 0x00ff0000) >> 16;
    green = (colors[i] & 0x0000ff00) >> 8;
    blue  = (colors[i] & 0x000000ff);
    gd_col[i] = gdImageColorAllocate(gd_img, red, green, blue);
  }
}

void view_gd()
{
  unsigned char *pix;
  int i, j;

  assert(gd_img);

  pix = view();

  for(j=0; j<img_h; j++)
    for(i=0; i<img_w; i++)
      gdImageSetPixel(gd_img, i, j, gd_col[*pix++]);
}

void view_dump_png(char *nom)
{
  FILE *fp;

  view_gd();
  fp = fopen(nom,"w");
  if ( ! fp ){
    ErrPrintf("erreur ouverture fichier png\n");
    exit(-4);
  }
  gdImagePng(gd_img, fp);
  fclose(fp);
}

void view_dump_jpeg(char *nom)
{
  FILE *fp;

  view_gd();
  fp = fopen(nom,"w");
  if ( ! fp ){
    ErrPrintf("erreur ouverture fichier jpg\n");
    exit(-4);
  }
  gdImageJpeg(gd_img, fp, 50); /* 50 -> medium quality and compression */
  fclose(fp);
}

void dump_image(char *filename, char *format)
{
  if (!strcmp(format,"png"))
    view_dump_png(filename);
  else
    view_dump_jpeg(filename);
}

#else

#include <gtk/gtk.h>
extern GtkWidget *drawingarea;

void dump_image(char *filename, char *format)
{
  //LogPrintf("dump_image: %s - %s\n",filename,format);
  GdkPixbuf *pixbuf = gdk_pixbuf_get_from_drawable(NULL,drawingarea->window,NULL,0,0,0,0,img_w,img_h);
  gdk_pixbuf_save(pixbuf, filename, format, NULL, NULL);
}

#endif //USE_GD

void dump_image_inter(int inter, char *format)
{
  static int cpt_inter=0;
  static int cpt_snap=0;
  static char nom[100];
  static char str[100];
  static char flag_sec=0;
  int nmin;
  //char *ext;

#ifdef PARALLEL
  extern int proc_id;
  sprintf(str,"_%d",proc_id);
#else
  *str = 0;
#endif

  //ext = (format==IMG_PNG)? "png" : "jpg";
  if (inter){
    if (!cpt_inter) flag_sec = inter % 60; //fmodf(inter,60);
    nmin = cpt_inter/60;
    if (!flag_sec)
      sprintf(nom,"%s%04d%s.%s", MOD_NAME, nmin, str, format);
    else
      sprintf(nom,"%s%04d-%02d%s.%s", MOD_NAME, nmin, cpt_inter - 60*nmin, str, format);
    cpt_inter += inter;
  }
  else{
    sprintf(nom,"SNAP%03d.%s", cpt_snap++, format);
  }

  dump_image(nom, format);
}

void view_quit()
{
  if (image) free(image);
#ifdef USE_GD
  if (gd_img) free(gd_img);
#endif
}

