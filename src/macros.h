/* ReSCAL - Global macro definitions
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

#include <inttypes.h>
#include <stdint.h>

#ifdef _MAIN_
int64_t total_memory = 0; //taille totale de memoire allouee
#ifdef LOG_FILE
FILE *log_file = NULL;
#endif
#else
extern int64_t total_memory;
#ifdef LOG_FILE
extern FILE *log_file;
#endif
#endif //_MAIN_


#ifdef LOG_FILE
#define LogPrintf(...) \
    if (log_file) \
      {fprintf(log_file, __VA_ARGS__); fflush(log_file);} \
    else \
      {printf(__VA_ARGS__); fflush(stdout);}
#else
#define LogPrintf(...) {printf(__VA_ARGS__); fflush(stdout);}
//#define Log(...) fprintf(stderr, __VA_ARGS__)
#endif //LOG_FILE

#define WarnPrintf(...) { \
    LogPrintf(__VA_ARGS__); \
    fprintf(stderr, __VA_ARGS__); \
    fflush(stderr); \
  }

#define ErrPrintf(...) { \
    LogPrintf(__VA_ARGS__); \
    LogPrintf("%s (%d)\n", __FILE__, __LINE__); \
    fprintf(stderr, __VA_ARGS__); \
    fprintf(stderr, "%s (%d)\n", __FILE__, __LINE__); \
    fflush(stderr); \
  }


#define AllocMemory(adr,type,taille) \
{     \
  if( (adr=(type *)malloc(sizeof(type)*taille)) == NULL ){  \
    ErrPrintf("malloc error : %d\n",(int)sizeof(type)*taille);    \
    exit(-1);         \
  }           \
  total_memory += sizeof(type)*taille;          \
}

#define ReallocMemory(adr,type,taille) \
{     \
  if( (adr=(type *)realloc(adr,sizeof(type)*taille)) == NULL ){ \
    ErrPrintf("realloc error: %d\n",(int)sizeof(type)*taille);  \
    exit(-1);         \
  }           \
}

#define ResetMemory(adr,type,taille) memset(adr,0,sizeof(type)*taille);

#define PrintMemory(var,type,taille,flag) \
  int64_t taille_mem = sizeof(type)*taille; \
  char *str = flag ? "re" : ""; \
  if (taille_mem >= 1000000) \
    LogPrintf("%sallocation %s : %.2f Mo\n", str, var, (float)taille_mem/1000000) \
  else \
    LogPrintf("%sallocation %s : %" PRId64 "\n", str, var, taille_mem)

#define AllocMemoryPrint(var,adr,type,taille) \
{ \
  PrintMemory(var,type,taille,0); \
  AllocMemory(adr,type,taille); \
}

#define ReallocMemoryPrint(var,adr,type,taille) \
{ \
  PrintMemory(var,type,taille,1); \
  ReallocMemory(adr,type,taille); \
}

#define PrintTotalMemory() \
{ \
  if (total_memory >= 1000000000) \
    LogPrintf("total de la memoire allouee : %.2f Go\n", (float)total_memory/1000000000) \
  else if (total_memory >= 1000000) \
    LogPrintf("total de la memoire allouee : %.2f Mo\n", (float)total_memory/1000000) \
  else \
    LogPrintf("total de la memoire allouee : %" PRId64 "\n", total_memory) \
}

#define FreeMemory(adr,type,taille) \
{ \
  if (adr == NULL){ \
    ErrPrintf("Error in FreeMemory\n"); \
    exit(-1);         \
  }           \
  free(adr); \
  total_memory -= sizeof(type)*taille;          \
}

#define Round(x) ((x)>=0 ? (int)((x)+0.5) : (int)((x)-0.5))

#define Calcule_xyz(index, x, y, z) \
{                           \
  int32_t res;                  \
  z = (int)index/HL;        \
  res = index - z*HL;       \
  y = (int)res/L;           \
  x = res - y*L;            \
}

#define Rotate_xz(x, z, angle) \
{                               \
  float alpha = angle*PI/180.0; \
  float co = cos(alpha);        \
  float si = sin(alpha);        \
  float dx0 = L/2 - 0.5;        \
  float dz0 = D/2 - 0.5;        \
  float dx = x - dx0;           \
  float dz = z - dz0;           \
  x = Round(dx0 + co*dx - si*dz); \
  z = Round(dz0 + si*dx + co*dz); \
}

#define Align4(n) n = (n+3) & 0xfffffffc

#define Check_min_max(_Val,_Min,_Max) _Val = ((_Val<_Min)? _Min : (_Val>_Max)? _Max : _Val)

#define Min(_A,_B) ((_A<=_B)? _A : _B)

#define Max(_A,_B) ((_A>=_B)? _A : _B)

#define sgn(x) ((x<0)?-1:((x>0)?1:0))

#define CheckAlign32(_addr) (((long)(_addr) & 0x1f)==0)

#ifdef GTK1
#define g_timeout_add_seconds(interval, function, data) \
  gtk_timeout_add(interval*1000, function, data)

#define g_timeout_add_seconds_full(pr, interval, function, data, notify) \
  gtk_timeout_add_full(pr, interval*1000, function, data, notify)
#endif

/*
#ifdef __CYGWIN__
# include <limits.h>
# define drand48() ((double)rand()/INT_MAX)
# define srand48(sv) (srand((unsigned)(sv)))
#endif
*/
