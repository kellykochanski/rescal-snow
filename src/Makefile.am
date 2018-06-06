bin_PROGRAMS = rescal genesis regenesis bin2png csp2png cspinfo csp2bin

VER_GTK = 2

AM_CFLAGS = -I/opt/local/include \
            -std=gnu99 -DGTK$(VER_GTK) $(GTK_CFLAGS) \
            -DRESCAL_DATA_DIR=\""$(PWD)"\" \
            -DHAVE_CONFIG_H \
            $(AVX_CFLAGS) \
	    $(OPENMP_CFLAGS) 

AM_LDFLAGS = -L/opt/local/lib 

rescal_LDADD = -lm $(GTK_LIBS) $(GMODULE_LIBS) $(GTHREAD_LIBS)
genesis_LDADD = -lm
regenesis_LDADD = -lm
bin2png_LDADD = -lm
csp2png_LDADD = -lm


common_SOURCES = \
  cells.c \
  cells.h \
  lgca.c \
  lgca.h \
  defs.h \
  doublets.c \
  doublets.h \
  format.c \
  format.h \
  macros.h \
  models.c \
  models.h \
  param.c \
  param.h \
  simul.c \
  simul.h \
  surface.c \
  surface.h \
  space.c \
  space.h \
  trace.c \
  trace.h \
  transitions.c \
  transitions.h \
  view.c \
  view.h


rescal_SOURCES = \
  $(common_SOURCES) \
  callbacks.c \
  callbacks.h \
  drand48.c \
  entry.c \
  entry.h \
  rescal.c \
  rescal.h

genesis_SOURCES = \
  defs.h \
  macros.h \
  cells.h \
  format.c \
  format.h \
  param.c \
  param.h \
  genesis.c \
  genesis.h

regenesis_SOURCES = \
  defs.h \
  macros.h \
  cells.h \
  format.c \
  format.h \
  regenesis.c
	
bin2png_SOURCES = \
  $(common_SOURCES) \
  bin2png.c

csp2png_SOURCES = \
  $(common_SOURCES) \
  csp2png.c

cspinfo_SOURCES = \
  format.c \
  format.h \
  cells.h \
  macros.h \
  cspinfo.c

csp2bin_SOURCES = \
  defs.h \
  format.c \
  format.h \
  cells.h \
  macros.h \
  csp2bin.c

EXTRA_DIST = \
  rescal-ui.xml \
  makefile.default \
  makefile.icc

