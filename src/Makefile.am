bin_PROGRAMS = $(top_builddir)/scripts/rescal \
               $(top_builddir)/scripts/genesis \
               $(top_builddir)/scripts/regenesis \
               $(top_builddir)/scripts/bin2png \
               $(top_builddir)/scripts/csp2png \
               $(top_builddir)/scripts/cspinfo \
               $(top_builddir)/scripts/csp2bin

VER_GTK = 2

AM_CFLAGS = -std=gnu99 -I/opt/local/include \
            -DRESCAL_DATA_DIR=\""$(PWD)"\" \
            -DHAVE_CONFIG_H \
            $(AVX_CFLAGS) \
	    $(OPENMP_CFLAGS) 

AM_LDFLAGS = -L/opt/local/lib

__top_builddir__scripts_rescal_LDADD = -lm $(GTK_LIBS) $(GMODULE_LIBS) $(GTHREAD_LIBS) ../lib/libpng-1.6.34/.libs/libpng16.a ../lib/zlib-1.2.11/libz.a
__top_builddir__scripts_genesis_LDADD = -lm
__top_builddir__scripts_regenesis_LDADD = -lm
__top_builddir__scripts_bin2png_LDADD = -lm ../lib/libpng-1.6.34/.libs/libpng16.a ../lib/zlib-1.2.11/libz.a
__top_builddir__scripts_csp2png_LDADD = -lm ../lib/libpng-1.6.34/.libs/libpng16.a ../lib/zlib-1.2.11/libz.a

AM_CPPFLAGS = -I$(top_builddir)/lib/libpng-1.6.34/ -I$(top_builddir)/lib/zlib-1.2.11/

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

__top_builddir__scripts_rescal_SOURCES = \
  $(common_SOURCES) \
  callbacks.c \
  callbacks.h \
  drand48.c \
  entry.c \
  entry.h \
  rescal.c \
  rescal.h

__top_builddir__scripts_genesis_SOURCES = \
  defs.h \
  macros.h \
  cells.h \
  format.c \
  format.h \
  param.c \
  param.h \
  genesis.c \
  genesis.h

__top_builddir__scripts_regenesis_SOURCES = \
  defs.h \
  macros.h \
  cells.h \
  format.c \
  format.h \
  regenesis.c

__top_builddir__scripts_bin2png_SOURCES = \
  $(common_SOURCES) \
  bin2png.c

__top_builddir__scripts_csp2png_SOURCES = \
  $(common_SOURCES) \
  csp2png.c

__top_builddir__scripts_cspinfo_SOURCES = \
  format.c \
  format.h \
  cells.h \
  macros.h \
  cspinfo.c

__top_builddir__scripts_csp2bin_SOURCES = \
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

