2019 update: 
This is the readme for ReSCAL v1.6, 2016, 
and the examples and dependencies are now outdated.
Refer to README.md for up-to-date information.

ReSCAL: a Real-Space Cellular Automaton Laboratory
--------------------------------------------------
Copyright (C) 2011-2016

(2019 note: this readme describes ReSCAL v1.6 in 2016.
The current software has been modified from this;
see README.md for a description of the current state
of this software)

Description
-----------

This is ReSCAL, a cellular automaton modeling tool in real 3d space, 
primarily designed for geomorphology. However we believe it is 
sufficiently generic in its approach to deal with a wide range of 
applications. 
Basically, it simulates asynchronous first-neighbour interactions 
(see the papers in the References section). We refer to Markov 
chains and Poisson processes for the evolution of the physical time.

We also implemented a number of additional mechanisms that are 
relevant in usual geophysical models
- correlations (coupled processes)
- avalanches (granular processes)
- lattice gas coupling (erosion)
- rotating table (multidirectional winds)
- map surrections (tectonic processes)


Content
-------

ReSCAL package comprises several programs:
- rescal: ReSCAL main application
- genesis: generation of a 3d cellular space (CSP or raw binary format)
- regenesis: modification of a cellular space
- csp2png: PNG image generator from CSP-file
- bin2png: PNG image generator from BIN-file
- cspinfo: CSP metadata reader
- csp2bin: CSP to raw binary conversion tool

Samples are present in the scripts directory, along with Matlab/octave
functions that support CSP format.


Requirements
------------

Needed libraries are glib-2.0, gtk+-2.0, gdk-pixbuf and pthread.
'genesis.py' uses NumPy, but it is possible to use C version of 
genesis instead.


Getting started
---------------

1) Build all programs by typing
./configure
make

(See also 'INSTALL' file for detailled instructions)

2) Go to the scripts directory and launch the 'run' script:
cd scripts
./run

This will open a window, with no menu, displaying a dunes simulation 
controlled by air flow. Note the initial phase of flow stabilization.

The environment variable OMP_NUM_THREADS lets you set the number of
OpenMP threads in the computation of the flow. You can get optimal 
speed by setting this variable according to the number of available 
CPU cores.

(See rescal command line options to customize display preferences
in the script)

3) Edit the 'run' script and change the value of PAR_FILE:
PAR_FILE="dun_rot.par"
Then launch rescal again:
./run

Now you have a similar execution, with rotations occuring regularly,
as if the wind was rotating. Again the flow needs to be stabilized
after each rotation.

4) Edit a parameter file, change some values and see what happens in 
your model.

For example, you can change the overall size by editing the 'H', 'L', 
'D' parameters in 'dun.par'.
H = heigth (vertical)
L = length (from east to west)
D = depth (from north to south)

Beware that a huge cellular space will slow down the simulation.

Another suggestion would be to choose a predefined CSP template for
the initial configuration of the cellular space. For example, you can 
get a conical pile of sand by setting:
Csp_template = RCONE 

5) Save the output files of your last simulation by using the 'xrep' 
script.


Advanced usage
--------------

It is possible to select and/or modify any of the predefined models.
So far, this requires to edit some of the sources in C language. 
Here we briefly explain how to proceed. Ideally, we assume 
that the user have some notion of programming.

1) Edition of genesis sources

The genesis program is provided in python and C languages. These two
standalone versions are respectively located in 'scripts' and 
'src' directories. The python version of genesis is easier to
read, however it does not support the CSP format yet, nor any of the
CSP templates.

Note that for historical reasons, a cellular space is composed of 
vertical plans in the east-west direction. Each plan is a 
concatenation of horizontal lines from left to right. The very first 
cell corresponds to the north-upper-left corner.

When lattice gas is used, a floor and a ceiling of DUM cells are 
required to prevent gas particles from escaping the cellular space.

2) Edition of rescal sources

Start by selecting a model and its specific options, in 'defs.h'. 

You may also edit the transitions in 'models.c'. Basically, the 
syntax looks like

trans(direction, cel1_start, cel2_start, cel1_end, cel2_end, rate).

Then you will get an oriented transition 

[cel1_start, cel2_start] -> [cel1_end, cel2_end] 

occuring at the specified rate. See USER'S GUIDE in 'models.c' for 
an overview of the available functions.

Use the command-line option '-info' to generate a number of log files.
For instance, the list of all transitions is present in 'TRANS.log'.

Feel free to have a look around and play with the existing mechanisms
in 'space.c' and 'surface.c', or to suggest any improvement you would 
like to see. Obviously, there is still a lot left to do ! (See the 
TODO file)

A modified version can be freely distributed provided that you keep 
intact all the notices that refer to the license (See 'COPYING' file).

We would be grateful if you mention the use of this program in a 
related work or publication.


References
----------

[1] On a small scale roughness of the core-mantle boundary
C. Narteau, J.L. Le Mouël, J.P. Poirier, E. Sepulveda and M. Shnirman
Earth and Planetary Science Letters, 191, 49-60 (2001).

[2] Setting the length and time scales of a cellular automaton dune 
model from the analysis of superimposed bedforms
C. Narteau, D. Zhang, O. Rozier and P. Claudin
Journal Geophysical Research, 114, F03006 (2009).

[3] Morphodynamics of barchan and transverse dunes using a cellular 
automaton model
D. Zhang, C. Narteau and O. Rozier
Journal Geophysical Research, 115, F03041 (2010).


========================================================================

ReSCAL is free software; you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation. It is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY WHATSOEVER.

Check the License for details. You should have received a copy of it, along
with the package; see the file 'COPYING'. If not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.


