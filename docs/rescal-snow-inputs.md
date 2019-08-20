
#Inputs to rescal-snow

Inputs controlling the behavior of rescal-snow simulations

Author: Kelly Kochanski, kelly.kochanski@gmail.com, Jun 12 2019

## Finding the input files

rescal-snow inputs are controlled in three files:
 - The run script (called run.run, snow_cone.run, etc) controls meta-parameters for the run, like the run length, the locations of input/output files, and the output frequency
 - The parameter file (called run.par, snow_cone.par, etc) controls the behavior of the model during the run, including the rates of the modelled physical processes
 - The real data parameters file (real_data/desert_earth.prop, real_data/sealevel_snow.prop) contains constants useful for converting model output to real dimensions of length and time

We have created a python script containing Parameter_file and Run_script classes, scripts/utilities/rescal_utilities.py, which can automatically generate both parameter files and run scripts. This may be useful for setting up large numbers of runs or avoiding syntax errors.

All three input files are sensitive to changes in whitespace. Sorry.

### The run script
The run script (run.run, snow_cone.run, etc) is a bash script. You can modify it to change the run in the following ways:
 - `./clean` : asks user if they want to remove previous output; use -f flag to skip request (useful for remote/queued/parallel runs)
 - `ln -s ../build/genesis .` : tells the run where to find the genesis executable; must be modified if you change the directory structure
 - `ln -s ../build/rescal .` : tells the run where to find the rescal executable; must be modified if you change the directory structure
 - `PAR_FILE` gives the location of the parameter file (see section below). *Expect to change this whenever you change model physics*
 - `export OMP_NUM_THREADS=1` : sets environment variable OMP_NUM_THREADS, we recommend using 1
 - `GENESIS_LOG_FILE` and `RESCAL_LOG_FILE` : set locations of two log files
 - `nice` : rescal-snow will run more slowly to avoid interfering with other applications on your computer, useful for local runs
 - run flags : 

### The parameter script

This file contains a list of parameters controlling simulation physics.
 - Output_directory : output goes here, useful for data management.
 - Model: rescal-snow can be compiled to allow or disallow certain physics. To ensure the compiled version works with the model you want to run, this must match the compiler flag in src/defs.h. Example: `Model = SNO` here, `#define MODEL_SNO` in src/defs.h.
 - Csp_template : gives a template name and arguments for initial conditions. These are defined and described in src/genesis.c.
 - Phys_prop_file : points to the real data file described in the section below.
 - Simulation size H, L, D : sets size of simulation in cells.
 - Boundary conditions : OPEN, PERIODIC, or REINJECTION. Determines what happens to cells that pass the downwind boundary of the simulation.
 - Lambda_X, Coef_X : control the rates of cell transitions, the main control on model physics. These are best shown graphically: see Fig. 2 of [Narteau et al., 2009](dx.doi.org/10.1029/2008JF001127).

The parameters you are likely to want to modify for snow studies include:
 - Csp_template : the initial conditions, see src/genesis.c
 - Tau_min : the control for wind speed - note that wind speed *decreases* non-linearly as you increase Tau_min
 - Lambda_I : the snowfall rate (must use `Csp_template = SNOWFALL`)
 - Lambda_S : the sintering rate
 - Lambda_F : the hardness of sintered snow, relative to the hardness of freshly-fallen snow
 - Coef_A : the ratio between vertical and horizontal transport of snow in air; controls the particle settling speed

### The real_data folder

All of the information we have given rescal-snow so far has been non-dimensional. This is useful insofar as it makes rescal simulations easy to generalize, but not useful for comparing the results to specific situations. 
This file contains information about the physical situation that rescal-snow represents. For more discussion on converting rescal-snow output to reality, see 'docs/calibration_and_validation.md'.

The file scripts/real_data/sealevel_snow.prop contains the following variables:
 - kappa: the von Karman constant, dimensionless, generally agreed to be near 0.4
 - z_0: the [aerodynamic roughness length](https://en.wikipedia.org/wiki/Roughness_length) of the surface, m
 - rhoair: the density of the fluid that moves the grains, kg/m^3
 - rhosand: the density of the particles being moved by the fluid, kg/m^3
 - g: gravitational acceleration, m/s^2, always close to 9.81 on Earth
 - d: the diameter of the particles, m

Take care while calculating these variables. Their values vary considerably for real snow, sand and air. We have endeavored to give reasonable defaults in sealevel_snow.prop but do not expect these to represent all real situations.
The accuracy of these numbers directly affects the accuracy of the dune sizes, speeds, etc produced by rescal-snow.
