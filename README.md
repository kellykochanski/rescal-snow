# Introduction to rescal-snow <a name="introduction"></a>

##Table of Contents
1. [Introduction](#introduction)
2. [Getting started](#starting)
    1. [Prerequisites](#Prerequisites)
    2. [Dependencies](#Dependencies)
    3. [Installation](#Installation)
    4. [Example 1: a snow cone](#test1)
    5. [Example 2: 5 snow cones](#test2)
3. [Controlling the simulations](#modifying)
    1. [Example 3: dune growth by snowfall](#test3)
    2. [Example 4: dune arrest by sintering](#test4)
4. [Setting up parallel runs](#parallel)
    1. [Example 5: parameter space exploration](#test5)
5. [Community guidelines](#community)
    1. [Citation](#Citation)
    2. [Support](#Support)
    3. [Reporting issues](#issues)
    4. [Contributing](#contributing)
6. [References and further reading](#references)
7. [Authors and contributors](#authors)

Snow bedforms are some of the most widespread surface textures on Earth, covering Antarctica, Greenland, and much tundra and sea ice around the poles.
They alter the albedo and thermal conductivity of snow surfaces, and thereby the climates of the polar regions and the world.

rescal-snow is designed to simulate the growth of snow dunes and ripples. It is adapted from a community-standard sand dune simulation, ReSCAL, that uses cellular automata to model the behavior of large dunes as the sum of interactions between small grains.

rescal-snow can model the growth of ripples from a flat bed of snow; the accumulation of dunes during snowfall events; and the solidification of dunes and snow-waves by sintering.

### Why use rescal-snow?

## Getting started <a name="starting"></a>

These instructions will get rescal-snow running on most linux environments; for additional installation options, tips on avoiding/installing missing dependencies, and MacOS installation instructions see `how_to_install.txt'.
If you would like to install on windows, we recommend you don't. It will almost certainly be easier to get access to a linux virtual machine (for small runs) or computing cluster (for larger runs) than to sort out the dependencies.

### Prerequisites

This build relies on autoconf, libpng, and gcc. It was tested on Ubuntu 18.04 with gcc 6.9. See `how_to_install.txt' for tips for installing autoconf or using intel compilers.

We assume you have reasonable familiarity with bash and terminal commands.
If you have never used bash, we recommend you stop and work through a short tutorial.
(Our favorite is 'The Unix Shell' from Software Carpentry: see http://swcarpentry.github.io/shell-novice/ )

### Dependencies
rescal-snow requires glib-2.0, zlib-1.6, libpng, and pthread, as well as a C compiler. We have included compatible versions of libpng and zlib in the 'lib' directory. The included configure script will automatically build them and point rescal toward those directories. 

Many of the auxiliary tools (see the 'analysis' and 'scripts/utilities' directories) are written in Python. These are written for Python3 and rely on libraries os, sys, numpy, scipy, csv and pandas; this is all packaged with the Anaconda Python distribution.

### Installation

In a terminal, navigate into the main rescal-snow directory (the one containing this readme, as well as 'scripts', 'src', etc'). Run:
  ./configure
  make
  
If this doesn't work, see 'how_to_install.txt'.

Unless stated otherwise, every command in this README starts from the top directory.

### Example 1: a snow cone <a name="test1"><a>

The first test run simulates a field of conical dunes under a strong wind. To run the example:
```bash
cd scripts
./snow_cone.run
```

This command may take several minutes to run, but it will produce output at intervals of a few seconds.

The easiest way to examine the output is to look at rescal-snow's natively generated png files.

```bash
eog scripts/*.png
```
on most Linux systems. If `eog` is not available, use `open` (MacOS) or any image viewer. This will display a series of png, including the examples below:

|Initial condition, SNO00000_t0.png 	|  SNO00003_t0.png   | SNO00009_t0.png  |
|------------------------|--------------------|------------------|
| ![](docs/example_images/snow_cone/00.png) | ![](docs/example_images/snow_cone/03.png) | ![](docs/example_images/snow_cone/09.png) |

These images show a shaded top-down view of a dune (top left), cross-sections through the dune, along the dashed lines (middle left, top right), and a cross-section showing the pressure intensity in the fluid (bottom left).

Here, the initial condition is a small cone of sand. As the simulation runs, the cone quickly self-organizes into an elongate dune, and moves downwind.
As time progresses, the dune gradually loses grains. These are not replaced, and the dune dwindles away.
Your output may not match the example images precisely due to deliberate stochastic behavior in rescal-snow, but the overall pattern should be the same.

### Example 2: 5 snow cones <a name="test2"></a>

This test run uses a different size and initial condition of the simulation. The simulation has a larger domain, and is likely to run slowly: this is a good test of the performance of rescal-snow on your machine.

```bash
rm scripts/out/*
./snow_cone_x5.run
```
and view as above.



## Modifying the simulation <a name="Modifications">

To use rescal-snow for scientific projects, you will likely wish to modify its behaviors.
In the above scripts, we used two pre-generated run files (snow_cone.run and snow_cone_x5.run). This section will show you how to modify the run files and model parameters yourself.
Here, we focus on the two most important snow parameters: snowfall and sintering.
Other parameters are described in more detail in docs/rescal-snow-input.md.

### Test run 3: dune growth by snowfall <a name="test3"></a>

### Example 4: dune arrest by sintering <a name="test4"></a>

## Setting up parallel runs

We believe that building robust, trustworthy models is much simpler when it's easy to make many model runs. This enables:
 - Parameter space exploration
 - Sensitivity analyses
 - Uncertainty quantification
and a general ability to run the model often enough to trust that our results are robust, reproducible, and difficult to skew by cherry-picking.

We have therefore added utilities that make it easy to run batches of dozens of runs of rescal-snow.

ReSCAL v1.6 is technically configured to run in parallel (see the OPENMP flag in src/defs.h; this allows the lattice gas the the cellular automaton to run on separate processors). 
We have not emphasized this feature in rescal-snow because we have not been able to achieve satisfying parallel efficiency; thus far, we have found it more useful to perform larger numbers of serial runs.

### Example 5: parameter space exploration

This test presumes you have access to parallel computing resources, such as a university computing cluster or a supercomputer. 
If you do not have access to a computing cluster, the Community Surface Dynamics Modelling System (CSDMS) organization provides free high-performance computing resources for Earth surfaces research.
See csdms.colorado.edu/wiki/HPC for details and to apply for an account.

For this run, we're going to run 10 instances of rescal. 
Rather than writing 10 parameter files by hand, we're going to use the python scripts in the scriptsutilities directory to generate them automaticlaly.
To set up the run, from the main rescal-snow directory, run: 

```bash
cd scripts/utilities
python param_space_exploration_example.py
cd ../..
```

You should see that this script has created a new directory called test_parallel_runs, containing ten subdirectories:

    ls test_parallel_runs
    >> tauMin0_lambdaI0.001  tauMin1000_lambdaI0.001  tauMin100_lambdaI0.001  tauMin200_lambdaI0.001  tauMin300_lambdaI0.001
    >> tauMin0_lambdaI0.01   tauMin1000_lambdaI0.01   tauMin100_lambdaI0.01   tauMin200_lambdaI0.01   tauMin300_lambdaI0.01
 
And each subdirectory contains the executables and input needed for a rescal-snow run:

    ls test_parallel_runs/tauMin0_lambdaI0.001
    >> genesis  real_data  rescal  run.par  run.run

Running `./run.run` in this test_parallel_runs/tauMin0_lambdaI0.001 directory would begin a rescal-snow run, using a snowfall rate (controlled by LambdaI) of 0.001 cells per unit time, with a threshold wind shear (controlled by tau Min) of 0.
Each of the nine other subdirectories is similar, but with different values of tauMin and LambdaI as expressed in the directory names.

We can submit all of these runs simultaneously using scripts/test_parallel_runs.msub (for Moab systems) or scripts/test_parallel_runs.sbatch (for slurm).
Note that both of these scripts contain many #MSUB or #SBATCH commands describing the user's email, the queue, the resources in terms of nodes or processors, etc.
These will need to be modified to match the structure of your computing cluster. Look for an example job script for your cluster.
Once the msub/sbatch script is modified, run:

```bash
cd scripts
msub test_parallel_runs.msub
```
on a Moab system, or `sbatch test_parallel_runs.sbatch` on slurm.

If the run is successful, each subdirectory will produce its own 'out' directory accompanied by log files, png output, and cellspace (.csp) files.

*Parallel troubleshooting*

Common reasons for parallel run failures include:
 - msub/sbatch commands not allowed by your computing cluster -> seek support from someone familiar with the cluster
 - permissions errors -> run `chmod u+rwx *` in test_parallel_runs; contact administrator if this is disallowed
 - script could not find run.run/run.par/sealevel_snow.prop -> rerun param_space_exploration_example.py and confirm that it produced the output above; check that you're running msub/sbatch from scripts; check relative directory references in msub/sbatch script and submit.sh
 - runs time out before producing useful output -> increase walltime; for large runs, expect hours of calculation

## Community guidelines <a name="community"></a>

We are excited to hear about your scientific work with ReSCAL.
This section lists the best way to bring your project --- and its successes, challenges, bugs, and developments --- to our attention.

We encourage you to interact with the project through github (see below) to encourage easy integration of your changes and issues into the major project.
However, we check our email more often than our github, and will respond faster if git issues are accompanied by emails.

### Citation

Do you want to incentivize developers to build more of the kinds of software you need?
Cite us!

We are working on getting a doi for this repository through the Journal of Open Source Software (JOSS). If you use this software in the meantime, please contact the developers for the most up-to-date reference.

This software inherits much of its functionality from the Real-Space Cellular Automaton Laboratory, ReSCAL. Please credit those developers by citing:
 - 'A real-space cellular automaton laboratory', O Rozier and C Narteau, Earth Surface Processes and Landforms 39(1) 98-109, 2013, doi=10.1002/esp.3479

### Support

If you have challenges or questions, look at the material under 'further information' or reach out to us.

Primary contact:
Kelly Kochanski
kelly.kochanski@gmail.com
www.github.com/kellykochanski

Note that questions pertaining to the backend function of the cellular automaton and lattice gas may be best addressed by the developers of ReSCAL v1.6.
ReSCAL v1.6 contact:
Clement Narteau
narteau@ipgp.fr

### Reporting issues <a name="issues"></a>

Issues occur when the software does not behave as stated in the documentation.
Issues may be reported using github's issue tracking function on this repository, www.github.com/kellykochanski/rescal-snow.

### Contributing

We have built a model of the basic growth and function of snow dunes, but expect that many users may want more detailed features.

If you wish to contribute a new feature, we recommend you fork our repository, commit your changes to a new feature or development branch, and create a pull request. An example contribution workflow, with git instructions, is outlined by the [LAVA software community project contribution
guide](https://docs.lavasoftware.org/lava/contribution.html)

## References and further reading <a name="references"></a>

For more information about the initial ReSCAL development and the backend function of the cellular automaton and lattice gas model, see:
 - ['Setting the length and timescales of a cellular automaton dune model', Narteau et al., 2009](https://doi.org/10.1029/2008JF001127)
 - ['A real-space cellular automaton laboratory', Rozier and Narteau, 2014](dx.doi.org/10.1002/esp.3479)
 - ['Transport capacity and saturation mechanism in a real-space cellular automaton dune model', Gao et al., 2014](dx.doi.org/10.5194/adgeo-37-47-2014)
To learn the underlying principles of the lattice gas cellular automaton (LGCA)  model (recommended before modifying the LGCA, the boundary conditions, or the aspect ratio of the simulation) see:
 - ['Lattice-gas automata for the Navier-Stokes equation', Frisch, Hasslacher and Pomeau, 1986'](https://doi.org/10.1103/PhysRevLett.56.1505)

## Authors and contributors <a name="authors"></a>
See AUTHORS.md for a full list of contributors.
