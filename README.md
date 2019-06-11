***rescal-snow***

Snow bedforms are some of the most widespread surface textures on Earth, covering Antarctica, Greenland, and much tundra and sea ice around the poles.
They alter the albedo and thermal conductivity of snow surfaces, and thereby the climates of the polar regions and the world.

rescal-snow is designed to simulate the growth of snow dunes and ripples. It is adapted from a community-standard sand dune simulation, ReSCAL, that uses cellular automata to model the behavior of large dunes as the sum of interactions between small grains.

rescal-snow can model the growth of ripples from a flat bed of snow; the accumulation of dunes during snowfall events; and the solidification of dunes and snow-waves by sintering.

**Getting started**

These instructions will get rescal-snow running on most linux environments; for additional installation options, tips on avoiding/installing missing dependencies, and MacOS installation instructions see `how_to_install.txt'.
If you would like to install on windows, we recommend you don't. It will almost certainly be easier to get access to a linux virtual machine (for small runs) or computing cluster (for larger runs) than to sort out the dependencies.

*Prerequisites*

This build relies on autoconf, libpng, and gcc. It was tested on Ubuntu 18.04 with gcc 6.9. See `how_to_install.txt' for tips for installing autoconf or using intel compilers.

We assume you have reasonable familiarity with bash and terminal commands.
If you have never used bash, we recommend you stop and work through a short tutorial.
(Our favorite is 'The Unix Shell' from Software Carpentry: see http://swcarpentry.github.io/shell-novice/ )

*Installation*

In a terminal, navigate into the main rescal-snow directory (the one containing this readme, as well as 'scripts', 'src', etc'). Run:
  ./autogen.sh
  ./configure
  make
  
If this doesn't work, see 'how_to_install.txt'.

*Test run 1: cones turn into dunes*

*Test run 2: dune growth by snowfall*

*Test run 3: dune shrinkage by sintering*

**Parallel runs**

We believe that building robust, trustworthy models is much simpler when it's easy to make many model runs. This enables:
 - Parameter space exploration
 - Sensitivity analyses
 - Uncertainty quantification
and a general ability to run the model often enough to trust that our results are robust, reproducible, and difficult to skew by cherry-picking.

ReSCAL v1.6 is technically configured to run in parallel (see the OPENMP flag in src/defs.h; this allows the lattice gas the the cellular automaton to run on separate processors). 
We have not emphasized this feature in rescal-snow because we have not been able to achieve satisfying parallel efficiency; thus far, we have found it more useful to perform larger numbers of serial runs.

*Test run 4: a small parameter space exploration* 

This test presumes you have access to parallel computing resources, such as a university computing cluster or a supercomputer. 
If you do not have access to a computing cluster, the Community Surface Dynamics Modelling System (CSDMS) organization provides free high-performance computing resources for Earth surfaces research.
See csdms.colorado.edu/wiki/HPC for details and to apply for an account.

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

**Developing new simulations**



**Authors and acknowledgements**
See AUTHORS.md for a full list of contributors.

**Contributing**
Email kelly.kochanski@gmail.com for snow-related developments; narteau@ipgp.fr for comments on sand dune simulation or ReSCAL.
