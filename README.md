# rescal-snow <a name="introduction"></a>
Simulating snow self-organization with cellular automata

![](docs/example_images/snowfall_example.gif)

1. [Introduction](#introduction)
2. [Getting started](#starting)
    1. [Prerequisites](#Prerequisites)
    2. [Dependencies](#Dependencies)
    3. [Installation](#Installation)
    4. [Example 1: a snow cone](#test-cone)
2. [Controlling the simulations](#modifying)
    1. [Example 2: sintering snow](#test-sinter)
    2. [Example 3: dune growth by snowfall](#test-snowfall)
    3. [Visualizing the output](#visualizing)
3. [Setting up parallel runs](#parallel)
    1. [Example 4: parameter space exploration](#test-parallel)
4. [Community guidelines](#community)
    1. [Citation](#Citation)
    2. [Support](#Support)
    3. [Contributing](#contributing)
5. [References and further reading](#references)
6. [Contributors](#authors)
7. [License](#License)

### Background

When wind blows over snow, it self-organizes. This forms surface features, such as ripples and dunes, that alter the reflectivity and thermal conductivity of the snow.

![](docs/example_images/field_examples.png)

These features have just begun to be studied by the snow and climate science communities
(see [1](https://doi.org/10.1002/2015JF003529), [2](https://doi.org/10.5194/tc-13-1267-2019), [3](https://doi.org/10.5194/tc-2019-45) for recent work). 

We created rescal-snow to provide a highly capable snow dune modelling toolkit, to enable snow scientists to study snow features in controlled numerical experiments, and produce high-quality quantitative output. 
We hope that this model will be useful to researchers in snow science, geomorphology, and polar climate.

### Features

 - Snow/sand grain erosion and deposition by wind
 - Snowfall
 - Time-dependent cohesion (snow sintering)
 - Avalanches of loose grains

Rescal-snow can model the growth of ripples from a flat bed of snow; the accumulation of dunes during snowfall events; and the solidification of dunes and snow-waves by sintering.



## Getting started <a name="starting"></a>

### Prerequisites

We assume you have reasonable familiarity with bash and terminal commands.
If you have never used bash, we recommend you stop and work through a short tutorial.
(Our favorite is ['The Unix Shell' from Software Carpentry](http://swcarpentry.github.io/shell-novice/).)
If you modify rescal-snow, you will need to modify and compile C code. We have also included some setup and analysis tools (used in Example 5) written in Python.



### Dependencies

 * C compiler (GCC and CLANG are known to work)
 * cmake>=3.1 (used for compiling)
 * make (used for compiling)
 * Optional packages used for analysis (see the [analysis](analysis) and [scripts/utilities](scripts/utilities) directories):
   * Python3 (used for analysis)
   * numpy (used by Python3 for analysis)
   * pandas (used by Python3 for analysis)
   * scipy (used by Python3 for analysis)

On a Debian-based/Ubuntu Linux machine, the dependencies can be acquired using: 

```bash
sudo apt install gcc cmake make python3 python3-numpy python3-pandas python3-scipy
```

On most machines, the Python packages can also be acquired using:

```bash
pip3 install numpy pandas scipy
```


### Download
Download rescal-snow by cloning this repository with git:

```bash
git clone https://github.com/kellykochanski/rescal-snow.git
cd rescal-snow
```

You may also download the repository manually from [Github](https://github.com/kellykochanski/rescal-snow).

### Installation

These instructions will get rescal-snow running on most Linux environments. For additional installation options, tips on avoiding/installing missing dependencies, and MacOS installation instructions see [docs/how_to_install.md](how_to_install.md).

In a terminal, navigate into the main rescal-snow directory (shown above). Run:
```bash
mkdir build
cd build
cmake -Wno-dev -DCMAKE_BUILD_TYPE=Release .. #Debug can be used instead of Release
make -j 4                                    #Adjust to the number of cores you have for a speedy build
```

### Example 1: a snow cone <a name="test-cone"><a>

The first test run simulates a conical dune shaped by a strong wind. To run the example:
```bash
cd scripts
./snow_cone.run
```
This command will take a few minutes to run, and produces terminal output like:
```bash
>> PAR_FILE=snow_cone.par
>> OMP_NUM_THREADS=1
```
> *Nothing happening?*
>
> This is mostly a good sign - if something is wrong, rescal-snow almost always crashes in the first few seconds. The simulations, however, are computationally expensive, and may take an unreasonably long time to run on slower machines. 
> If this is the case for you, we've put some tips for improving performance in [docs/performance_and_parallelization.md](docs/performance_and_parallelization.md).
> Fortunately, the substantive output from rescal-snow is saved at regular intervals. You can view intermediate output while the simulation is still running, and you can stop the simulation (usually `Ctrl-C`) without losing that output. 

The easiest way to examine the output is to look at rescal-snow's natively generated png files:

```bash
eog scripts/*.png
```
If `eog` is not available, use `open` (MacOS) or any image viewer.

These files are numbered sequentially. The files are saved at regular intervals set by the `-dpng` flag in [scripts/snow_cone.run](scripts/snow_cone.run). For this run, we used `-dpng 10t0`, so file 0 is saved at simulation time t0=0, file 1 (`SNO00001_t0.png`) is saved at
simulation time t0=10, and so on.

|Initial condition, SNO00000_t0.png         |  SNO00003_t0.png                          | SNO00009_t0.png                           |
|-------------------------------------------|-------------------------------------------|-------------------------------------------|
| ![](docs/example_images/snow_cone/00.png) | ![](docs/example_images/snow_cone/03.png) | ![](docs/example_images/snow_cone/09.png) |

Each of the three images above shows a shaded top-down view of a dune (top left), cross-sections through the dune, along the dashed lines (middle left, top right), and a cross-section showing the pressure intensity in the fluid (bottom left).

The initial condition (`t0=0`) is a small cone of sand. As the simulation runs, the cone quickly self-organizes into an elongate dune, and moves downwind.
As time progresses, the dune gradually loses grains. These are not replaced, and the dune dwindles away: this dune is not stable in these (high-wind, no resupply) conditions.

## Controlling the simulation <a name="modifying"></a>

In the next two examples, we walk through the `.run` and `.par` scripts that control the behavior of the simulation, and present example simulations that show you how to modify two important snow parameters: sintering and snowfall.

### Example 2: sintering snow <a name="test-sinter"></a>

The [first example dune](#test-cone) we discussed disappeared quickly as grains blew away in the wind.
Real snow, however, hardens over time. This allows real dunes to persist for days, or even months.
We mimic this process in rescal-snow by adding sintered grains, which erode less easily than regular grains.

You can run an example snow feature with sintering with:
```bash
cd scripts
./sintering.run
```

This produces the following series of png images, with sintered grains shown in light purple and non-sintered grains shown in tan:

![](docs/example_images/sintering/sintering.gif)

This image shows a triangular wave of snow. It quickly leans over to the right. Like the first dune, it begins to blow away in the wind. Some grains at the bottom of the dune, however, sinter (turning light purple in the top left panel of the image). These grains persist for
much longer than the non-sintered grains, but eventually they too blow away.

*Why is this run different from the first one?*

The only change between `snowfall.run` and the previous script, `snow_cone.run` (excepting a few changes to the comments) is the parameter file:
```bash
diff scripts/sintering.run scripts/snow_cone.run
>> < PAR_FILE="sintering.par"
>> > PAR_FILE="snow_cone.par"
```

If we look at the differences between parameter files, we see:
```bash
diff scripts/sintering.par scripts/snow_cone.par
>> < Csp_template = WAVE(15)
>> > Csp_template = CONE(20,40,50)
>> < Lambda_S = 0.01
>> > Lambda_S = 0.00
>> < Lambda_F = 0.05
>> > Lambda_F = 0
```
*What do these changes do?*
 - `Csp_template` selects one of several initial conditions defined in [src/genesis.c](src/genesis.c): the initial condition is now a 15 cell high triangular wave instead of a 20 cell high cone
 - `Lambda_S` controls the rate of sintering: we increased it from 0 to 0.01/t0.
 - `Lambda_F` controls the relative erodibility of the sintered grains: we set it to 0.05/t0, or 5% of the erodibility of the non-sintered grains

**All rescal-snow parameters are all given brief descriptions in the .par file, and in [docs/rescal-snow-inputs.md](docs/rescal-snow-inputs.md).**

### Example 3: dune growth by snowfall <a name="test-snowfall"></a>

Sintering increases the lifespan of snow features, but they still blow away eventually. In reality, snow features don't just disappear - they also appear, and grow. This can occur when the dunes gain mass, for example from snowfall.

We have set up a snowfall simulation in the scripts directory that you can run with:
```bash
cd scripts
./snowfall.run
```
Again, this run script calls a new parameter file:
```bash
diff scripts/snowfall.run snow_cone.run
>> > PAR_FILE="snowfall.par"
>> < PAR_FILE="snow_cone.par"
```

The parameter file for this run has numerous differences from the ones before. The critical ones are:
```bash
diff snowfall.run snow_cone.run
>> > Csp_template = SNOWFALL(4)
>> < Csp_template = CONE(20,40,50)
>> > Lambda_I = 0.001
>> < Lambda_I = 0
```

 - The SNOWFALL template sets an initial condition to be a flat layer of cells of thickness 4, and creates a layer of injection cells to generate snowfall on the simulation ceiling.
 - `Lambda_I` controls the rate of snow injection; its behavior depends on the type and location of the injection cells, and thus on the template.

### Visualizing the simulation output <a name="visualizing"></a>

The snowfall example, above, produces png files among its outputs. Unfortunately, the default rescal-snow rendering does not capture the behavior of these dunes as well as it captured the cone; the falling snow obscures the surface (left-most picture in image below; airborne grains are red).

To see the output clearly, we will make custom images with matplotlib (Python). We will make topographic maps showing the dunes' elevation. This data is output in `scripts/out/ALTIxxxxx.log` files.

We have included a script to recolor three example output frames from the snowfall simulation:

```bash
cd docs/example_images/snowfall
python recolor.py
```
This produces the following images:

| Default image, SNO00050_t0.png | ALTI00050_recolored.png | ALTI00100_recolored.png | ALTI00150_recolored.png |
|--------------------------------|-------------------------|-------------------------|-------------------------|
| ![](docs/example_images/snowfall/SNO00050_t0.png) | ![](docs/example_images/snowfall/ALTI00050_t0_recolored.png) | ![](docs/example_images/snowfall/ALTI00100_t0_recolored.png) | ![](docs/example_images/snowfall/ALTI00150_t0_recolored.png) |

The full evolution of this simulation over a longer time period is shown in the gif at the top of [README.md](README.md).

> Additional scripts for analysing and visualizing the runs are available in the scripts/utilities and analysis directories.

## Setting up parallel runs <a name="parallel"></a>

We believe that building robust, trustworthy models is much simpler when it's easy to make many model runs. This enables:
 - Parameter space exploration
 - Sensitivity analyses
 - Uncertainty quantification
and a general ability to run the model often enough to trust that our results are robust, reproducible, and difficult to skew by cherry-picking.

We have therefore added utilities to help you set up batches of many runs of rescal-snow.

### Example 4: parameter space exploration <a name="test-parallel"></a>

This example requires access to parallel computing resources, such as a university computing cluster. 

> If you do not have access to a computing cluster, the Community Surface Dynamics Modelling System (CSDMS) organization provides free high-performance computing resources for Earth surfaces research.
See [csdms.colorado.edu/wiki/HPC](https://csdms.colorado.edu/wiki/HPC) for details.


In this example, we're going to run 10 instances of rescal. Begin by downloading and installing rescal-snow on your computing cluster, and running one of the above examples to test the installation.

> *Why are we running 10 instances of rescal-snow, not one instance on 10 cores?*
> Rescal-snow doesn't parallelize well: see [docs/performance_and_parallelization.md](docs/performance_and_parallelization.md)

*Writing large numbers of input files*
[scripts/utilities/rescal_utilities.py](scripts/utilities/rescal_utilities.py) contains tools for writing .run and .par scripts automatically. This is much easier, and more bug-free, than writing 10 parameter files by hand.

To set up the runs:

```bash
cd scripts/utilities
python param_space_exploration_example.py
cd ../..
```

You should see that this script has created a new directory called test_parallel_runs, containing ten subdirectories:
```bash
ls test_parallel_runs
>> tauMin0_lambdaI0.001  tauMin1000_lambdaI0.001  tauMin100_lambdaI0.001  tauMin200_lambdaI0.001  tauMin300_lambdaI0.001
>> tauMin0_lambdaI0.01   tauMin1000_lambdaI0.01   tauMin100_lambdaI0.01   tauMin200_lambdaI0.01   tauMin300_lambdaI0.01
```
Each subdirectory contains the executables and input needed for a rescal-snow run:
```bash
ls test_parallel_runs/tauMin0_lambdaI0.001
>> genesis  real_data  rescal  run.par  run.run
```
Running `./run.run` in this `test_parallel_runs/tauMin0_lambdaI0.001` directory would begin a rescal-snow run (you may wish to test this now).

The directory names list the run parameters. For example, `test_parallel_runs/tauMin0_lambdaI0.001` uses a snowfall rate of 0.001/t0, and has the highest possible wind strength (tauMin=0).
The other nine subdirectories contain material for runs with different snowfall rates and wind strengths.

We can submit all of these runs simultaneously using scripts/test_parallel_runs.msub (for Moab systems) or scripts/test_parallel_runs.sbatch (for slurm).
These scripts contain #MSUB or #SBATCH commands describing the user's email, the queue, the resources in terms of nodes or processors, etc.
You will need to modify these commands to match the structure of your computing cluster, or an example job script that works with your account on your computing cluster.

Once you have a working `msub/sbatch` script, run:

```bash
cd scripts
msub test_parallel_runs.msub
```
on a Moab system, or `cd scripts; sbatch test_parallel_runs.sbatch` on slurm.

Each subdirectory should now produce its own `out` directory accompanied by log files, png output, and cellspace (`.csp`) files. After four output intervals pass, it looks like this:

```bash
ls test_parallel_runs/tauMin0_lambdaI0.01
>> DUN.csp             SNO00001_t0.csp.gz    SNO00003_t0.csp.gz  genesis    rescal-ui.xml
>> GENESIS.log	       SNO00001_t0.png	     SNO00003_t0.png     out        run.par
>> SNO00000_t0.csp.gz  SNO00002_t0.csp.gz    SNO00004_t0.csp.gz  real_data  run.run
>> SNO00000_t0.png     SNO00002_t0.png       SNO00004_t0.png     rescal

ls test_parallel_runs/tauMin0_lambdaI0.01/out
>> ALTI00000_t0.log  CELL.log                 LGCA.log      TRANSITIONS.log
>> ALTI00001_t0.log  CELLSPACE_SIGNATURE.log  MVT_IO.log    VEL.log
>> ALTI00002_t0.log  CGV_COEF.log             PROB_CGV.log
>> ALTI00003_t0.log  DENSITE.log              SIGN_HPP.log
>> ALTI00004_t0.log  DOUBLETS.log             TIME.log
```

The following phase diagram shows the 100th images produced by each of these runs (`t0 = 1000; snowfall rate = lambda_I; wind strength = Tau_min`). 

![snowfall-wind phase diagram](docs/example_images/phase_space_exploration/phase_diagram1.png)

The runs with the higher snowfall rate have a much deeper average snow depth than the runs with lower snowfall rate.
The runs with high wind speed (low `Tau_min`) have less even snow cover, with better-defined ripples and dunes.

>*Parallel run errors?*
>
> - msub/sbatch commands not allowed by your computing cluster -> seek support from someone familiar with the cluster
> - permissions errors -> run `chmod u+rwx \*` in test_parallel_runs; contact administrator if this is disallowed
> - script could not find input (run.run or run.par or sealevel_snow.prop) -> rerun `param_space_exploration_example.py` and confirm that it produced the output above; check that you're running `msub/sbatch` from scripts; check relative directory references in msub/sbatch script and submit.sh
> - runs time out before producing useful output -> increase walltime; large simulations on slow machines may ultimately take a few hours



## Community guidelines <a name="community"></a>

We are excited to hear about your scientific work with ReSCAL.
This section lists the best way to bring your project - and its successes, challenges, bugs, and developments - to our attention.

We encourage you to interact with the project through Github (see below). This will allow easy integration of your changes and prevent rescal-snow from fragmenting excessively.


### Citation

Do you want to incentivize developers to build and maintain the software you need?
Cite us!

We are working on getting a doi for this repository through the Journal of Open Source Software (JOSS). 

This software inherits much of its functionality from the Real-Space Cellular Automaton Laboratory, ReSCAL. Please credit those developers by citing:
 - ['A real-space cellular automaton laboratory', O Rozier and C Narteau, Earth Surface Processes and Landforms 39(1) 98-109, 2013, doi=10.1002/esp.3479](https://onlinelibrary.wiley.com/doi/abs/10.1002/esp.3479)


### Support

If you have challenges or questions, look at the material under 'further information' or reach out to us.

Issues should be reported using Github's issue tracking function on this repository, [here](https://github.com/kellykochanski/rescal-snow/issues).

Issues which cannot be handled via Github can be addressed to

    Kelly Kochanski
    kelly.kochanski@gmail.com
    www.github.com/kellykochanski


### Contributing

We have built a model of the basic growth and function of snow dunes, but expect that many users may want more detailed features.

If you wish to contribute a new feature, we recommend you fork our repository, commit your changes to a new feature or development branch, and create a pull request. An example contribution workflow, with git instructions, is outlined by the [LAVA software community project contribution
guide](https://docs.lavasoftware.org/lava/contribution.html)

Rescal-snow is distributed under the GNU GPL 3.0 license; all contributions must be made under this license or a later version.

## References and further reading <a name="references"></a>

The [docs](docs) folder contains additional information on 

 - [Alternate installations](docs/how_to_install.md), 
 - [Performance and parallelization issues](docs/performance_and_parallelization.md), 
 - [Model inputs and configuration](docs/rescal-snow-inputs.md),
 - [Model calibration and validation](docs/calibration_and_validation.md),
 - [Development history](docs/NEWS.md)

For more background on snow dunes, sintering, and self-organization, see:
 - ['Snow bedforms: A review, new data, and a formation model', Filhol and Sturm, 2015](https://doi.org/10.1002/2015JF003529)
 - ['The evolution of snow bedforms in the Colorado Front Range', Kochanski, Anderson and Tucker, 2019](https://doi.org/10.5194/tc-13-1267-2019)
 - ['Statistical classification of self-organized snow surfaces', Kochanski, Anderson and Tucker, 2018](https://doi.org/10.1029/2018GL077616)
 - ['Studies on interaction between wind and dry snow surface', Kobayashi, 1980](https://eprints.lib.hokudai.ac.jp/dspace/bitstream/2115/20242/1/A29_p1-64.pdf)

For models of snow dunes based on ReSCAL, see:
 - ['Understanding snow bedform formation by adding sintering to a cellular automata model', Sharma, Braud and Lehning, preprint](https://doi.org/10.5194/tc-2019-45)

For more information about the initial ReSCAL development and the backend function of the cellular automaton and lattice gas model, see:
 - ['Setting the length and timescales of a cellular automaton dune model', Narteau et al., 2009](https://doi.org/10.1029/2008JF001127)
 - ['A real-space cellular automaton laboratory', Rozier and Narteau, 2014](dx.doi.org/10.1002/esp.3479)
 - ['Transport capacity and saturation mechanism in a real-space cellular automaton dune model', Gao et al., 2014](dx.doi.org/10.5194/adgeo-37-47-2014)
To learn the underlying principles of the lattice gas cellular automaton (LGCA)  model (recommended before modifying the LGCA, the boundary conditions, or the aspect ratio of the simulation) see:
 - ['Lattice-gas automata for the Navier-Stokes equation', Frisch, Hasslacher and Pomeau, 1986'](https://doi.org/10.1103/PhysRevLett.56.1505)

## Contributors <a name="authors"></a>
See [AUTHORS.md](AUTHORS.md).

## License
GNU GPL 3.0 or any later version. See [docs/LICENSE.md](docs/LICENSE.md).

SPDX-License-Identifier: GPL-3.0-or-later

Release: LLNL-CODE-785837. See [NOTICE](NOTICE) for details.