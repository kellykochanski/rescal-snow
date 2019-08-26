# rescal-snow <a name="introduction"></a>
Simulating the evolution of dunes and snow waves

![](docs/example_images/snowfall_example.gif)

1. [Motivation](#introduction)
2. [Getting started](#starting)
    1. [Prerequisites](#Prerequisites)
    2. [Dependencies](#Dependencies)
    3. [Installation](#Installation)
3. [Features and examples](#examples)
    1. [Example 1: a snow cone](#test-cone)
    2. [Example 2: sintering snow](#test-sinter)
    3. [Example 3: dune growth by snowfall](#test-snowfall)
    4. [Example 4: parameter space exploration](#test-parallel)
4. [References and further reading](#references)
5. [Community guidelines](#community)
    1. [Citation](#Citation)
    2. [Support](#Support)
    3. [Contributing](#contributing)
6. [Contributors](#authors)
7. [License](#License)

### 1. Motivation

When wind blows over snow, it self-organizes. This forms surface features, such as ripples and dunes, that alter the reflectivity and thermal conductivity of the snow.

![](docs/example_images/field_examples.png)

These features have just begun to be studied by the snow and climate science communities
(see [1](https://doi.org/10.1002/2015JF003529), [2](https://doi.org/10.5194/tc-13-1267-2019), [3](https://doi.org/10.5194/tc-2019-45) for recent work). 

We created rescal-snow to provide a highly capable snow dune modelling toolkit, to enable snow scientists to study snow features in controlled numerical experiments, and produce high-quality quantitative output. 
Rescal-snow is able to simulate:

 - Snow/sand grain erosion and deposition by wind
 - Snowfall
 - Time-dependent cohesion (snow sintering)
 - Avalanches of loose grains

Rescal-snow also contains tools for data management and analysis, including:
 - Variable levels of run-time analysis
 - Workflow tools for generating and running many simulations in parallel

We hope that this model will be useful to researchers in snow science, geomorphology, and polar climate.

## 2. Getting started <a name="starting"></a>

### 2.1 Prerequisites

We assume you have reasonable familiarity with bash and terminal commands.
If you have never used bash, we recommend you stop and work through a short tutorial.
(Our favorite is ['The Unix Shell' from Software Carpentry](http://swcarpentry.github.io/shell-novice/).)
If you modify rescal-snow, you will need to modify and compile C code. We have also included some setup and analysis tools (used in Example 5) written in Python.

### 2.2 Dependencies
Rescal-snow requires glib-2.0 and pthread, as well as a C compiler and the [lib](lib) directory, which contains compatible versions of zlib and libpng.

Many of the auxiliary tools (see the [analysis](analysis) and [scripts/utilities](scripts/utilities) directories) are written in Python3. These were developed in an Anaconda environment and require os, sys, numpy, scipy, csv and pandas.

### 2.3 Download
Download rescal-snow by cloning this repository with git:
```bash
git clone https://github.com/kellykochanski/rescal-snow.git
cd rescal-snow
```
You may also download the repository manually from [https://github.com/kellykochanski/rescal-snow](https://github.com/kellykochanski/rescal-snow).

**All further blocks of bash instructions start from this directory.**

### Installation

These instructions will get rescal-snow running on most linux environments; for additional installation options, tips on avoiding/installing missing dependencies, and MacOS installation instructions see [docs/how_to_install.md](how_to_install.md).

In a terminal, navigate into the main rescal-snow directory (shown above). Run:
```bash
  ./configure
  make
```

## 3. Use and features <a name="examples"></a>

To start running, configuring and analyzing Rescal-snow simulations, go through the tutorial: [docs/rescal-snow-tutorial.md](docs/rescal-snow-tutorial.md).

The tutorial fully describes the examples below. The [docs](docs) folder also contains descriptions of additional configuration and analysis options; check these if you're looking for functionality not found in the tutorial.

### Sand and snow dunes <a name="test-cone"></a>
The following examples are described fully in [docs/rescal-snow-tutorial.md](docs/rescal-snow-tutorial.md).

The default configuration for Rescal-snow simulates snow (or sand) dune formation. This simulates processes including air flow; grain entrainment, saltation, suspension and deposition; and granular avalanches.

In these conditions, a pile of sand/snow (left) evolves into a dart-shaped barchan dune (middle), then dwindles as grains blow away without being resupplied (right).

|Initial condition, t = 0t0 	|  t = 30t0   | t = 90t0  |
|------------------------|--------------------|------------------|
| ![](docs/example_images/snow_cone/00.png) | ![](docs/example_images/snow_cone/03.png) | ![](docs/example_images/snow_cone/09.png) |

Each of the three images above shows a shaded top-down view of a dune (top left), cross-sections through the dune, along the dashed lines (middle left, top right), and a cross-section showing the pressure intensity in the fluid (bottom left).

### Sintering snow <a name="test-sinter"></a>

Snow cohesion increases over time: this is called sintering.
Rescal-snow is able to simulate the transition of loose (beige) grains into sintered (light purple) grains within waves and dunes.

![](docs/example_images/sintering/sintering.gif)

### Dune growth by snowfall <a name="test-snowfall"></a>

Rescal-snow simulates snow by adding loose grains to the top of the simulation. The gif at the top of this page shows a height-map of a field of dunes and waves growing during snowfall.

### Parallel instances and parameter space exploration <a name="test-parallel"></a>
We believe that building robust, trustworthy models is much simpler when it's easy to make many model runs.

Rescal-snow contains a series of tools for running many simulation instances in parallel, and managing the associated flows of input and output data.

The following phase diagram shows images produced by ten parallel runs simulating different rates of snowfall (lambda\_I) and wind speed (Tau\_min). 

![snowfall-wind phase diagram](docs/example_images/phase_space_exploration/phase_diagram1.png)

## 4. References and further reading <a name="references"></a>

The [docs](docs) folder contains additional information on 
[alternate installations](docs/how_to_install.md), 
[performance and parallelization issues](docs/performance_and_parallelization.md), 
[model inputs and configuration](docs/rescal-snow-inputs.md),
[model calibration and validation](docs/calibration_and_validation.md),
and [development history](docs/NEWS.md)

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

## 5. Community guidelines <a name="community"></a>

We are excited to hear about your scientific work with ReSCAL.
This section lists the best way to bring your project - and its successes, challenges, bugs, and developments - to our attention.

We encourage you to interact with the project through github (see below). This will allow easy integration of your changes and prevent rescal-snow from fragmenting excessively.

### Citation

Do you want to incentivize developers to build and maintain the software you need?
Cite us!

We are working on getting a doi for this repository through the Journal of Open Source Software (JOSS). 

This software inherits many features from the Real-Space Cellular Automaton Laboratory, ReSCAL. Please credit those developers by citing:
 - ['A real-space cellular automaton laboratory', O Rozier and C Narteau, Earth Surface Processes and Landforms 39(1) 98-109, 2013, doi=10.1002/esp.3479](https://onlinelibrary.wiley.com/doi/abs/10.1002/esp.3479)

### Support

If you have challenges or questions, look at the material under 'further information' or reach out to us.

Primary contact:
Kelly Kochanski
kelly.kochanski@gmail.com
www.github.com/kellykochanski

Issues may be reported using github's issue tracking function on this repository, [www.github.com/kellykochanski/rescal-snow](www.github.com/kellykochanski/rescal-snow).

### Contributing

We have built a model of the basic growth and function of snow dunes, but expect that many users may want more detailed features.

If you wish to contribute a new feature, we recommend you fork our repository, commit your changes to a new feature or development branch, and create a pull request. An example contribution workflow, with git instructions, is outlined by the [LAVA software community project contribution
guide](https://docs.lavasoftware.org/lava/contribution.html)

Rescal-snow is distributed under the GNU GPL 3.0 license; all contributions must be made under this license or a later version.


## 6. Contributors <a name="authors"></a>
See [AUTHORS.md](AUTHORS.md).

## 7. License
GNU GPL 3.0 or any later version. See [docs/LICENSE.md](docs/LICENSE.md). SPDX-License-Identifier: GPL-3.0-or-later

Release: LLNL-CODE-785837. See [NOTICE](NOTICE) for details.
