# Rescal-snow utilites
This folder contains scripts for configuring, running and analyzing Rescal-snow simulations. These are run through Python3 and bash.
They make the Rescal-snow workflow more consistent and reproducible by making it possible to do most of the work in Python and store both input and output in Python classes.

## [rescal_utilities.py](rescal_utilities.py)
Utilities to set up Rescal-snow runs quickly and easily.
These tools are designed to aid parameter space explorations, sensitivity analyses, and large batches of runs.
The following classes manage most Rescal-snow options and inputs (excepting those in the "real_data" file):

  *Parameters* : Class to hold, update, change, read, or write all the parameters that ReSCAL needs to run

  *RunScript*  : Class to hold, update, change, read, or write a ReSCAL run script with appropriate flags

  *DesignRun*  : Umbrella class to hold, update, change, read or write all parameters, sorting them into Parameters and RunScript

Example usage (python3):
```python
# Create the DesignRun object and (optional but advised) a basic description
my_run     = rescal_utilities.DesignRun()
my_run.set_name("my_run")
my_run.set_header("An example run varying L, D and Lambda_S, and looking for the rescal executable in build")

# Describe values for some subset of Rescal-snow parameterescal.sbatch
# Others will take defaults from rescal_utilities.Parameters._default_parameters() and rescal_utilities.RunScript._default_options()
parameters = {'L':500, 'D':50, 'Lambda_S':0.01, 'rescallocation':'../../build'}
my_run.set_parameters(parameters)
for the output should be created. The default is RESCAL_SNOW_ROOT/data_runs.
# Create a run script named my_run.run and a parameter script named my_run.par
my_run.write()
```
See also the usage in [scripts/utilities/param_space_exploration_example.py](param_space_exploration_example.py)

### Example: [param_space_exploration_example.py](param_space_exploration_example.py)
This is an example file showing how to run a parameter space exploration. In this example, most simulation parameters are set to fixed values, then 30 runs are created to vary:

- five snowfall rates, controlled by parameter Lambda_I
- six wind speeds, controlled by parameter Tau_min,
- and all 30 combinations

The script creates 30 new directories, one for each run. Each directory is seeded with a run script, a parameter file, and rescal and genesis executables. This allows all 30 simulations to run in parallel withou
t interferance.

Usage is described in example 4 of [../../docs/rescal-snow-tutorial.md](../../docs/rescal-snow-tutorial.md).

## [datarun.py](datarun.py)
A python interface to Rescal-snow. This sets Rescal-snow parameters for one or more runs, runs the simulation, and keeps the results in pythonic arrays for ongoing processing.
A DataRun object takes in the parameters and meta-parameters required to run Rescal-snow. The DataRun can receive and processthe output of Rescal-snow while Rescal-snow is running.
To run Rescal-snow using a DataRun instance, the environment variable RESCAL_SNOW_ROOT should be defined and be the path of a Rescal-snow installation. Also, a directory for the output should be created. The default is RESCAL_SNOW_ROOT/data_runs.

```python
# Describe some subset of Rescal-snow parameters and a directory for this run
parameters = {'L':500, 'D':50, 'stop after':'200_t0', 'output interval':'50_t0'}
dirname    = 'myrundirectory'

# Create and run the run
thisrun    = datarun.DataRun(parameters, dirname)
thisrun.setup()
thisrun.run()
```

Example usage is shown in [example_pyrescal.py](example_pyrescal.py)

### Example: [example_pyrescal.py](example_pyrescal.py)
This example tests the functionality of [datarun.py](datarun.py) in a high-performance computing environment. Usage is described in [../../docs/rescal-in-python.md](../../docs/rescal-in-python.md)

## [heightmap.py](heightmap.py)
Utilities to create, read and visualize Rescal-snow height maps, a 2D surface height map which is standard Rescal output. These are managed through a HeightMap class.

Utilities for creating Rescal-snow heightmaps (these may be used as initial conditions for the Rescal-snow simulation; see the `INPUT_ELEVATION CSP_TEMPLATE` in `src/genesis.c`).

- `invader_template`
 - `gaussian_hill`
 - `make_sinusoid`
 - `scale`

Utilities for analyzing Rescal-snow heightmaps (most of these use fourier transforms, as snow/sand self-organization has strong emergent wavelengths):

 - `fft2d_analyze`
 - `fft2d_analyze_map_pic`
 - `fft2d_crop_blur`
 - `fft2d_center_blur`

 Utilities for visualizing height maps:

  - `make_surface`

## [cellspace.py](cellspace.py)
Utilities to read, write, modify and visualize Rescal-snow cell spaces, the full state of the cellular automaton. Managed through `CellSpace` class. Used by [datarun.py](datarun.py).
May be used for visualizing and modifying internal sand/snow structures that are not adequately represented by 2D HeightMaps.
