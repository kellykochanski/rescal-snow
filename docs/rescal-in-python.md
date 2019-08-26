
# Running Rescal-snow in python

Run Rescal-snow within python to process output data concurrently

Author: Gian-Carlo DeFazio, defaziogiancarlo@gmail.com

## Running Rescal-snow in python

The parallel run options described in `README.md` may generate a large number of files.
The .csp files, which are snapshots of the entire 3D cell space, can take up large amounts
of memory. Also, the ALTI files, which are height maps of the cell space, may take some
time for ReSCAL to create.

To speed up the simulation and reduce the number of files produced, ReSCAL can be run
as the child process of a python process. Instead of writing .csp and ATLI files,
the contents of the .csp file are send to the parent python process, where they can be
processed, aggregated, and compressed before being written to the output directory.

The python code that calls ReSCAL expects the environment variable RESCAL_SNOW_ROOT to be defined.
This code is in the `datarun.py` and a DataRun instance is created for each data run.
You can check if it is already defined:
```bash
echo $RESCAL_SNOW_ROOT
```

You should see the path to the top directory of you Rescal-snow installation.
If RESCAL_SNOW_ROOT is not defined:
```bash
export RESCAL_SNOW_ROOT='<path to your rescal install>'
```

To set RESCAL_SNOW_ROOT permanently (not recommended if you are doing development with multiple Rescal-snow instances), put
```bash
export RESCAL_SNOW_ROOT='<path to your rescal install>'
```
in your `~/.bashrc` file.

From now on RESCAL_SNOW_ROOT refers to the top directory of Rescal-snow (or $RESCAL_SNOW_ROOT for bash commands)

You will also need a directory to store the data files.
The default is `data_runs` which will be in RESCAL_SNOW_ROOT.
A different location can be specified when creating a DataRun instance.
For each run a subdirectory will be created. These subdirectories will be similar to
those created in the <b>Setting up parallel runs</b> section of `README.md`.
However, some files will be missing, and in the output directory, there will be aggregated
and processed outputs, such as `height_maps.npz` and `ffts.npz` instead of `.csp` and `ALTI`
files. There should also be a `meta_data` file which holds all the parameters used for the data run.

There is an example script `example_pyrescal` which can be run individually,
but is meant to be run using `pyrescal.sbatch`.
To run these files:

First go to the utilities directory
```bash
cd scripts/utilities
```

run an individual example
```bash
./example_pyrescal.py
```

or run several in parallel using sbatch

```bash
sbatch pyrescal.sbatch
```

if your system has a debug partition you may get scheduled sooner using

```bash
sbatch -p pdebug pyrescal.sbatch
```

If you used the sbatch option, a log.txt file will appear.
Check the log.txt file for errors during the run.

If you received an error about the `data_runs` directory, make the directory:
```bash
cd $RESCAL_SNOW_ROOT
mkdir data_runs
```

Now go back to the `utilities` directory and try to run it again.

Once you've successfully completed either a single or parallel run,
check the `data_runs` directory.<br>
If you did a single run you should see:
```bash
ls $RESCAL_SNOW_ROOT/data_runs
>> exp0
```

If you did the sbatch run you should see
```bash
ls $RESCAL_SNOW_ROOT/data_runs
>> exp0 exp1 exp2 exp3 exp4 exp5 exp6 exp7
```

Look at one of the directories
```bash
ls $RESCAL_SNOW_ROOT/data_runs/exp0
>> DUN.csp  meta_data  out  run.par
ls $RESCAL_SNOW_ROOT/data_runs/exp0/out
>> CELL.log  CGV_COEF.log	DENSITE.log  TRANSITIONS.log  VEL.log  ffts.npz  height_maps.npz
```

The `.npz` files contain the processed outputs.
You can read them in using the numpy module in python
```python
cd $RESCAL_SNOW_ROOT/data_runs/exp0/out
python3
>>> import numpy as np
>>> h = np.load('height_maps.npz')['height_maps']
>>> f = np.load('ffts.npz')['ffts']
>>> h.shape
(5, 200, 400)
>>> f.shape
(5, 16, 33)
```
In this case, there are 5 timesteps (0_t0, 100_t0, 200_t0, 300_t0, 400_t0) and the dimensions of the 3D space
were (200, 80, 400) so the heightmaps are (200,400). The ffts are smaller because by default they only show the portion that has been found to have the 
dominant frequencies.

The default processing will make height map and fft files.

The arrays may be viewed with standard processing; see the `visualization` tutorial, or use the tool of your choice (e.g. matplotlib.pyplot.imshow()).

