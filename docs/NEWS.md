# News and development history

## rescal-snow 1.0
#### Models, templates and functionality

Added SNO model

Added snowfall ability to SNO model, with accompanying Csp_template SNOWFALL

Added sintering ability to SNOW model

Removed non-dune-related models including RIV, CMB, etc

#### Tools for set-up, post-processing, and parameter space exploration

Added 'scripts/utilities' folder containing python tools for setting up large numbers of runs for parameter space explorations

Added cross-correlation analysis tool

Added fft analysis and visualization tool

#### Organization

Removed many regions of commented-out code

Model output now funnelled through 'write_output' files for standardization

Model output may now be organized in output directories

#### Performance and compatibility

Converted all instances of 'int' to 'int32_t' for cross-platform compatibility

Added additional AVX optimizations to the lattice gas

Removed GTK dependency; graphics now handled entirely through LIBPNG

Modified configure and make pathways to account for new dependencies

## ReSCAL Version 1.6

AVX optimization of lattice-gas propagation.

Optimization of flow interpolation.

New command-line option (-pat) for the concurrent execution of stochastic engine 
and lattice gas, with non-deterministic algorithm.

Command-line help improved with new options -h and -hm.

More CSP templates and parameters for DUN model.

Various fixes.

## ReSCAL Version 1.5
Relief light-shading with various colors, for colored sand and topographic 
obstacles (DUM cells) in DUN and AVA models, and for fluid cells in RIV model.

More simple syntax for graphical options in command-line.

The graphical display readapts automatically as the window gets larger.

Better interactivy with the positioning of cross sections.

New TRACE_FLUX mode in DUN model for the calculation of the average sand flux
over time, compatible with rotations.

New option for stability analysis in DUN model.

New parameters and transtions for the mobility of sand in DUN model.

Introduction of vegetated cells in the DUN model, with specific transitions
and parameters.

Time metadata added in CSP format.

More CSP templates, and arguments added in previous templates.

Many tweaks and fixes.

## ReSCAL Version 1.4
New FULL rotation mode for rotating the whole cellular space with support of
periodic boundary conditions and lattice gas coupling.

More CSP templates and samples.

RIV model improved.

More graphical options and parameters.

Scripts added for video encoding.

Various fixes and optimizations.

## ReSCAL Version 1.3.2
OpenMP support in lgca module.

New trace mode TRACE_PAR_COL for the tracking of colored particles.

## ReSCAL Version 1.3.1
New parameter for boundary conditions.

Library gd not needed anymore.

## ReSCAL Version 1.3
The time scale is calculated in DUN model, and t0 time unit can be 

used for all delays on the command-line.

Cellular space is automatically reoriented on images, after rotations.

New parameters with CELL_COLOR option.

Support for CSP templates in genesis.

New LIFE model implemented.

## ReSCAL Version 1.2
CSP format with metadata has been implemented.

New parameters available for avalanches and lattice gaz.

Interactive display and new disposition in drawing area.

## ReSCAL Version 1.1
Integration of glade interface with toolbar and statusbar.

Compatibility with gtk+ 1.x has been abandoned.

## ReSCAL Version 1.0
This is the first release of ReSCAL, under the GNU General Public 
License.  It is far from complete. But it provides a number of 
predefined models with some parameters and mechanisms, included 
lattice gaz coupling, and real-time display with several graphical 
options.

