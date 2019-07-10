# Performance and parallelization of rescal-snow

This is a quick description of how rescal run-time varies with domain size and parallelization options.

Author: Kelly Kochanski, 18 June 2019

## Major observations and recommendations

### Rescal-snow performance is largely limited by the cost of doublet transitions between grains.

The simulation speed decreases as the number of available grain transitions increases. This can occur when:
 - many transition types are enabled
 - the horizontal domain (set by L and D) is large
 - the surface is completely covered by grains
 
 The simulation is not highly sensitive to changes in the lattice gas. 
 Changes to the fluid depth (set by height H) has negligible effect on run time.
 Parallelizing the lattice gas component is also largely ineffective.
 
### Rescal-snow is highly resistant to parallelization

During the course of this development, we attempted many parallelizations for rescal-snow.
The details are beyond the scope of this document, but our attempts included
 - putting the lattice gas and cellular automaton on separate processors (developed in ReSCAL v1.6)
 - exporting random number generation to GPUs
 - slicing the domain into separate pieces controlled by separate processors
 None of these efforts led to particularly efficient parallel computation.
 
 We therefore emphasise applications, such as parameter space explorations, that leverage many processors to perform many serial runs, rather than applications that require parallelizing rescal-snow.

*Why is parallelizing rescal-snow hard?*

In algorithmic terms, rescal-snow looks a lot like a [discrete event simulation](https://en.wikipedia.org/wiki/Discrete-event_simulation). These are a notoriously hard class of simulations to parallelize.

In rescal-snow terms, the simulation is very hard to parallelize because it has no 'maximum speed'. Some events, like avalanches, cover large distances nearly instantaneously. It is not easy to identify parts of the simulation that are independent over a given time scale.

Parallelizing cellular automata is an interesting HPC problem. We're still working on it.

### Recommendations for increasing performance

In order of ease and likely effectiveness:

 - Turn off unnecessary transitions
      - Set the transition rate to 0 in the parameter file, or
      - Disable the transition in [src/models.c](src/models.c)
 - Decrease the horizontal domain size, L and D
      - For some applications, the effect of a long domain can be approximated with periodic boundaries
 - Run rescal-snow on a dedicated processor of the highest available speed
 - Decrease output frequency using the `-dcsp` and `-dpng` options in the run script
 - Try using AVX optimizations (see [docs/how_to_install.md](how_to_install.md) )
 - Try compiling with intel compiler vs gcc
 
 We suspect that the largest performance gain would be obtained by replacing some stochastic cellular automaton transitions with deterministic transitions.
 This would reduce the frequency with which random numbers must be generated and random doublets must be selected.
 Such changes could be implemented in [src/models.c](src/models.c), though they would change the structure of rescal-snow,
 and might lead to systematic biases in its behavior.
 
 ## Scaling experiments
 
 Below are select scaling experiments. These are run with the configuration from the [snowfall example](../scripts/snowfall.run).
 'Wall time' is the time taken to reach t0=400.

 
 ### Effects of domain size on run time
 
 The run-time of rescal-snow increases linearly with both horizontal dimensions, D and L:
 
 | L    | walltime (MM:SS) |
|------|----------|
| 150  | 06:29    |
| 300  | 12:30    |
| 600  | 22:05    |
| 1200 | 48:43    |

| D   | wall time  |
|-----|-----------------|
| 50  | 06:37           |
| 100 | 12:30           |
| 200 | 22:01           |
| 400 | 47:20           |
 
 The run-time is not very sensitive to changes in the vertical dimension H:
 
 | H   | walltime |
|-----|----------|
| 50  | 12:30    |
| 100 | 14:10    |
| 200 | 13:49    |

In each table, the dimensions not explicitly varied are held constant at H=50, D=100, L=300.

### Effects of native parallelization

ReSCAL v1.6 implemented a parallelization scheme in which the cellular automaton and lattice gas run on separate processors.
Unfortunately, this parallelization has a very limited effect on run-time in this configuration,
presumably because the cellular automaton takes much more processing power than the lattice gas in this setup.
This parallelization might be more effective on a mostly-empty simulation with few grains and few grain transitions.

For these experiments, H=100, L=300, D=100.

| Cores | walltime |
|-------|----------|
| 1     | 14:10    |
| 2     | 12:38    |
| 2 (again)  | 12:47    |

This shows that adding a second core creates only a 13% speedup, for a parallel efficiency of 54%.
 
