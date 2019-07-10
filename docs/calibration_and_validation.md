#Calibrating and validating rescal-snow

How to compare rescal-snow results with real dunes and snow features

Author: KK, 10 June 2019

## Natural length, time, and stress scales in rescal-snow

rescal-snow is a cellular automaton. It has a natural unit length : the size of one cell. This length has no meaning whatsoever in the real world. 
To reproduce real snow features with rescal-snow, you will need to calculate the real lengths, times, and speeds of the features in rescal-snow.

rescal-snow has three non-dimensional scales:

 - A length scale, $l_0$, equal to the length of one cell
 - A time scale, $t_0$, equal to the time taken for one cell to erode in the wind (all other rates and times are relative to this)
 - A stress scale, $\tau _0$, which controls the strength of the coupling between the lattice gas (wind) and cellular automata (grains)

## Calibration approach

[Narteau et al, 2009](dx.doi.org/10.1029/2008JF001127) claim that each of these lengths can be calculated by matching an emergent property of the cellular automaton to an emergent property of real dunes:

 - The length scale can be calibrated by matching the maximum unstable wavelength of the simulation to the real maximum unstable wavelength as calculated from geomorphological theory
 - The stress scale can be calibrated by matching the ratio between wind speed and threshold wind speed from the simulation with the same ratio as estimated in reality
 - The time scale can be calculated by matching the simulated saturated flux of grains to the real saturated flux of grains

This process is described in detail in [Narteau et al, 2009](dx.doi.org/10.1029/2008JF001127).  

## Calibration warnings and challenges

This calibration presents two major model validation difficulties:
 - Changing the model parameters will change the simulated saturated fluxes, and require recalibration of the time scale.
 - The calibration relies on measurements of physical properties, such as the threshold wind speed, that may be near-impossible to measure in the field. This can add significant uncertainty to rescal-snow output.
