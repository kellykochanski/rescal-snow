---
title: 'Rescal-snow: Simulating snow dunes with cellular automata'
tags:
  - C
  - snow
  - geomorphology
  - Earth surface
  - dunes
authors:
  - name: Kelly Kochanski
    orcid: 0000-0003-0754-0184
    affiliation: "1,2,3"
  - name: Gian-Carlo Defazio
    affiliation: "1,4"
  - name: Eric Green
    affiliation: 1
  - name: Richard Barnes
    affiliation: "5,6"
    orcid: 0000-0002-0204-6040
  - name: Carlos Downie
    affiliation: 1
  - name: Adam Rubin
    affiliation: "7"
  - name: Barry Rountree
    orcid: 0000-0002-0087-4301
    affiliation: 1
affiliations:
 - name: Computing Division, Lawrence Livermore National Laboratory, Livermore, CA, USA
   index: 1
 - name: Department of Geological Sciences, University of Colorado at Boulder, Boulder, CO, USA
   index: 2
 - name: Institute for Arctic and Alpine Research, University of Colorado at Boulder, Boulder, CO, USA
   index: 3
 - name: Department of Computer Science, Western Washington University, Bellingham, WA, USA
   index: 4
 - name: Energy and Resources Group, University of California Berkeley, Berkeley, CA, USA
   index: 5
 - name: Department of Electrical Engineering and Computer Sciences, University of California Berkeley, Berkeley, CA, USA
   index: 6
 - name: Department of Environmental Engineering, University of Colorado at Boulder, Boulder, CO, USA
   index: 7
date: 21 July 2019
bibliography: paper.bib
---

# Summary

When wind blows over dry snow, it creates shapes known as snow bedforms.
These features, which include snow dunes, waves, snow-steps and sastrugi (figure below), ornament Antarctica, Arctic sea ice, tundra, and mountain ridges [@Filhol2015; @Kochanski2018; @Kobayashi1980].
They change the reflectivity and average thermal conductivity of snow, and may change the patterns of snow accumulation and transport.
Despite these effects, however, they are poorly understood and their effects are not yet included in major snow or climate models.

![](../docs/example_images/field_examples.png)

_Snow bedforms on Niwot Ridge, Colorado. From left to right, small dunes and sastrugi (looking upwind), snow dunes (looking downwind), snow-waves (looking upwind)._

## Requirements of a snow bedform model
Recent field studies have identified three new computational components needed for a good understanding of snow bedforms [@Kochanski2018; @Kochanski2019; @Filhol2015].
First, most snow bedforms (e.g ripples, barchan dunes, snow-steps and sastrugi) are between 0.1 and 2 m in length, with select bedforms (e.g. snow-waves, some whaleback dunes) extending from 5 to 30 m.
These length scales are based on physical phenomena such as the hop length of blowing snow grains [@Kobayashi1980] and the length scales of turbulent structures in the wind [@Kobayashi1980, @Kochanski2019].
Existing models of wind-blown snow, however, are designed to model snow redistribution over mountainous [@Liston2007; @Lehning2002; @Marsh2018] or continental [@Gallee2012] scales
(with grid sizes of 0.1-10 km)
and are unable to resolve the processes that lead to the formation of snow bedforms.

Second, snow bedforms are the result of a balance between sand-like granular motion and cohesive resistance [@Kochanski2019; @Filhol2015].
Existing sand dune models (e.g. @Narteau2014, @Lammel2012) successfully model a wide range of granular processes, including sand avalanches and wind-driven saltation, suspension, and creep.
These models, however, have not had reason to include the kinds of time-dependent cohesion that occur in snow.
Similarly, the simulations that have advanced our understanding of cohesion and sintering in snow (e.g. @Colbeck1983, @Lehning2002a) have focused on immobile snow.

Third, snow bedforms tend to move during and immediately after snowfall events [@Kochanski2018], which are also not included in existing dune models.
Thus, computational studies of snow bedforms will require a model that incorporates snow processes, such as cohesion and snowfall, into a granular physics framework at a relatively high (<0.1 m) resolution.

## Features of Rescal-snow
Rescal-snow is designed to enable the quantitative study of snow bedforms.
It simulates 10-100 m domains at 0.05-0.20 m resolution, allowing it to capture all but the smallest snow bedforms.
The simulation is adapted from a cellular automaton sand dune model, ReSCAL [@Narteau2014], and inherits ReSCAL's granular transport capacities.
We have added features to simulate processes unique to snow, including snowfall and time-dependent cohesion.
Fluid processes are modeled with a lattice gas cellular automaton, a method chosen to provide a good approximation to the Navier-Stokes equation at reasonable computational cost.
Finally, the backbone of the simulation is a cellular automaton, an algorithm known to be particularly good for modeling self-organization.

![](../docs/example_images/rescal-snow_transitions.png)

This simulation will allow snow scientists to translate field studies, which are location-specific, into general terms.
It will also make it easier to investigate the effects of snow bedforms on (1) surface roughness, (2) snow cover fractions, and (3) accumulation rates. This will allow us to describe the effects of meter-scale bedforms in terms of variables that affect regional climates.

### Example simulations
We have been able to use Rescal-snow to simulate the formation and movement of snow dunes and snow-waves under a range of wind, snowfall, and sintering conditions.
We illustrate these results through the example simulations in our `README` file,
and in the [associated tutorial](../docs/rescal-snow-tutorial.md), which also leads readers through a previously-unmodeled scientific question: how do bedforms affect the accumulation of snow?

### Limitations
The natural length and time scales of Rescal-snow are set by the configuration of the cellular automata.
Although these can be related to real length and time scales, it requires careful calibration, as described fully in @Narteau2009.
The length scale of Rescal-snow cells, for reasonable model configurations, is 0.05-0.15 m. This limits the model's ability to resolve centimeter-scale snow features, such as snow-steps [@Kochanski2019], and to successfully model sharp-edged features, such as sastrugi.

As a cellular automata, Rescal-snow cells also have discrete states. Our sintering process, for example, is binary: it includes 'loose' and 'sintered' snow grains, but does not represent intermediate states. (Although it is possible to implement an arbitrarily large number of cell states in the cellular automaton, this incurs a significant performance cost.) Thus, Rescal-snow only approximates the effects of snow metamorphosis and cannot represent those processes with the level of subtlety available in continuous models such as @Lehning2002a.

### Good practices in computational snow science
We aim to demonstrate good practices that will encourage robust, reproducible science by releasing Rescal-snow through Journal of Open Source Science.
Our work is aimed at geomorphologists and snow scientists, and we use this model frequently while working with students.
We expect that our work will be many users' first introduction to some subset of bash, git, C, Python or high-performance computing, and we aim to make this a positive learning experience.
Therefore, we designed our examples around scientific applications of all of these skills, and we punctuate them with references to relevant tutorials on git, bash, etc.

We also believe that good computational science is easier when users are able to make large numbers of model runs.
This allows users order to test the stability of the model, explore a wide range of physical parameters, and gain an accurate understanding of the model uncertainty.
We have set up structures for configuring, running, and analyzing parallel simulation instances to enable users to run high-quality numerical experiments with Rescal-snow.

# Acknowledgments

This work was performed under the auspices of the U.S. Department of Energy by Lawrence Livermore National Laboratory under Contract DE-AC52-07NA27344. This paper is released under LLNL-JRNL-786878-DRAFT.

We thank Clement Narteau and Oliver Rozier (IPGP) for advice and support in development beginning with ReSCAL 1.6,
Robert Anderson and Gregory Tucker (CU) for advice on the scientific direction of this software,
and Tapasya Patki, Divya Mohan, Jeff Booher-Kaeding and Aaron Robeson (LLNL) for minor contributions to Rescal-snow.

The project was supported by a Department of Energy Computational Science Graduate Fellowship (DE-FG02-97ER25308), by support from the Data Science Summer Institute at Lawrence Livermore National Laboratory, and by a UROP award from the University of Colorado.

# References
