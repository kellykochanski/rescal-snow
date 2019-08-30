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
    affiliation: "1,2"
  - name: Gian-Carlo Defazio
    affiliation: "2,3"
  - name: Eric Green
    affiliation: 2
  - name: Richard Barnes
    affiliation: "4"
    orcid: 0000-0002-0204-6040
  - name: Carlos Downie
    affiliation: 2
  - name: Adam Rubin
    affiliation: "5"
  - name: Barry Rountree
    orcid: 0000-0002-0087-4301
    affiliation: 2
affiliations:
 - name: Department of Geological Sciences, University of Colorado at Boulder, Boulder, CO, USA
   index: 1
 - name: Computing Division, Lawrence Livermore National Laboratory, Livermore, CA, USA
 - index: 2
 - name: Department of Computer Science, Western Washington University, Bellingham, WA, USA
   index: 3
 - name: Department of Electrical Engineering and Computer Sciences, University of California Berkeley, Berkeley, CA, USA
   index: 4
 - name: Department of Civil Engineering, University of Colorado at Boulder, Boulder, CO, USA
   index: 5
date: 21 July 2019
bibliography: paper/paper.bib
---

# Summary

One of the most intriguing effects of wind on snow is the creation of mysterious and captivating shapes known as snow bedforms. These features ornament Antarctica, Arctic sea ice, tundra, and many mountain tops [@filhol:2015; @kochanski:2018; @kobayashi:1980].
They change the reflectivity and average thermal conductivity of snow, and may well change the patterns of snow accumulation and transport.
They also present a range of interesting self-organized shapes that might interest geomorphologists and Earth surface process scientists.
Despite these effects, however, they are poorly understood and not yet included in major snow models.

Several recent field studies have identified new computational components needed for a good understanding of snow bedforms.
First, most snow bedforms (e.g ripples, barchan dunes, snow-steps and sastrugi) are between 0.1 and 2 m in length, with select bedforms (e.g. snow-waves, some whaleback dunes) extending from 5 to 30 m.
These length scales are based on physical phenomena such as the hop length of blowing snow grains [@kobayashi:1980] and the length scales of turbulent structures in the wind [@kobayashi:1980, @kochanski:2019].
Existing models of wind-blown snow, however, are designed to model snow redistribution over mountainous [@liston:2007; @lehning:2002; @marsh:2018] or continental [@gallee:2013] scales
(with grid sizes of 0.1-10 km)
and are unable to resolve the processes that lead to the formation of snow bedforms.
Second, snow bedforms are the result of a balance between sand-like granular motion and cohesive resistance [@kochanski:2019; @filhol:2015].
Existing sand dune models (e.g. @narteau:2014, @lammel:2012 ) successfully model a wide range of granular processes, including sand avalanches and wind-driven saltation, suspension, and creep.
These models, however, have not had reason to include the kinds of time-dependent cohesion that occur in snow.
Similarly, the simulations that have advanced our understanding of cohesion and sintering in snow (e.g. @herwignen:2013 , @colbeck:1983 ) have focused on immobile snow, without considering the effects of interrupting metamorphosis with motion.
Third, snow bedforms tend to move during and immediately after snowfall events [@kochanski:2018], which are also not included in existing dune models.
These features collectively indicate a need for a new computational model that maps snow processes, such as cohesion and snowfall, into a granular physics framework.

Rescal-snow is designed to enable the quantitative study of snow bedforms.
It simulates 10-100 m domains at 0.05-0.20 m resolution, allowing it to capture all but the smallest snow bedforms.
The simulation is adapted from a sand dune model, ReSCAL [@narteau:2014], and inherits ReSCAL's granular transport capacities.
We have added features to simulate processes unique to snow, including snowfall and time-dependent cohesion.
Fluid processes are modelled with a lattice gas cellular automaton, a method chosen to provide a good approximation to the Navier-Stokes equation at reasonable computational cost.
Finally, the backbone of the simulation is a cellular automaton, an algorithm known to be particularly good for modelling self-organization.

We have been able to use Rescal-snow to simulate the formation and movement of snow dunes and snow-waves under a range of wind, snowfall, and sintering conditions.
We illustrate these results through the example simulations in our `README` file, 
and in the [associated tutorial](docs/rescal-snow-tutorial.md), which also leads readers through a previously-unmodelled scientific question: how do bedforms affect the accumulation of snow?

The present limitations of Rescal-snow mostly concern length scales, time scales, and the representation of turbulence.
The output of Rescal-snow is disretized in length, which limits its ability to resolve tiny (<0.05 m) features, and by cell type, which limits its ability to represent the full spectrum of gradual changes in the nature of snow, such as those that occur
during sintering.
The natural length and time scales of Rescal-snow are set by the configuration of the cellular automata.
Although these can be related to real length and time scales, it requires careful calibration: this is described fully in @narteau:2009.
Moreover, we have not yet been able to model particularly sharp-edged snow features, such as snow-steps or sastrugi.
We attribute this to Rescal-snow's relatively simple fluid dynamics algorithm and inability to resolve the flow detachments [@kochanski:2019] that likely form at the corners of snow-steps and sastrugi.

We aim to demonstrate good practices that will encourage robust, reproducible science by releasing Rescal-snow through Journal of Open Source Science.
For example, we expect that our work will be many users' first introduction to some subset of bash, C, Python or high-performance computing, and we endeavored to make this a positive learning experience. 
Our examples include scientifically interesting applications of all of these skills, and our tutorial is punctuated with references to relevant tutorials.
As a second example, we aim to make it easier for geoscientists to run robust numerical experiments which explore many relevant variables and allow for trustworthy uncertainty quantification.
We have therefore set up structures for configuring, running, and analyzing parallel simulation instances.
We hope that the above practices will enable users to run high-quality numerical experiments with Rescal-snow.

# Acknowledgements

This work was performed under the auspices of the U.S. Department of Energy by Lawrence Livermore National Laboratory under Contract DE-AC52-07NA27344. This paper is released under LLNL-JRNL-786878-DRAFT.

We thank Clement Narteau and Oliver Rozier (IPGP) for advice and support in development beginning with ReSCAL 1.6,
Robert Anderson and Gregory Tucker (CU) for advice on the scientific direction of this software,
and Tapasya Patki, Divya Mohan, Jeff Booher-Kaeding and Aaron Robeson (LLNL) for minor contributions to Rescal-snow.

The project was supported by a Department of Energy Computational Science Graduate Fellowship (DE-FG02-97ER25308), by support from the Data Science Summer Institute at Lawrence Livermore National Laboratory, and by a UROP award from the University of Colorado.
