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
    affiliation: "1,2" # (Multiple affiliations must be quoted)
  - name: Eric Green
    affiliation: 2
  - name: Gian-Carlo Defazio
    affiliation: "2,3"
  - name: Carlos Downie
    affiliation: 2
  - name: Adam Rubin
    affiliation: 1
  - name: Barry Rountree
    orcid: 0000-0002-0087-4301
    affiliation: 2
affiliations:
 - name: University of Colorado at Boulder, Boulder, CO, USA
   index: 1
 - name: Computing Division, Lawrence Livermore National Laboratory, Livermore, CA, USA
 - index: 2
 - name: Western Washington University, Bellingham, WA, USA
   index: 3
date: 21 July 2019
bibliography: docs/paper.bib
---

# Summary

*Abstract*

*Paper*

One of the most intriguing effects of wind on snow is the creation of mysterious and captivating shapes known as snow bedforms. These features ornament Antarctica, Arctic sea ice, tundra, and many mountain tops [@filhol:2015; @kochanski:2018; @kobayashi:1980].
They change the reflectivity and average thermal conductivity of the snow, and may well change the patterns of snow accumulation and transport.
They also present a range of interesting self-organized shapes that might interest geomorphologists and Earth surface process scientists.
Despite these effects, however, they are poorly understood and not yet included in major snow models.

Several recent field studies have pointed to the missing theoretical components needed for a good understanding of snow bedforms.
First, most snow bedforms (e.g ripples, barchan dunes, snow-steps and sastrugi) are between 0.1 and 2 m in length, with select bedforms (e.g. snow-waves, some whaleback dunes) extending from 5 to 30 m.
These length scales are not arbitrary. They are tied to physical features, such as the hop length of blowing snow grains [@kobayashi:1980] and the length scales of turbulent structures in the wind [@kobayashi:1980, @kochanski:2019].
Existing models of wind-blown snow, however, are designed to model snow redistribution over mountainous [@liston:2007; @lehning:2002; @marsh:2018] or continental [@gallee:2013] scales,
and are unable to resolve the processes that lead to the formation of snow bedforms.
Second, snow bedforms are the result of a balance between sand-like granular motion and cohesive resistance [@kochanski:2019; @filhol:2015].
Existing sand dune models (e.g. @narteau:2014, @lammel:2012 ) successfully model a wide range of granular processes, including sand avalanches and wind-driven saltation, suspension, and creep.
These models, however, have not had reason to include the kinds of time-dependent cohesion that occur in snow.
Similarly, the simulations that have advanced our understanding of cohesion and sintering in snow (e.g. @herwignen:2013 , @colbeck:1983 ) have focused on immobile snow, without considering the effects of interrupting metamorphosis with motion.
Finally, none of the models above has a mechanism for capturing the effects of self-organization, or for representing any emergent surface features in snow. 
 
Rescal-snow is designed to enable the quantitative study of snow bedforms.
It is designed to simulate 10-100 m domains at 0.05-0.20 m resolution, allowing it to capture all but the smallest snow bedforms.
The simulation is adapted from a sand dune model, ReSCAL [@narteau:2014], and inherits ReSCAL's granular transport capacities.
We have added features to simulate processes unique to snow, including snowfall and time-dependent cohesion.
Fluid processes are modelled with a lattice gas cellular automaton, a method chosen to provide a good approximation to the Navier-Stokes equation at reasonable computational cost.
Finally, backbone of the simulation is a cellular automaton, an algorithm known to be particularly good for modelling self-organization.

We have been able to use Rescal-snow to simulate the formation and movement of snow dunes and snow-waves under a range of wind, snowfall, and sintering conditions.
We illustrate these results through the example simulations in our `README` file, which also leads readers through an answer to a previously-unmodelled scientific question: how do bedforms affect the accumulation of snow?

The present limitations of Rescal-snow mostly concern length scales, time scales, and the representation of turbulence.
Thus far, we have not been able to model particularly sharp-edged snow features, such as snow-steps or sastrugi.
We attribute this to Rescal-snow's relatively simple fluid dynamics algorithm and inability to resolve the flow detachments [@kochanski:2019] that likely form at the corners of snow-steps and sastrugi.
As a cellular automaton, the output of Rescal-snow is disretized in time, which limits its ability to resolve tiny (<0.05 m) features, and by cell type, which limits its ability to represent the full spectrum of gradual changes in the nature of snow, such as those that occur
during sintering.
Both of these limits can be adjusteded; here we present a configuration with a practical balance between precision and computational expense.
Finally, the natural length and time scales of Rescal-snow are set by the configuration of the cellular automata.
Although these can be related to real length and time scales, it requires careful calibration: this is described fully in @narteau:2009.

By releasing Rescal-snow through the Journal of Open Source Science, we aim to demonstrate good practices that will encourage robust, reproducible science.
For example, we expect that our work may be the first introduction many of our readers have to some subset of bash, C, python, and high-performance scientific computing.
We have endeavored to make this a positive learning experience. Our examples walk through example applications of all of these skills, and our README is punctuated with references to relevant tutorials.
We also want to encourage our users to design and run rigorous numerical experiments.
To aid this, we have set up structures that make it straightforward to run parallel instances of Rescal-snow (and we demonstrate a parallel workflow in the examples).
We hope that this will encourage users to get more than the minimum number of data points needed for their phase-space experiments or inverse model runs, 
and will generally enable rigorous science and routine quantification of uncertainty. 

---


# Citations

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

For a quick reference, the following citation commands can be used:
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"

# Acknowledgements

We thank Clement Narteau and Oliver Rozier (IPGP) for advice and support in development beginning with ReSCAL 1.6,
Robert Anderson and Gregory Tucker (CU) for advice on the scientific direction of this software,
and Tapasya Patki, Divya Mohan, Jeff Booher-Kaeding and Aaron Robeson (LLNL) for minor contributions to Rescal-snow.

This work was supported by a Department of Energy Computational Science Graduate Fellowship (DE-FG02-97ER25308), by support from the Data Science Summer Institute at Lawrence Livermore National Laboratory, and by an UROP award from the University of Colorado.

*Add DSSI grant #*

*Add Barry support*

*Add computational resources*

