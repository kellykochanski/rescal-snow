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

*Snow redistribution by wind is important (emph on med length scales)*

Snow blows in the wind. This movement has many effects on polar climates and ecosystems.
It dulls the reflectivity of the snow that insulates the polar regions; it increases the rate of sublimation of snow in dry air; and it strips the winter insulation from tundra.
One of the most intriguing effects of the wind on snow, however, is the creation of mysterious and captivating shapes known as snow bedforms (Fig 1a-c). These features ornament Antarctica, Arctic sea ice, tundra, and many mountain tops.
They change the reflectivity and average thermal conductivity of the snow, and may well change the patterns of snow accumulation and transport.
They also present a range of interesting self-organized shapes that might interest geomorphologists and Earth surface process scientists.
Despite these effects, however, they are poorly understood and not yet included in major snow models.

Several recent field studies have pointed to the missing theoretical components needed for a good understanding of snow bedforms.
Most known bedforms are between 0.1 and 5 m in length, with a few long-wavelength features (such as snow-waves) reaching wavelengths of 5-30 m.
The features are clearly shaped by their interaction with the wind and may well be sustained by turbulent eddies or large-scale turbulent structures.
A significant fraction of snow bedforms resemble features found in sand, such as dunes and ripples, and experience some granular properties, such as avalanches of dry grains,
 but other bedforms show evidence of time-dependent cohesion, a feature not yet well understood in the Earth surface processes community.
The bedforms also tend to grow during blizzards and snowfall events.

Rescal-snow is designed to enable the quantitative study of snow bedforms.
In contrast to previous snow models, it covers 0.1-100 m length scales. 
The simulation is adapted from a sand dune model, which has been used successfully in geomorphology to model granular processes including aeolian transport and avalanching, plus processes that are unique to snow,
like time-dependent cohesion and snowfall.
The backbone of the simulation is a cellular automaton, an algorithm known to be particularly good for modelling self-organization.
Fluid processes are modelled with a lattice gas cellular automaton, a method chosen to provide a good approximation to the Navier-Stokes equation at reasonable computational cost.

We illustrate each of these features in the example simulations in our `README` file, which is designed to walk readers through a previously-unmodelled scientific question: how does falling snow collect into bedforms, and how does loose blowing snow end up stuck to the
ground?

*Previous models do not cover small length scales or self-organization -> Statement of need*

Previous snow models do not cover <100m length scales. This means they miss important small scale heterogeneities, which are known to have important climate effects.
They also entirely miss snow dunes and bedforms.
Several authors have shown that these features are shaped by processes such as time-dep cohesion that have not been included in previous geomorph models, and that are not easy to modify in existing snow models --- need space for process-based geomorph

Also there are models of snow microstructure but these are too small. However, they do inform some of our physics choices like sensible rates of cohesion increase.

*Good practices - make review easy for JOSS*

We also go for good practices, hope useful for snow sci and E surf proc:
open source code, licensed under GNU GPL 3 or later version (see LICENSE.md),
installation instructions,
new files carefully documented (though some inherited from older work, clearly stated at top of each file),
community guidelines to encourage future contribution, inc. tutorials designed to be accessible to researchers with diverse computational experience.

Also designed to be accessible to scientists with different levels of computational background by including clear prerequisites and pointers towards useful tutorials where those may be necessary.

Also have done our best to set an example of good contribution for the snow science community: continued dev of open-source software, clear acknowledgement of previous work and authors, developing in a direction that encourages good, robust, reproducible computational science by updating performance, giving directions for parallel runs for UQ/rigorous param space exploration, etc.


TODO: Figure with maybe real dune and wave, then rescal dune and wave 
---


# Citations

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

For a quick reference, the following citation commands can be used:
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"

# Figures

Figures can be included like this: ![Example figure.](figure.png)

# Acknowledgements

We thank Clement Narteau and Oliver Rozier (IPGP) for their advice and support in development beginning with ReSCAL 1.6,
Robert Anderson and Gregory Tucker (CU) for their advice on the scientific direction of this software,
and Tapasya Patki, Divya Mohan, Jeff Booher-Kaeding and Aaron Robeson (LLNL) for contributions to the quality and performance of rescal-snow.

This work was supported by a Department of Energy Computational Science Graduate Fellowship (DE-FG02-97ER25308), by support from the Data Science Summer Institute at Lawrence Livermore National Laboratory, and by an UROP award from the University of Colorado.

*Add DSSI grant #*

*Add Barry support*

*Add computational resources*

# References

