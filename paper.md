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
    affiliation: "1, 2, 3" # (Multiple affiliations must be quoted)
  - name: Eric Green
    affiliation: 4
  - name: Carlos Downie
    affiliation: 4
  - name: Adam Rubin
    affiliation: 5
  - name: Barry Rountree
    orcid: 0000-0002-0087-4301
    affiliation: 4
affiliations:
 - name: Department of Geological Sciences, University of Colorado Boulder, Boulder, CO, USA
   index: 1
 - name: Institute for Arctic and Alpine Research, University of Colorado Boulder, Boulder, CO, USA
   index: 2
 - name: Cooperative Institute for Research in Environmental Sciences, University of Colorado Boulder, Boulder, CO, USA
 - index: 3
 - name: Computing Division, Lawrence Livermore National Laboratory, Livermore, CA, USA
 - index: 4
 - name: Environmental Engineering Program, University of Colorado Boulder, Boulder, CO, USA
 - index: 5
date: 21 June 2019
bibliography: paper.bib
---

# Summary

Snow cover and accumulation are vital to many ecosystems.
Each winter, snow drifts from the sky onto the polar ice sheets and sea ice.
It adds mass, brightness, and insulation that cool the Earth in the face of climate change.
In southern regions, snow whirls over prairies and tundra. Where it lands, it covers plants and animals against the winter cold, and melts into their fresh water in the spring.
Most snow that falls, however, does not stay on the ground: it blows away.

Previous snow models do not cover <100m length scales. This means they miss important small scale heterogeneities, which are known to have important climate effects.
They also entirely miss snow dunes and bedforms.
Several authors have shown that these features are shaped by processes such as time-dep cohesion that have not been included in previous geomorph models

We address all of these difficulties:
rescal-snow covers 0.1-100m length scales
it is designed to capture self-org w cellular automata
we have added time-dep cohesion, in a framework easy for user to alter

We also go for good practices, hope useful for snow sci and E surf proc:
open source code etc
building on prev contributions
encourage future contributions
easy parallel runs, scaling analysis, aim to make UQ and param space exploration standard feature of snow sci and geomorph analysis

---

Rescal-snow is designed to give snow scientists new tools for studying the redistribution of snow by wind. In particular, it is designed to capture
snow self-organization and redistribution on 0.1-10m length scales.
Previous snow models of wind-blown snow have focused on
0.1-10 km scales [@Liston:2007; Lehning2002], continental scales [@Gallee2012; @Lenaerts2012] or global scales [@Blanchard2015].
These models have been useful for calculating the large-scale effect of snow on ice sheet mass balances, but still have significant discrepancies when compared to real Antarctic snow measurements (e.g. @Gallee2012)
We know the small length scales are important because...
So we address this with a model focused on 0.1-10m length scales.

Also these models dont get snow dunes...
Which are important because...
The closest models to get those are Liston2018, new Cryosphere paper... 
These are applications of models with a smaller capability, they show the need for a new model that scientists can apply
To address this, we built a model using a cellular automaton, well known for being good for self-organization...
We also add features (snowfall, sintering) that several field authors studying snow dunes have identified as important...

Our examples walk through 'how does snow stick to ground'.

Much of the surface of the Earth is covered by snow. In winter, this cover reaches out from the poles to cover ice sheets, sea ice, and tundra; it collects in 

Statement of need

rescal-snow is...

This contribution describes new work...

Ongoing problems from rescal-snow

# Mathematics

Single dollars ($) are required for inline mathematics e.g. $f(x) = e^{\pi/x}$

Double dollars make self-standing equations:

$$\Theta(x) = \left\{\begin{array}{l}
0\textrm{ if } x < 0\cr
1\textrm{ else}
\end{array}\right.$$


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

We thank Clement Narteau and Oliver Rozier for their advice and support.

This work was supported by a Department of Energy Computational Science Graduate Fellowship (DE-FG02-97ER25308).

*Add Barry support*

*Add computational resources*

# References

---
references:
- id: Clark2012
  author: Clark, M. P., J. Hendrikx, A. G. Slater, D. Kavetski, B. Anderson, N. J. Cullen, T. Kerr, E. O  Hreinsson, and R. A. Woods 
  issued:
  year: 2011 
  title: Representing spatial variability of snow water equivalent in hydrologic and land-surface models: A review
  container-title: Water Resources Research
  volume: 47
  issue: W07539
  DOI: 10.1029/2011WR010745

- id:Lenaerts2012
  author: Lenaerts, J. T. M., M. R. van den Broeke, W. J. van de Berg, E. van Meijgaard, P. Kuipers Munneke
  title: A new, high-resolution surface mass balance map of Antarctica (1979-2010) based on regional atmospheric climate modeling
  container-title: Geophysical Research Letters
  volume: 39
  issue: 4
  year: 2012
  DOI: 10.1029/2011GL050713

- id: Blanchard2015
  author: Blanchard-Wrigglesworth, E., S. L. Farrell, T. Newman, and C. M. Bitz
  title: Snow cover on Arctic sea ice in observations and an Earth System Model
  container-title: Geophysical Research Letters
  volume: 42
  page: 10342-10348
  DOI: 10.1002/2015GL066049
  year: 2015

- id: Lehning2002
  author: Lehning, M., P. Bartelt, B. Brown, C. Fierz
  title: A physical SNOWPACK model for the Swiss avalanche warning: Part III: meteorological forcing, thin layer formation and evaluation
  container-title: Cold Regions Science and Technology
  volume: 35
  issue: 3
  page: 169-184
  DOI: 10.1016/S0165-232X(02)00072-1
---

Filhol, S., and M. Sturm (2015), Snow bedforms: A review, new data, and a formation model, *Journal of Geophysical Research: Earth Surface*, *120*, 9, doi:10.1002/2015JF003529.

Gall\'ee, H., A. Trouvilliez, C. Agosta, C. Genthon, V. Favier, F. Naiim-Bouvet (2012), Transport of snow by the wind: A comparison between observations in Ad\'elie Land, Antarctica, and simulations made with the regional climate model MAR, *Boundary-Layer Meteorology*,
*146*, 1, 133-147.

Kobayashi, S. (1980), Studies on interaction between wind and dry snow surface, *Contributions from the Institute of Low Temperature Science*, *A29*, 1-64.

Kochanski, K., R. S. Anderson, and G. E. Tucker (2019), The evolution of snow bedforms in the Colorado Front Range and the processes that shape them, *The Cryosphere*, *13*, 1267-1281, doi:10.5194/tc-13-1267-2019.

Kochanski, K., R. S. Anderson, and G. E. Tucker (2018), Statistical classification of self-organized snow surfaces, *Geophysical Research Letters*, *45*, 13, 6532-6541, doi:10.1029/2018GL077616.

Liston, G. E., R. B. Haehnel, M. Sturm, C. A. Hiemstra, S. Berezovskaya, R. D. Tabler (2007), Simulating complex snow distributions in windy environments using SnowTran-3D, *Journal of Glaciology*, *53*, 181, 241-256, doi:10.3189/172756507782202865.
