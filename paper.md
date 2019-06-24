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
    affiliation: "1" # (Multiple affiliations must be quoted)
  - name: Eric Green
    affiliation: 2
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
date: 21 June 2019
bibliography: docs/paper.bib
---

# Summary

*Abstract*

*Snow redistribution by wind is important (emph on med length scales)*

Snow cover and accumulation are vital to many ecosystems.
Each winter, snow drifts from the sky onto the polar ice sheets and sea ice.
It adds mass, brightness, and insulation that cool the Earth in the face of climate change.
In southern regions, snow whirls over prairies and tundra. Where it lands, it covers plants and animals against the winter cold, and melts into their fresh water in the spring.
Most snow that falls, however, does not stay on the ground: it blows away.

*Previous models do not cover small length scales or self-organization -> Statement of need*

Previous snow models do not cover <100m length scales. This means they miss important small scale heterogeneities, which are known to have important climate effects.
They also entirely miss snow dunes and bedforms.
Several authors have shown that these features are shaped by processes such as time-dep cohesion that have not been included in previous geomorph models, and that are not easy to modify in existing snow models --- need space for process-based geomorph

*rescal-snow covers wind redistribution + snow self-orgazation in a process-based framework on <100m length scales*

We address all of these difficulties:
rescal-snow covers 0.1-100m length scales
it is designed to capture self-org w cellular automata
we have added time-dep cohesion, in a framework easy for user to alter

*Ongoing work with rescal-snow*

My stuff
Authors interested in rescal for snow dunes; we'll make their work easier for them

*Good practices - make review easy for JOSS*

We also go for good practices, hope useful for snow sci and E surf proc:
open source code, licensed under GNU GPL 3 or later version (see LICENSE.md),
installation instructions,
new files carefully documented (though some inherited from older work, clearly stated at top of each file),
community guidelines to encourage future contribution, inc. tutorials designed to be accessible to researchers with diverse computational experience.

Examples: walk through an obv. problem this software aims to solve: 'how does snow stick through the ground?', unique new features of this software. 
Also designed to be accessible to scientists with different levels of computational background by including clear prerequisites and pointers towards useful tutorials where those may be necessary.

Also have done our best to set an example of good contribution for the snow science community: continued dev of open-source software, clear acknowledgement of previous work and authors, developing in a direction that encourages good, robust, reproducible computational science by updating performance, giving directions for parallel runs for UQ/rigorous param space exploration, etc.

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

We thank Clement Narteau, Oliver Rozier for their advice and support in development beginning with ReSCAL 1.6,
Robert Anderson and Gregory Tucker for their advice on the scientific direction of this software,
and Tapasya Patki, Divya Mohan, Jeff Booher-Kaeding and Aaron Robeson for contributions to the quality and performance of rescal-snow.

This work was supported by a Department of Energy Computational Science Graduate Fellowship (DE-FG02-97ER25308), by support from the Data Science Summer Institute at Lawrence Livermore National Laboratory, and by an UROP award from the University of Colorado.

*Add DSSI grant #*

*Add Barry support*

*Add computational resources*

# References

