InVivo Toolbox
===========

A selection of MATLAB code from my PhD years for analyzing neural spike trains (and related signals). Main features include oscillation detection and analysis, burst and pause detection, synchrony quantification, and some time series regression.

Functions in the main directory are suitable for general usage. The batch_functions subdirectory contains code for my more specific processing pipelines, but may serve as a useful base for your own pipelines when using code from this toolbox. Namely, these functions maintain all data in a single struct which also includes results from analyses and the parameters used to obtain those results for easier replication. External dependencies not written by me are in the ext_depend subdirectory.

See also the related repository tcwhalen/Whalen2020 to see examples of many of these functions used to generate the figures in Whalen et al. 2020, "Delta oscillations are a robust biomarker of dopamine depletion severity and motor dysfunction in awake mice". That repository includes disconnected copies of code from this repo which are frozen in time to replicate those figures. In contrast, this repository undergoes changes, so code may not be identical between the two.

This code (except external dependencies) was written by Tim C. Whalen in the labs of Aryn Gittis (Carnegie Mellon University) and Jon Rubin (Univeristy of Pittsburgh). Thanks also to Rob Turner (University of Pittsburgh) who provided some dependencies and initial versions of functions upon which my oscillation and burst detection functions are built.

Questions? Contact timcwhalen@gmail.com