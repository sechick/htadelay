# Sequential Experimentation with Delay.

ALL RIGHTS RESERVED BY AUTHORS. Code provided as is for noncommercial, academic usage only.

This repository contains Matlab code which implements several contributions from the paper:

'A Bayesian Decision-Theoretic Model of Sequential Experimentation with Delayed Response'

by Stephen E. Chick, Martin Forster, and Paolo Pertile (accepted to appear in JRSS B, online since 2017 at http://onlinelibrary.wiley.com/doi/10.1111/rssb.12225/abstract), with preprint of earlier draft available at https://ideas.repec.org/p/yor/yorken/15-09.html.

The code computes optimal stopping times for sequential experiments when samples are observed with delay. The code includes:

 - constructor and validator methods for creating problem instances
 - computation of the optimal PDE solution for stage II sampling, and for stage I sampling
 - drivers for monte carlo simulations
 - the code computes quite a few other things which are not included in the paper.

 Analysis assumes a known sampling variance and other assumptions from the paper.

# Workflow

1. Download the matlab code from this repo. When installed locally, there should be files DelayDriver.m and two directories, 'delaycore' and 'delaypaper', and potentially several other files.
2. [Optional] Update the LocalDelaySetPaths.m file to contain local directory information to set the path to the matlab code for the computations if you intend to run the code in a nonstandard directory structure.
3. Run matlab and set the current directory to the local copy of the repo you have created (e.g. /htadelay-master).
4. Open the file DelayExperimentsForPaper.m. It contains chunks of code which generates the plots in the paper, and many more related plots.
 -- If you edit a few parameters in this file, you can change the accuracy of the PDE derivations of optimal stopping boundaries and of other quantities, as well as of Monte Carlo results.
5. Open the file DelayDriver.m. It contains a sequence of chunks of code which illustrate how the interface to the routines works.
 -- Copy the chunks of code identified in DelayDriver.m into the matlab console.
6. Enjoy the plots.

# Comment on simulations

1. Smooth curves (potentially with a few jumps or discontinuties) in plots are computed from the free boundary PDE of the continuous time approximation.
2. Jagged curves represent output from Monte Carlo (MC) results.
    a. MC results labeled 'Bayesian' means that sample paths were generated by assuming the given value of mu0 on the plot is the prior mean for sampling unknown means W_i from the prior distrubution.
    b. MC results label 'Frequentist' or 'Freq.' means that sample paths were run by fixing all W_i to the value of mu0 on the axis
    The values of mu0 are then varied for both of those types of plots.

# Updates

1. 2018 04 13: Updated certain plotting functions. As of Matlab 2017a, the 'legend' functionality has changed. Code was updated to avoid certain spurious entries in the legends of some of the graphs. Need to include 'AutoUpdate', 'Off' in legends to avoid the unwanted legend entries.

Some code under development is built to allow for certain functionality when assuming sampling variances are unknown. This is in progress....
