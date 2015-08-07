% Routine currently in use to produce plots for the paper. Includes both
% definition of stopping boundaries and simulations

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Set up so matlab's path includes the core routines for the delayed
%%%%%% sample optimizations. Clean up the matlab environment. Define
%%%%%% parameters which set up the PDE and MC accuracy. Please set up
%%%%%% options in this section based on what works best for your
%%%%%% experimental purposes. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%delaycoredir = 'd:\users\papers\Forster\hta\trunk\delaycore\'
%delaycoredir = 'c:\martin\chick\hta\trunk\delaycore\' ;
%addpath(genpath(delaycoredir))
%addpath(genpath('..\delaycore\'))
LocalDelaySetPaths;  % Default case: Sets up PATH of Matlab to get code for running. This file may be localized.
%LocalDelaySetPaths('mydirectory') % Example usage if you have installed
%delaycore and delaypaper directories into the directory 'mydirectory' on
%your machine. Usually, LocalDelaySetPaths should not need this
%localization, and the default usage (calling with no aguments, should be
%adequate.
close all hidden
clear 

% CAN SET FOLLOWING VALUES DEPENDING ON ACCURACY IN PLOTS DESIRED.
doProductionRuns = true;    % TO BE SET BY END USER: Set to true if lots of replications desired, false if short test is desired.
doProductionRuns = false;    % TO BE SET BY END USER: Set to true if lots of replications desired, false if short test is desired.

pause on
pause off

TESTREPS = 100;              % Allow end user to configure number of simulation replications
PRODUCTIONREPS = 15000;
TESTNUMPERSTD = 50;          % and fineness of grid for PDE computations
PRODUCTIONNUMBERSTD = 150;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This sets up the correct number of replications for the several values

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Set up data for Section 4 graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
if ~exist('fignum','var'), fignum = 20; end;

[basic, advanced] = SetIllustration_SECTION4();
if doProductionRuns
    ProductionReps = PRODUCTIONREPS;
    advanced.MinGridPerStdev = PRODUCTIONNUMBERSTD;
else
    ProductionReps = TESTREPS;
    advanced.MinGridPerStdev = TESTNUMPERSTD;
end
basic.tau = 1000;
advanced.saveplot = false;           % set to true to save plots in files, false to not save files automatically
advanced.DOPLOT= false ;
%basic.PPatients = 10000;
advanced.keepAllOutput = true ;

numinsimFreqDeltaVec=280;  % can change the number of points in freqdeltavec. The 11/3 in the next line is to get \pm 11/3 standard errors in width for vector's range
advanced.simFreqDeltaVec = (basic.sigma/sqrt(basic.t0)) * (11/3) * (-ceil(numinsimFreqDeltaVec/2):ceil(numinsimFreqDeltaVec/2)) / ceil(numinsimFreqDeltaVec/2); % MF 27/03 changed to 15/03 because otherwise the range was too short (did not extend to the optimal `do nothing' values of mu0)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Run graphics for base case
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
graphicsuffix = '';
PathReps = 10;      % number of replications for generating sample paths for demonstration
[fignum, mat] = DoSectionFourPlots(fignum, basic,advanced,PathReps,ProductionReps,graphicsuffix);
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Run graphics for comparator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
graphicsuffix = '_comp';
basiccomp = basic;
basiccomp.tau = 500;
advancedcomp = advanced;
PathReps = 10;      % number of replications for generating sample paths for demonstration
[fignum, matcomp] = DoSectionFourPlots(fignum, basiccomp,advancedcomp,PathReps,ProductionReps,graphicsuffix);

fignum = DoReversalPlot(fignum,basic,advanced,mat,basiccomp,advancedcomp,matcomp,graphicsuffix);
minutestocompute = toc/60

save sectionfour


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Set up data for Section 5's graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


DelayDriver_Stents_subm.m








