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

%ORIGIN OF THIS FILE IS FILE "DelayExperimentsForPaper". This is my version

LocalDelaySetPaths;

close all hidden
clear 

TESTREPS = 100;              % Allow end user to configure number of simulation replications
PRODUCTIONREPS = 15000; % was 15000
TESTNUMPERSTD = 50;          % and fineness of grid for PDE computations
PRODUCTIONNUMBERSTD = 120; % was 150

doProductionRuns = true ;    % TO BE SET BY END USER: Set to true if lots of replications desired, false if short test is desired.

pause on
pause off


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% HIP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic
if ~exist('fignum','var'), fignum = 20; end;
[basic, advanced] = SetHip_OSM();
advanced.dirstring = 'hip';

advanced.saveplot = true;           % set to true to save plots in files, false to not save files automatically
advanced.DOPLOT= false ;
advanced.keepAllOutput = true ;
if doProductionRuns
    ProductionReps = 15000; 
    advanced.MinGridPerStdev = 200;
else
    ProductionReps = 150;
    advanced.MinGridPerStdev = 20;
end

numinsimFreqDeltaVec=350;  % can change the number of points in freqdeltavec. The 11/3 in the next line is to get \pm 11/3 standard errors in width for vector's range
advanced.numinsimFreqDeltaVec=400; 
%advanced.simFreqDeltaVec = (basic.sigma/sqrt(basic.t0)) * (11/3) * (-ceil(numinsimFreqDeltaVec/2):ceil(numinsimFreqDeltaVec/2)) / ceil(numinsimFreqDeltaVec/2);  %this is now overwritten within the file "DoSectionFivePlotsStents

graphicsuffix = '';
PathReps = 10;      % number of replications for generating sample paths for demonstration
[fignum, mat] = DoSectionFivePlotsHip(fignum, basic,advanced,PathReps,ProductionReps,graphicsuffix);
%fignum = DelaySimSurvival( fignum, basic, advanced, mat ) ;
minutestocompute = toc/60
save('Hip.mat') ;

