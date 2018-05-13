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

%ORIGIN OF THIS FILE IS FILE "DelayExperimentsForPaper". 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% STENTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic
if ~exist('fignum','var'), fignum = 20; end;
[basic, advanced] = SetStents_SECTION5();
advanced.dirstring = 'stents';

advanced.saveplot = true;           % set to true to save plots in files, false to not save files automatically
advanced.DOPLOT= false ;
advanced.keepAllOutput = true ;
if doProductionRuns
    ProductionReps = PRODUCTIONREPS;  %15000
    advanced.MinGridPerStdev = PRODUCTIONNUMBERSTD; %200;
else
    ProductionReps = TESTREPS;
    advanced.MinGridPerStdev = TESTNUMPERSTD; %20;
end

%numinsimFreqDeltaVec=350;  % can change the number of points in freqdeltavec. The 11/3 in the next line is to get \pm 11/3 standard errors in width for vector's range
%advanced.numinsimFreqDeltaVec=400;  % 200 probably enough, or 350 is more than enough
advanced.numinsimFreqDeltaVec=200;  % 200 probably enough, or 350 is more than enough
%advanced.simFreqDeltaVec = (basic.sigma/sqrt(basic.t0)) * (11/3) * (-ceil(numinsimFreqDeltaVec/2):ceil(3/2)) / ceil(numinsimFreqDeltaVec/2);  %this is now overwritten within the file "DoSectionFivePlotsStents

graphicsuffix = '';
PathReps = 10;      % number of replications for generating sample paths for demonstration
[fignum, mat] = DoSectionFivePlotsStents(fignum, basic,advanced,PathReps,ProductionReps,graphicsuffix);
%fignum = DelaySimSurvival( fignum, basic, advanced, mat ) ;
minutestocompute = toc/60


