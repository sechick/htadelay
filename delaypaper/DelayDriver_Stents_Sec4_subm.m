% Routine currently in use to produce plots for the paper. Includes both
% definition of stopping boundaries and simulations

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Set up so matlab's path includes the core routines for the delayed
%%%%%% sample optimizations. Clean up the matlab environment. Define
%%%%%% parameters which set up the PDE and MC accuracy. Please set up
%%%%%% options in this section based on what works best for your
%%%%%% experimental purposes. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% STENTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic
if ~exist('fignum','var'), fignum = 20; end;
[basic, advanced] = SetIllustration_SECTION4();
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


numinsimFreqDeltaVec=350;  % can change the number of points in freqdeltavec. The 11/3 in the next line is to get \pm 11/3 standard errors in width for vector's range
advanced.numinsimFreqDeltaVec=400; 
%advanced.simFreqDeltaVec = (basic.sigma/sqrt(basic.t0)) * (11/3) * (-ceil(numinsimFreqDeltaVec/2):ceil(numinsimFreqDeltaVec/2)) / ceil(numinsimFreqDeltaVec/2);  %this is now overwritten within the file "DoSectionFivePlotsStents

graphicsuffix = '';
PathReps = 10;      % number of replications for generating sample paths for demonstration
[fignum, mat] = DoSectionFourPlots(fignum, basic,advanced,PathReps,ProductionReps,graphicsuffix);
%fignum = DelaySimSurvival( fignum, basic, advanced, mat ) ;
minutestocompute = toc/60

