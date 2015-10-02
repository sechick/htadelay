%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This gives an alternate version for computing figures 2 and 3 in main
% paper, and for producing graphs which display probability of 'reversing'
% the decision (not same decision after all samples observed as compared to
% decision which is optimal on basis of only samples seen through time of
% decision to stop).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Set up data for Section 4 graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
if ~exist('fignum','var'), fignum = 20; end;

[basic, advanced] = SetIllustration_SECTION4();
if doProductionRuns
    ProductionReps = TESTREPS; % NB: don't do lots of sims for the color plots... PRODUCTIONREPS;
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
graphicsuffix = '_colour';
PathReps = 10;      % number of replications for generating sample paths for demonstration
[fignum, mat] = DoSectionFourPlots(fignum, basic,advanced,PathReps,ProductionReps,graphicsuffix,true);
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Run graphics for comparator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
graphicsuffix = '_comp_colour';
basiccomp = basic;
basiccomp.tau = 500;
advancedcomp = advanced;
PathReps = 10;      % number of replications for generating sample paths for demonstration
[fignum, matcomp] = DoSectionFourPlots(fignum, basiccomp,advancedcomp,PathReps,ProductionReps,graphicsuffix,true);

minutestocompute = toc/60

