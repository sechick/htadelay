% Routine currently in use to produce plots for the paper. Includes both
% definition of stopping boundaries and simulations

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Set up so matlab's path includes the core routines for the delayed
%%%%%% sample optimizations. Clean up the matlab environment. Define
%%%%%% parameters which set up the PDE and MC accuracy. Please set up
%%%%%% options in this section based on what works best for your
%%%%%% experimental purposes. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all hidden; clear;

% bring in files from current directory and from two other locations,
%   delaycore\ contains routines which are generally useful
%   delaypaper\ contains routines specific to this paper
LocalDelaySetPaths;  % Default case: Sets up PATH of Matlab to get code for running. This file may be localized.
%LocalDelaySetPaths('mydirectory') % Example usage if you have installed
%delaycore and delaypaper directories into the directory 'mydirectory' on
%your machine. Usually, LocalDelaySetPaths should not need this
%localization, and the default usage (calling with no aguments, should be
%adequate.

pause on      % set to 'on' if you want extra pauses at certain times during computation to view things, to 'off' otherwise
pause off

% CAN SET FOLLOWING VALUES DEPENDING ON ACCURACY IN PLOTS DESIRED.
doProductionRuns = true;    % TO BE SET BY END USER: Set to true if lots of replications desired, false if short test is desired.
doProductionRuns = false;    % TO BE SET BY END USER: Set to true if lots of replications desired, false if short test is desired.
%%%% the following can be customized to increase to to decrease accuracy of
%%%% both test run mode and production run mode. test run is fast but not
%%%% accurate, production mode is for generating higher resolution graphs.
TESTREPS = 100;              % Allow end user to configure number of simulation replications
PRODUCTIONREPS = 15000;
TESTNUMPERSTD = 50;          % and fineness of grid for PDE computations
PRODUCTIONNUMBERSTD = 150;

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all hidden;
clear;
LocalDelaySetPaths;

TESTREPS = 100;              % Allow end user to configure number of simulation replications
PRODUCTIONREPS = 15000;
TESTNUMPERSTD = 30;          % and fineness of grid for PDE computations
PRODUCTIONNUMBERSTD = 150;

doProductionRuns = false;    % TO BE SET BY END USER: Set to true if lots of replications desired, false if short test is desired.

pause on
pause off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Plots for Section 3: Illustration of features of optimal policy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%•	DelayDriver_Stents_Sec4_subm : plots fig 2a and 3
clear basic; clear advanced; clear mat;
DelayDriver_Stents_Sec4_subm
movefile('.\Figure\paths.eps','.\paths.eps')
movefile('.\Figure\diffs_baseline.eps','.\diffs_baseline.eps')
movefile('.\Figure\pr_select_best.eps','.\pr_select_best.eps')
save('StentsFig2a3.mat');

%•	DelayDriver_Stents_Sec4_comp_subm: plots fig 2b
clear basic; clear advanced; clear mat;
DelayDriver_Stents_Sec4_comp_subm
movefile('.\Figure\paths.eps','.\paths_comp.eps')
save('StentsFig2b.mat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Plots for Section 4: example for drug eluting stents
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%•	 DelayDriver_Stents_subm: plots fig 4
clear basic; clear advanced; clear mat;
DelayDriver_Stents_subm
movefile('.\Figure\pathsStents.eps','.\pathsStents.eps')
movefile('.\Figure\diffsStents.eps','.\diffsStents.eps')
movefile('.\Figure\E_num_seenStents.eps','.\E_num_seenStents.eps')
movefile('.\Figure\pr_select_bestStents.eps','.\pr_select_bestStents.eps')
save('StentsFig4.mat');

%•	Drive_stents_sensit_tau: plots fig 6 and 7
clear basic; clear advanced; clear mat;
Drive_stents_sensit_tau
movefile('.\Figure\FigureOptVsOneshot.eps','.\FigureOptVsOneshot.eps')
movefile('.\Figure\FigureNumStarted.eps','.\FigureNumStarted.eps')
movefile('.\Figure\FigureDiffValueFunctionRR.eps','.\FigureDiffValueFunctionRR.eps')
save('StentsFig67.mat');

%•	Drive_stents_sensit_c: plots fig 5
clear basic; clear advanced; clear mat;
Drive_stents_sensit_c
movefile('.\Figure\FigureBounds.eps','.\FigureBounds.eps')
save('StentsFig5.mat');

MISSING SetStents_sensit
