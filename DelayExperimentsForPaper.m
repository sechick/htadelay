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
doSaveMatFile = false;

%%%% the following can be customized to increase to to decrease accuracy of
%%%% both test run mode and production run mode. test run is fast but not
%%%% accurate, production mode is for generating higher resolution graphs.
PRODUCTIONREPS = 15000;
PRODUCTIONNUMBERSTD = 200;
PRODUCTIONNUMBERSTD = 120;

TESTREPS = 200;              % Allow end user to configure number of simulation replications
TESTNUMPERSTD = 30;          % and fineness of grid for PDE computations

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
if ~exist('fignum','var'), fignum = 20; end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Plots for Section 3: Illustration of features of optimal policy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% RUNTIMES:
% * On Win 8.1 x64 machine with 2.9Ghz clock, 8Gb RAM:
%   > when doProductionRuns = false
%       - First 3 clusters of code take 1.7-2min to compute, each
%       - cluster for figures 6 and 7 requires 2.7 minutes
%       - cluster for figure 5 requires 1.9 minutes
%   > when doProductionRuns = true and PRODUCTIONNUMBERSTD=100
%       - First cluster of code takes 41 min to compute
%
%•	DelayDriver_Stents_Sec4_subm : plots fig 2a and 3
clear basic; clear advanced; clear mat;
DelayDriver_Stents_Sec4_subm
movefile('.\Figure\paths.eps','.\paths.eps')
movefile('.\Figure\diffs_baseline.eps','.\diffs_baseline.eps')
movefile('.\Figure\pr_select_best.eps','.\pr_select_best.eps')
if doSaveMatFile
    save('StentsFig2a3.mat');
end

pause; close all hidden;

%•	DelayDriver_Stents_Sec4_comp_subm: plots fig 2b
clear basic; clear advanced; clear mat;
DelayDriver_Stents_Sec4_comp_subm
movefile('.\Figure\paths.eps','.\paths_comp.eps')
if doSaveMatFile
    save('StentsFig2b.mat');
end

pause; close all hidden;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Plots for Section 4: example for drug eluting stents
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pause; close all hidden;

%•	 DelayDriver_Stents_subm: plots fig 4
clear basic; clear advanced; clear mat;
DelayDriver_Stents_subm
movefile('.\Figure\pathsStents.eps','.\pathsStents.eps')
movefile('.\Figure\diffsStents.eps','.\diffsStents.eps')
movefile('.\Figure\E_num_seenStents.eps','.\E_num_seenStents.eps')
movefile('.\Figure\pr_select_bestStents.eps','.\pr_select_bestStents.eps')
if doSaveMatFile
    save('StentsFig4.mat');
end

pause; close all hidden;

%•	Drive_stents_sensit_tau: plots fig 6 and 7
clear basic; clear advanced; clear mat;
Drive_stents_sensit_tau
movefile('.\Figure\FigureOptVsOneshot.eps','.\FigureOptVsOneshot.eps')
movefile('.\Figure\FigureNumStarted.eps','.\FigureNumStarted.eps')
movefile('.\Figure\FigureDiffValueFunctionRR.eps','.\FigureDiffValueFunctionRR.eps')
movefile('.\Figure\FigureBounds.eps','.\SensitTauFigureBounds.eps')
movefile('.\Figure\FigureBoundsLessTau.eps','.\SensitTauFigureBoundsLessTau.eps')
if doSaveMatFile
    save('StentsFig67.mat');
end

pause; close all hidden;

%•	Drive_stents_sensit_c: plots fig 5
clear basic; clear advanced; clear mat;
Drive_stents_sensit_c
movefile('.\Figure\FigureBounds.eps','.\FigureBounds.eps')
if doSaveMatFile
    save('StentsFig5.mat');
end

pause; close all hidden;

%•	DelayColors: alternative method for producing figures 2 and 3 (illustration)
% in paper, in color, for example for use in seminar presentations
clear basic; clear advanced; clear mat;
doProductionRuns = true;
DelayColors
%movefile('.\Figure\FigureBounds.eps','.\FigureBounds.eps')
movefile('.\Figure\paths_colour.eps','.\paths_colour.eps')
movefile('.\Figure\paths_comp_colour.eps','.\paths_comp_colour.eps')
if doSaveMatFile
    save('DelayAlternative.mat');
end


