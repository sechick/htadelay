close all; 
clear;
tic
LocalDelaySetPaths % Edit a local copy of this MACRO, don't check in on top of installed file.

pause on      % set to 'on' if you want extra pauses at certain times during computation to view things, to 'off' otherwise
pause off

% CAN SET FOLLOWING VALUES DEPENDING ON ACCURACY IN PLOTS DESIRED.
doProductionRuns = false;    % TO BE SET BY END USER: Set to true if lots of replications desired, false if short test is desired.
%doProductionRuns = true;    % TO BE SET BY END USER: Set to true if lots of replications desired, false if short test is desired.
doSaveMatFile = false;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Set up data for Section 5's graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
if ~exist('fignum','var'), fignum = 20; end;

%%%%%%%%%% SET UP THE EXPERIMENT %%%%%%%%%%%%%

advanced.saveplot = true;           % set to true to save plots in files, false to not save files automatically
advanced.DOPLOT= false ;
advanced.keepAllOutput = true ;


PRODUCTIONREPS = 15000;
PRODUCTIONNUMBERSTD = 200;
TESTREPS = 200;              % Allow end user to configure number of simulation replications
TESTNUMPERSTD = 40;          % and fineness of grid for PDE computations

if doProductionRuns
    ProductionReps = PRODUCTIONREPS;  %15000
    PathReps = PRODUCTIONNUMBERSTD; %200;
else
    ProductionReps = TESTREPS;
    PathReps = TESTNUMPERSTD; %20;
end


[basic, advanced] = SetStents_Unk();

advanced.numinsimFreqDeltaVec=300;  % can change the number of points in freqdeltavec. the 11/3 in the next line is to get \pm 11/3 standard errors in width for vector's range
advanced.simFreqDeltaVec = (basic.sigma/sqrt(basic.t0)) * (11/3) * (-ceil(advanced.numinsimFreqDeltaVec/2):ceil(advanced.numinsimFreqDeltaVec/2)) / ceil(advanced.numinsimFreqDeltaVec/2);
%advanced.saveplot = false;           % set to true to save plots in files, false to not save files automatically
advanced.MinGridPerStdev = 120; % MF 28/02: increased this to 120 to increase accuracy of free boundary calcs and to show more dots for Stage 0
advanced.DOPLOT= false ;
%basic.PPatients = 10000;
advanced.keepAllOutput = true;



advanced.numinsimFreqDeltaVec=300;  % can change the number of points in freqdeltavec. the 11/3 in the next line is to get \pm 11/3 standard errors in width for vector's range
advanced.simFreqDeltaVec = (basic.sigma/sqrt(basic.t0)) * (11/3) * (-ceil(advanced.numinsimFreqDeltaVec/2):ceil(advanced.numinsimFreqDeltaVec/2)) / ceil(advanced.numinsimFreqDeltaVec/2);

[basic, advanced] = SetStents_Unk();
advanced.DOPDE = false;              % will vary dopde and unkvariance fields below
advanced.UnkVariance = true;
advanced.NumGridsForContours = PathReps;
graphicsuffix = '_prova';
%advanced.simNumReps = ProductionReps;    

[fignum, mat] = DoSectionFivePlotsStentsUnk(fignum, basic,advanced,PathReps,ProductionReps,graphicsuffix);

toc/60
