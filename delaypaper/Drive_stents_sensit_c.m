% For delay health tech assessment optimal stopping project with Martin, Paolo, Steve.  
% This script demonstrates the calling conventions for creating problem
% structures, doing an analysis on those structures, and then plotting the
% results. Cut and paste sections of this file in order to run the
% analysis.
%
% Adapted from DelaySimExperiment. Used for generating a bunch of plots
% where the analysis shows how 'tau' can influence the trial design.
%
% RUNS:
%   SINUS: Varies tau over several values.
%   
%  Date: 2015 04 06 SEC
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%          TESTS DIFFERENT VALUES OF c WITH STENTS              %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic

if ~exist('fignum','var'), fignum = 20; end;
[basic, advanced] = SetStents_sensit();
[basic, advanced] =  DelayStructureCompute (basic, advanced);

advanced.dirstring = 'Stents_sensit_c';
advanced.graphicextension = 'eps';
if doProductionRuns
    ProductionReps = PRODUCTIONREPS;  %15000
    advanced.MinGridPerStdev = PRODUCTIONNUMBERSTD; %200;
    advanced.NumGridsForContours = 70; 
else
    ProductionReps = TESTREPS;
    advanced.MinGridPerStdev = TESTNUMPERSTD; %20;
    advanced.NumGridsForContours = 40; 
end

advanced.NumPointsQuadrature = 100;
advanced.saveplot = true;           % set to true to save plots in files, false to not save files automatically
basic.online = false;
[basic, advanced, rval, msgs] = DelayInputValidator( basic, advanced );

numstd=5;
advanced.numinsimFreqDeltaVec = ceil(1 + numstd * (basic.sigma/sqrt(basic.t0))/advanced.dmu);
advanced.simFreqDeltaVec = (-advanced.numinsimFreqDeltaVec:advanced.numinsimFreqDeltaVec) * advanced.dmu; % SEC: to avoid some interpolation issues with mismatch between dmu grid and simFreqDelta grid
advanced.DOPLOT = false;

if ~rval, msgs, end;
fieldname = 'c';       % pick the name of the field whose values are to be changed
basicflag = true;           % set to true if the field is in the 'basic' structure, false if it is in the 'advanced' structure
%fieldvec = [0 20 40 80 250 500 1000];    % create an array of values to rotate through for that field 
fieldvec = [200 5000];    % create an array of values to rotate through for that field 
subtitle = 'Vary delay \tau';    % give a short subtitle name for the figures
fmodifier = fieldname;      % used to diferentiate file name if it is used for saving plots
[fignum, basicvec, advancedvec, legendvec, matvec] = TestDelayIterate_sensit_tau(fignum,basic,advanced,basicflag,fieldname,fieldvec,subtitle,fmodifier);
%[fignum] = UtilExperimentVectorPlot( fignum, basicvec, advancedvec, legendvec, matvec, subtitle,fmodifier);
%save('Stents_sensit_c.mat') ;

MULOW =  advancedvec(1).plot_lower;
MUHIGH = advancedvec(1).plot_upper;
VALDIFMAX = 2e7;

sameheightforallmaxloss = false;
for i=1:length(advancedvec)
    advanced.dirstring = 'Stents_sensit_c';
    advancedvec(i).filestring='stentzoom';
end
Sensitivity_SEC_5_c()

toc/60

