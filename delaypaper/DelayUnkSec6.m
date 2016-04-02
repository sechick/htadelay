close all; clear;

LocalDelaySetPaths % Edit a local copy of this MACRO, don't check in on top of installed file.

pause on      % set to 'on' if you want extra pauses at certain times during computation to view things, to 'off' otherwise
pause off

% CAN SET FOLLOWING VALUES DEPENDING ON ACCURACY IN PLOTS DESIRED.
doProductionRuns = false;    % TO BE SET BY END USER: Set to true if lots of replications desired, false if short test is desired.
doProductionRuns = true;    % TO BE SET BY END USER: Set to true if lots of replications desired, false if short test is desired.
doSaveMatFile = false;

%%%% the following can be customized to increase to to decrease accuracy of
%%%% both test run mode and production run mode. test run is fast but not
%%%% accurate, production mode is for generating higher resolution graphs.
PRODUCTIONREPS = 15000;
PRODUCTIONREPS = 15;
PRODUCTIONNUMBERSTD = 200;
PRODUCTIONNUMBERSTD = 60;

TESTREPS = 200;              % Allow end user to configure number of simulation replications
TESTNUMPERSTD = 40;          % and fineness of grid for PDE computations

%%%%%%%%%% SET UP THE EXPERIMENT %%%%%%%%%%%%%

advanced.saveplot = true;           % set to true to save plots in files, false to not save files automatically
advanced.DOPLOT= false ;
advanced.keepAllOutput = true ;
if doProductionRuns
    ProductionReps = PRODUCTIONREPS;  %15000
    ProductionGridPerStdev = PRODUCTIONNUMBERSTD; %200;
else
    ProductionReps = TESTREPS;
    ProductionGridPerStdev = TESTNUMPERSTD; %20;
end

%%%%%%%% TEST KNOWN VARIANCE PDE VERSUS UNKNOWN VARIANCE WITH KG*-TYPE COMPUTATION %%%%%%%%
if ~exist('fignum','var'), fignum = 20; end;
[basic, advanced] = SetStents_Unk();
advanced.filestring='UnkUvarKGs'
advanced.DOPDE = true;              % will vary dopde and unkvariance fields below
advanced.UnkVariance = false;
advanced.NumGridsForContours = ProductionGridPerStdev;
advanced.simNumReps = ProductionReps;    
%basic.theta = 1;
[basic, advanced, rval, msgs] = DelayInputValidator( basic, advanced );
if ~rval, msgs, end;
fieldname = 'UnkVariance';       % pick the name of the field whose values are to be changed
basicflag = false;           % set to true if the field is in the 'basic' structure, false if it is in the 'advanced' structure
fieldvec = [false true];    % create an array of values to rotate through for that field
subtitle = 'PDE (known var) and KG* (unknown var)';    % give a short subtitle name for the figures
advanced.verbose = true;
fmodifier = 'PDEKGunk';      % used to diferentiate file name if it is used for saving plots
if advanced.verbose % the 'true' part of this keeps the various variables around. the 'else' part has it implemented in a subroutine, and thus does not keep the output for futher use
    % create vectors of basic and advanced structures
    [basicvec, advancedvec, legendvec, veclen] = UtilExperimentVectorCreate( basic, advanced, basicflag, fieldname, fieldvec);
    % finish setting up the parameters for the runs
    i=2;        % for second item in vector of experiments
    advancedvec(i).UnkVariance = true;      % set unknown variance to true - assumes that the t0 and other params had been set already
    advancedvec(i).DOPDE = false;           % do KG* rather than PDE for stopping rule computation    
    % run the analysis for each of the structures
    for i=1:veclen
        [~, mat] = DelayCurvesRecur(basicvec(i), advancedvec(i));        % Do stage II of DP
        [mat] = DelayStageOne(basicvec(i), advancedvec(i), mat );           % Do stage I of DP
%        if advancedvec(i).UnkVariance 
            [ fignum, mat ] = DelaySimOverviewUnkPaper( fignum, basicvec(i), advancedvec(i),mat );
%        else
%            [ fignum, mat ] = DelaySimOverview( fignum, basicvec(i), advancedvec(i), mat );
%        end
        fignum = DelaySimOutput( fignum, basicvec(i), advancedvec(i), mat );     % requires DelaySimOverview to have been called
        if i==1 matvec = mat; else matvec = [matvec mat]; end
    end
    % do the plots based on the structures
    [fignum] = UtilExperimentVectorPlot( fignum, basicvec, advancedvec, legendvec, matvec, subtitle, fmodifier);
    for i=1:veclen
        advancedvec(i).subtitle = legendvec{i};
        advancedvec(i).filestring = [advancedvec(i).filestring num2str(i)];
%        fignum = DelayDoContours(fignum, basicvec(i), advancedvec(i), matvec(i));        % Generate a bunch of contour plots
    end
else
    [fignum, basicvec, advancedvec, legendvec, matvec] = TestDelayIterate(fignum,basic,advanced,basicflag,fieldname,fieldvec,subtitle,fmodifier);
    [fignum] = UtilExperimentVectorPlot( fignum, basicvec, advancedvec, legendvec, matvec, subtitle,fmodifier);
end

save Unk15script.mat
