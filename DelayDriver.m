% For delay health tech assessment optimal stopping project with Martin, Paolo, Steve.  
% This script demonstrates the calling conventions for creating problem
% structures, doing an analysis on those structures, and then plotting the
% results. Cut and paste sections of this file in order to run the
% analysis.
%
% Code is 'as is' and no guarantees for correctness.  
%
% (c) 2014, S Chick
% Created: 17 April 2014
% Last touched: 17 April 2014
% 
% NEXT STEPS?
% (*) Do functions for more frequentist statistics (be it PDE or be it
% via simulation).
% (*) Do parking lot stuff to improve API/GUI after we know better what is 
% functionality we can do, and what end user usage cases might look like.
% (*) Fix apparent bug for Prob pick new / prob hit top when tau is quite large
% (*) See README-TODO.txt
%
% FIX: PARKING LOT 'TO DO'
% (*) Fix the list of parameters which are created by other routines for the
% advanced and mat structures.
%

%'note: if DOPLOT is set to true, then NumPlots should be small '
%'otherwise it takes long to run, but then the contour plots are '
%'a bit wiggly. For smoother/better contour plots, set NumPlots '
%'bigger and set DOPLOT to false'
%'Also, move plots 4, 5, 6, 7 to different parts of screen if DOPLOT is on'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % SECTION 1: Standard example analysis flow with several
            % problem examples

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

LocalDelaySetPaths % Edit a local copy of this MACRO, don't check in on top of installed file.

% Test of calling structure using the optional arguments to
% DelayInputConstructor. The following two example options result 
% in identical values for basic and advanced.
% OPTION 1:
[basic, advanced] = DelayInputConstructor();
basic.c = 2.1; basic.theta = 0.999;
advanced.filestring = 'test'; advanced.dirstring = 'wow'; advanced.simNumReps = 500; 
advanced.saveplot = 1; advanced.DOPLOT = false;
% OPTION 2:
somebasicflags = {'c', 2.1, 'theta', 0.999};
someadvancedflags = { 'filestring', 'test', 'dirstring', 'wow', 'simNumReps', 500, 'saveplot', 1, 'DOPLOT', false };
[basic, advanced] = DelayInputConstructor(somebasicflags,someadvancedflags);
% now, validate the structure
[basic, advanced, rval, msgs] = DelayInputValidator( basic, advanced )
if rval
    [~, mat] = DelayCurvesRecur(basic, advanced);
    [mat] = DelayStageOne(basic, advanced, mat ); 
    if ~exist('fignum','var'), fignum = 20; end;
    [fignum, mat] = DelayCurvesTheoretical(fignum, basic, advanced, mat );
    fignum = DelayDoContours(fignum, basic, advanced, mat);        % Generate a bunch of contour plots
    fignum = DelayPlotBayesInfo( fignum, basic, advanced, mat );
    [ fignum, mat ] = DelaySimOverview( fignum, basic, advanced, mat );
    fignum = DelaySimOutput( fignum, basic, advanced, mat );     % requires DelaySimOverview to have been called
else
    msgs
end


% Test of calling structure using the optional arguments to check
% computed boundaries with theoretical boundaries
[basic, advanced] = DelayInputConstructor();
basic.online = 0;
basic.c = 260; basic.theta = 0.999; basic.ICost = 200000; basix.TMax = 3000; basic.tau = 2;
basic.c =  0; basic.theta = 0.999; basic.ICost = 200000; basix.TMax = 3000; basic.tau = 2;
basic.c = 25; basic.theta = 1.000; basic.ICost = 100000; basix.TMax = 3000; basic.tau = 2;
basic.c = 80; basic.theta = 1.000; basic.ICost = 100000; basix.TMax = 3000; basic.tau = 2;
advanced.filestring = 'test'; advanced.dirstring = 'wow'; advanced.simNumReps = 500; advanced.saveplot = 1; advanced.DOPLOT = false;
advanced.MinGridPerStdev = 30;
% now, validate the structure
[basic, advanced, rval, msgs] = DelayInputValidator( basic, advanced );
shiftfromI = basic.ICost/basic.PPatients
shiftfromc = basic.c/(1-basic.theta)/basic.PPatients

[~, mat] = DelayCurvesRecur(basic, advanced);
[mat] = DelayStageOne(basic, advanced, mat ); 
if ~exist('fignum','var'), fignum = 20; end;
[fignum, mat] = DelayCurvesTheoretical(fignum, basic, advanced, mat );
advanced.dmu

%    fignum = DelayDoContours(fignum, basic, advanced, mat);        % Generate a bunch of contour plots
%    fignum = DelayPlotBayesInfo( fignum, basic, advanced, mat );
%    [ fignum, mat ] = DelaySimOverview( fignum, basic, advanced, mat );
%    fignum = DelaySimOutput( fignum, basic, advanced, mat );     % requires DelaySimOverview to have been called

% TEST
% OPTION 1:
[basic, advanced] = DelayInputConstructor();
basic.c = 2.1; basic.theta = 0.999;
advanced.filestring = 'test'; advanced.dirstring = 'wow'; advanced.simNumReps = 800; 
advanced.saveplot = 1; advanced.DOPLOT = false;
advanced.MinGridPerStdev = 25; advanced.NumGridsForContours = 50;
[basic, advanced, rval, msgs] = DelayInputValidator( basic, advanced )
if rval
    [~, mat] = DelayCurvesRecur(basic, advanced);
    [mat] = DelayStageOne(basic, advanced, mat ); 
    if ~exist('fignum','var'), fignum = 20; end;
    [fignum, mat] = DelayCurvesTheoretical(fignum, basic, advanced, mat );
    fignum = DelayDoContours(fignum, basic, advanced, mat);        % Generate a bunch of contour plots
    fignum = DelayPlotBayesInfo( fignum, basic, advanced, mat );
    [ fignum, mat ] = DelaySimOverview( fignum, basic, advanced, mat );
    fignum = DelaySimOutput( fignum, basic, advanced, mat );     % requires DelaySimOverview to have been called
else
    msgs
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% A standard example of calling routine
if ~exist('fignum','var'), fignum = 20; end;
[basic, advanced] = SetBaseParameters();
basic.online = true; % will do ONLINE learning
%basic.theta = 0.999;
advanced.DOPLOT = false; % comment this out or set DOPLOT to false to suppress the dynamic graphs during the run
advanced.MinGridPerStdev = 30; advanced.NumGridsForContours = 60; 
advanced.DoRegretIntegral = true;       % false to do quadrature in regret estimation, true to do integration
advanced.DoRegretIntegral = false;       % false to do quadrature in regret estimation, true to do integration
%basic.tau = 1200;
%advanced.RegretPenalty = 0;
advanced.saveplot = true;           % set to true to save plots in files, false to not save files automatically
[basic, advanced, rval, msgs] = DelayInputValidator( basic, advanced );
if ~rval, msgs, end;
[~, mat] = DelayCurvesRecur(basic, advanced);
[mat] = DelayStageOne(basic, advanced, mat ); 

fignum = DelayPlotBayesInfo( fignum, basic, advanced, mat );

[fignum, mat] = DelayCurvesTheoretical(fignum, basic, advanced, mat );
fignum = DelayDoContours(fignum, basic, advanced, mat);        % Generate a bunch of contour plots


%advanced.keepAllOutput=true;
advanced.simNumReps = 800;
[ fignum, mat ] = DelaySimOverview( fignum, basic, advanced, mat );
fignum = DelaySimOutput( fignum, basic, advanced, mat );     % requires DelaySimOverview to have been called

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% An example with parameters for the hip arthro paper, mostly to get
% contour plots
if ~exist('fignum','var'), fignum = 20; end;
[basic, advanced] = SetHipArthro();
'note: fix the params of hip arthro, e.g. PPatients and ICost'
advanced.saveplot = true;           % set to true to save plots in files, false to not save files automatically
[basic, advanced, rval, msgs] = DelayInputValidator( basic, advanced );
if ~rval, msgs, end;
[~, mat] = DelayCurvesRecur(basic, advanced);
[mat] = DelayStageOne(basic, advanced, mat ); 
[fignum, mat] = DelayCurvesTheoretical(fignum, basic, advanced, mat );
fignum = DelayDoContours(fignum, basic, advanced, mat);        % Generate a bunch of contour plots
fignum = DelayPlotBayesInfo( fignum, basic, advanced, mat );

[ fignum, mat ] = DelaySimOverview( fignum, basic, advanced, mat );
fignum = DelaySimOutput( fignum, basic, advanced, mat );     % requires DelaySimOverview to have been called

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% A different example with parameters for the hip arthro paper, which does
% several plots, but not the contours
if ~exist('fignum','var'), fignum = 20; end;
[basic, advanced] = SetHipArthro();
'note: fix the params of hip arthro, e.g. PPatients and ICost'
basic.TMax = 200;   % try a larger value just for kicks.
advanced.MinGridPerStdev = 30;
advanced.saveplot = true;           % set to true to save plots in files, false to not save files automatically
[basic, advanced, rval, msgs] = DelayInputValidator( basic, advanced );
if ~rval, msgs, end;
[~, mat] = DelayCurvesRecur(basic, advanced);
[mat] = DelayStageOne(basic, advanced, mat ); 
[fignum, mat] = DelayCurvesTheoretical(fignum, basic, advanced, mat );

%fignum = DelayDoThreeTestPlots(fignum, basic, advanced);
fignum = DelayPlotBayesInfo( fignum, basic, advanced, mat );

[ fignum, mat ] = DelaySimOverview( fignum, basic, advanced, mat );
fignum = DelaySimOutput( fignum, basic, advanced, mat );     % requires DelaySimOverview to have been called

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% An example with 'test' parameters
if ~exist('fignum','var'), fignum = 20; end;
[basic, advanced] = SetTestExample();
[basic, advanced, rval, msgs] = DelayInputValidator( basic, advanced );
if ~rval, msgs, end;
[~, mat] = DelayCurvesRecur(basic, advanced);
[mat] = DelayStageOne(basic, advanced, mat ); 
fignum = DelayDoContours(fignum, basic, advanced, mat);        % Generate a bunch of contour plots
fignum = DelayPlotBayesInfo( fignum, basic, advanced, mat );

[ fignum, mat ] = DelaySimOverview( fignum, basic, advanced, mat );
fignum = DelaySimOutput( fignum, basic, advanced, mat );     % requires DelaySimOverview to have been called

fignum = DelayDoThreeTestPlots(fignum, basic, advanced);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% An example with 'stent' parameters
if ~exist('fignum','var'), fignum = 20; end;
[basic, advanced] = SetStents();
advanced.DOPLOT = 0;
%basic.tau = 10; 
advanced.saveplot = true;
[basic, advanced, rval, msgs] = DelayInputValidator( basic, advanced );
if ~rval, msgs, end;
[~, mat] = DelayCurvesRecur(basic, advanced);
[mat] = DelayStageOne(basic, advanced, mat ); 

% Test of length of optimal one-stage trial and its length, considering
% lengths which are found in the vector toneshot (can make that vector
% sparser if one wishes to improve speed, or one might try to optmize
% DelayOptimalBayesOneStage more.
toneshot = 5*(0:basic.TMax);
advanced.DOPLOT = 1;
[mat] = DelayOptimalBayesOneStage(basic, advanced, toneshot, mat ); 

% Do some monte carto simulations and print some output for it
[ fignum, mat ] = DelaySimOverview( fignum, basic, advanced, mat );
fignum = DelaySimOutput( fignum, basic, advanced, mat );     % requires DelaySimOverview to have been called

fignum = DelayDoContours(fignum, basic, advanced, mat);        % Generate a bunch of contour plots
fignum = DelayPlotBayesInfo(fignum, basic, advanced, mat);

fignum = DelayDoThreeTestPlots(fignum, basic, advanced);
fignum = UtilExperimentVectorPlot( fignum, basic, advanced, legend, mat, '', 'tst' )

%====
% An example with 'stent' parameters
if ~exist('fignum','var'), fignum = 20; end;
[basic, advanced] = SetStents_c0();
advanced.DOPLOT = 0;
basic.tau = 10;
[basic, advanced, rval, msgs] = DelayInputValidator( basic, advanced );
if ~rval, msgs, end;
[~, mat] = DelayCurvesRecur(basic, advanced);
[mat] = DelayStageOne(basic, advanced, mat ); 
[ fignum, mat ] = DelaySimOverview( fignum, basic, advanced, mat );
fignum = DelaySimOutput( fignum, basic, advanced, mat );     % requires DelaySimOverview to have been called

fignum = DelayDoContours(fignum, basic, advanced, mat);        % Generate a bunch of contour plots
fignum = DelayPlotBayesInfo(fignum, basic, advanced, mat);

%====
% An example with 'stent' parameters
if ~exist('fignum','var'), fignum = 20; end;
[basic, advanced] = SetStents_rho0();
advanced.DOPLOT = 0;
basic.tau = 10;
[basic, advanced, rval, msgs] = DelayInputValidator( basic, advanced );
if ~rval, msgs, end;
[~, mat] = DelayCurvesRecur(basic, advanced);
[mat] = DelayStageOne(basic, advanced, mat ); 
[ fignum, mat ] = DelaySimOverview( fignum, basic, advanced, mat );
fignum = DelaySimOutput( fignum, basic, advanced, mat );     % requires DelaySimOverview to have been called

fignum = DelayDoContours(fignum, basic, advanced, mat);        % Generate a bunch of contour plots
fignum = DelayPlotBayesInfo(fignum, basic, advanced, mat);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % SECTION 2: Call some functions which do tests calls, which
            % may be useful in debugging the code.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% NOTE: Old test functions no longer exist - they are obsolete with the new
% code. Need a new test function to demonstrate the terminal rewards.
% However, new data structure (basic, advanced) to describe problem and the
% computation means that the old calling convention for the old test
% function scripts is not longer valid. 
%
% Instead, use the following structures, which has a special set of
% routines to create vectors of numerical experiments with one parameter
% varied, for example tau, discount rate, etc etc. See the following
% examples.

%%%%%%%%% TESTS DIFFERENT VALUES OF TAU %%%%%%%%%%%%%%%
if ~exist('fignum','var'), fignum = 20; end;
[basic, advanced] = SetBaseParameters();
advanced.NumGridsForContours = 60; advanced.MinGridPerStdev = 20; % these are smaller than usual as we iterate multiple graphs
advanced.NumPointsQuadrature = 80;
advanced.saveplot = true;           % set to true to save plots in files, false to not save files automatically
advanced.dirstring = 'tau';
basic.online = false;
[basic, advanced, rval, msgs] = DelayInputValidator( basic, advanced );
if ~rval, msgs, end;
fieldname = 'tau';       % pick the name of the field whose values are to be changed
basicflag = true;           % set to true if the field is in the 'basic' structure, false if it is in the 'advanced' structure
fieldvec = 100*2.^(8/3:1/3:4);    % create an array of values to rotate through for that field 
fieldvec = [0 1 10 100 1000];    % create an array of values to rotate through for that field 
subtitle = 'Vary delay \tau';    % give a short subtitle name for the figures
fmodifier = fieldname;      % used to diferentiate file name if it is used for saving plots
[fignum, basicvec, advancedvec, legendvec, matvec] = TestDelayIterate(fignum,basic,advanced,basicflag,fieldname,fieldvec,subtitle,fmodifier);
[fignum] = UtilExperimentVectorPlot( fignum, basicvec, advancedvec, legendvec, matvec, subtitle,fmodifier);

%' for fun, try varying tau above by trying with and without the "regret" penalty, which can be done'
%' by setting "RegretPenalty" above to 0 to have no regret in objective function, and to 1.0 to have
%' a penalty in the regret. 
%' interestingly, the upper boundary gets "higher" when the regret is included: that is, it is important'
%' to sample for longer if there is an additional penalty for regret over a potentially incorrect decision'

%%%%%%%%% TESTS DIFFERENT VALUES OF GRID SIZE %%%%%%%%%%%%%%%
if ~exist('fignum','var'), fignum = 20; end;
[basic, advanced] = SetBaseParameters();
advanced.NumGridsForContours = 40; advanced.MinGridPerStdev = advanced.NumGridsForContours; % these are smaller than usual as we iterate multiple graphs
[basic, advanced, rval, msgs] = DelayInputValidator( basic, advanced );
if ~rval, msgs, end;
fieldname = 'MinGridPerStdev';       % pick the name of the field whose values are to be changed
basicflag = false;           % set to true if the field is in the 'basic' structure, false if it is in the 'advanced' structure
fieldvec = 15*2.^(0:3);    % create an array of values to rotate through for that field 
subtitle = 'Vary grid size';    % give a short subtitle name for the figures
fmodifier = fieldname;      % used to diferentiate file name if it is used for saving plots
[fignum, basicvec, advancedvec, legendvec, matvec] = TestDelayIterate(fignum,basic,advanced,basicflag,fieldname,fieldvec,subtitle,fmodifier);
[fignum] = UtilExperimentVectorPlot( fignum, basicvec, advancedvec, legendvec, matvec, subtitle,fmodifier);

%%%%%%%%% TESTS KNOWN VERSUS UNKNOWN VARIANCE %%%%%%%%%%%%%%%
if ~exist('fignum','var'), fignum = 20; end;
[basic, advanced] = SetStents();
[basic, advanced] = SetBaseParameters();
advanced.MinGridPerStdev = 20; % these are smaller than usual as we iterate multiple graphs
advanced.NumGridsForContours = 2 * advanced.MinGridPerStdev; 
basic.theta=1.0;    % enforce no discounting
basic.mumin = 1.5*basic.mumin;
basic.mumax = 1.5*basic.mumax;
[basic, advanced, rval, msgs] = DelayInputValidator( basic, advanced );
if ~rval, msgs, end;
fieldname = 'UnkVariance';       % pick the name of the field whose values are to be changed
basicflag = false;           % set to true if the field is in the 'basic' structure, false if it is in the 'advanced' structure
fieldvec = [true, false];    % create an array of values to rotate through for that field 
subtitle = 'Unknown and Known Variance';    % give a short subtitle name for the figures
fmodifier = fieldname;      % used to diferentiate file name if it is used for saving plots
[fignum, basicvec, advancedvec, legendvec, matvec] = TestDelayIterate(fignum,basic,advanced,basicflag,fieldname,fieldvec,subtitle,fmodifier);
[fignum] = UtilExperimentVectorPlot( fignum, basicvec, advancedvec, legendvec, matvec, subtitle,fmodifier);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % SECTION 3: More complicated examples of how to put several
            % graphs for several analysis as one parameter of the model is
            % varied.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%% TEST ONLINE VERSUS OFFLINE LEARNING %%%%%%%%
if ~exist('fignum','var'), fignum = 20; end;
[basic, advanced] = SetHipArthro();
[basic, advanced, rval, msgs] = DelayInputValidator( basic, advanced );
if ~rval, msgs, end;
fieldname = 'online';       % pick the name of the field whose values are to be changed
basicflag = true;           % set to true if the field is in the 'basic' structure, false if it is in the 'advanced' structure
fieldvec = [true false];    % create an array of values to rotate through for that field 
subtitle = 'Online and offline';    % give a short subtitle name for the figures
if advanced.verbose % the 'true' part of this keeps the various variables around. the 'else' part has it implemented in a subroutine, and thus does not keep the output for futher use
    % create vectors of basic and advanced structures    
    [basicvec, advancedvec, legendvec, veclen] = UtilExperimentVectorCreate( basic, advanced, basicflag, fieldname, fieldvec);
    % run the analysis for each of the structures
    for i=1:veclen
        [~, mat] = DelayCurvesRecur(basicvec(i), advancedvec(i));        % Do stage II of DP
        [mat] = DelayStageOne(basicvec(i), advancedvec(i), mat );           % Do stage I of DP
        [ fignum, mat ] = DelaySimOverview( fignum, basic, advanced, mat );
        fignum = DelaySimOutput( fignum, basic, advanced, mat );     % requires DelaySimOverview to have been called
        if i==1 matvec = mat; else matvec = [matvec mat]; end
    end
    % do the plots based on the structures
    fmodifier = 'OnOffVector';      % used to diferentiate file name if it is used for saving plots
    [fignum] = UtilExperimentVectorPlot( fignum, basicvec, advancedvec, legendvec, matvec, subtitle, fmodifier);
    for i=1:veclen
        advancedvec(i).subtitle = legendvec{i};
        advancedvec(i).filestring = [advancedvec(i).filestring num2str(i)];
        fignum = DelayDoContours(fignum, basicvec(i), advancedvec(i), matvec(i));        % Generate a bunch of contour plots
    end
else
    fmodifier = 'OnOf';      % used to diferentiate file name if it is used for saving plots
    [fignum, basicvec, advancedvec, legendvec, matvec] = TestDelayIterate(fignum,basic,advanced,basicflag,fieldname,fieldvec,subtitle,fmodifier);
    [fignum] = UtilExperimentVectorPlot( fignum, basicvec, advancedvec, legendvec, matvec, subtitle,fmodifier);
end


fignum = DelayDoContours(fignum, basic, advanced, mat);        % Generate a bunch of contour plots
fignum = UtilExperimentVectorPlot( fignum, basic, advanced, legend, mat, '', 'tst' )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPARE discounting versus no discounting  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('fignum','var'), fignum = 20; end;
[basic, advanced] = SetHipArthro();
[basic, advanced, rval, msgs] = DelayInputValidator( basic, advanced );
if ~rval, msgs, end;
basic.online=false;
advanced.saveplot = true;
fieldname = 'theta';       % pick the name of the field whose values are to be changed
basicflag = true;           % set to true if the field is in the 'basic' structure, false if it is in the 'advanced' structure
fieldvec = [basic.theta 1.0];    % create an array of values to rotate through for that field 
subtitle = 'Discounted and undiscounted reward (offline)';    % give a short subtitle name for the figures
fmodifier = fieldname;      % used to diferentiate file name if it is used for saving plots
[fignum, basicvec, advancedvec, legendvec, matvec] = TestDelayIterate(fignum,basic,advanced,basicflag,fieldname,fieldvec,subtitle,fmodifier);
[fignum] = UtilExperimentVectorPlot( fignum, basicvec, advancedvec, legendvec, matvec, subtitle,fmodifier);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPARE different sampling costs  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('fignum','var'), fignum = 20; end;
[basic, advanced] = SetHipArthro();
advanced.NumGridsForContours = 100;
[basic, advanced, rval, msgs] = DelayInputValidator( basic, advanced );
if ~rval, msgs, end;
basic.online=false;
fieldname = 'c';       % pick the name of the field whose values are to be changed
basicflag = true;           % set to true if the field is in the 'basic' structure, false if it is in the 'advanced' structure
fieldvec = [0 100*2.^(0:2)];    % create an array of values to rotate through for that field 
subtitle = 'Discounted and undiscounted reward (offline)';    % give a short subtitle name for the figures
fmodifier = fieldname;      % used to diferentiate file name if it is used for saving plots
[fignum, basicvec, advancedvec, legendvec, matvec] = TestDelayIterate(fignum,basic,advanced,basicflag,fieldname,fieldvec,subtitle,fmodifier);
[fignum] = UtilExperimentVectorPlot( fignum, basicvec, advancedvec, legendvec, matvec, subtitle,fmodifier);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPARE with and without regret  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('fignum','var'), fignum = 20; end;
[basic, advanced] = SetHipArthro();
[basic, advanced] = SetBaseParameters();
advanced.MinGridPerStdev = 60; advanced.NumGridsForContours = 60; 
advanced.DOPLOT = false;
[basic, advanced, rval, msgs] = DelayInputValidator( basic, advanced );
if ~rval, msgs, end;
%advanced.RegretPenalty=1;
basic.online=true;
fieldname = 'RegretPenalty';       % pick the name of the field whose values are to be changed
basicflag = false;           % set to true if the field is in the 'basic' structure, false if it is in the 'advanced' structure
fieldvec = [0 1 ];    % create an array of values to rotate through for that field ; negative gives a faster result, but approximati
subtitle = 'Regret Penalty (online)';    % give a short subtitle name for the figures
fmodifier = fieldname;      % used to diferentiate file name if it is used for saving plots
[fignum, basicvec, advancedvec, legendvec, matvec] = TestDelayIterate(fignum,basic,advanced,basicflag,fieldname,fieldvec,subtitle,fmodifier);
[fignum] = UtilExperimentVectorPlot( fignum, basicvec, advancedvec, legendvec, matvec, subtitle,fmodifier);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPARE different total population sizes for contract  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('fignum','var'), fignum = 20; end;
[basic, advanced] = SetHipArthro();
advanced.NumGridsForContours = 120;
[basic, advanced, rval, msgs] = DelayInputValidator( basic, advanced );
if ~rval, msgs, end;
basic.online=false; 
advanced.saveplot=true;
advanced.filestring = sprintf('on%s',num2str(basic.online));       % pick the name of the field whose values are to be changed
basicflag = true;           % set to true if the field is in the 'basic' structure, false if it is in the 'advanced' structure
fieldvec = 10.^(3:1/5:4);    % create an array of values to rotate through for that field 
subtitle = sprintf('Population size for decision (online = %s)',num2str(basic.online));    % give a short subtitle name for the figures
fmodifier = fieldname;      % used to diferentiate file name if it is used for saving plots
[fignum, basicvec, advancedvec, legendvec, matvec] = TestDelayIterate(fignum,basic,advanced,basicflag,fieldname,fieldvec,subtitle,fmodifier);
[fignum] = UtilExperimentVectorPlot( fignum, basicvec, advancedvec, legendvec, matvec, subtitle,fmodifier);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The following code tests the UtilInterSVec function. 
% It assumes that the basic, advanced structures are defined, and that the
% mat structure has been computed by running the stage II and stage I
% computations.
if ~exist('fignum','var'), fignum = 20; end;
fignum = fignum+1;
figure(fignum);
sbasic = zeros(size(advanced.simFreqDeltaVec)); %preallocate
sfix = zeros(size(advanced.simFreqDeltaVec));
for i = 1:length(advanced.simFreqDeltaVec)
    mu0 = advanced.simFreqDeltaVec(i);
    sbasic(i) = UtilInterpSVec(mat.muvec,mu0,mat.bestsvec,basic.tau,basic.TMax); % compute optimal number of stage 1 samples: This might be 0, tau or more, or something in between
    % for fix4, try to get the optimal one-shot duration
    sfix(i) = UtilInterpSVec(mat.muvec,mu0,mat.OptOneShotLength,basic.tau,basic.TMax); % compute optimal number of stage 1 samples: This might be 0, tau or more, or something in between
end
plot(advanced.simFreqDeltaVec,sbasic,'o',mat.muvec,mat.bestsvec,'x')
fignum = fignum+1;
figure(fignum);
plot(advanced.simFreqDeltaVec,sfix,'o',mat.muvec,mat.OptOneShotLength,'x')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The following code tests a KG*-type approximation for stage II when
% sampling variance is unknown.

% Setup structures
if ~exist('fignum','var'), fignum = 20; end;
fignum = fignum+1;
[basic, advanced] = SetHipArthro();
basic.c = 2.1; basic.theta = 0.999; basic.t0 = 10;
basic.mumin = 2*basic.mumin; basic.mumax=2*basic.mumax;
advanced.UnkVariance = true;
advanced.MinGridPerStdev = -05; advanced.NumGridsForContours = 50;
advanced.filestring = 'test'; advanced.dirstring = 'wow'; advanced.simNumReps = 500; 
advanced.saveplot = 1; advanced.DOPLOT = false;
% now, validate the structure
[basic, advanced, rval, msgs] = DelayInputValidator( basic, advanced );
[~, mat] = DelayCurvesRecur(basic, advanced);
[mat] = DelayStageOne(basic, advanced, mat ); 
if ~exist('fignum','var'), fignum = 20; end;

tmpstruct = mat;    % verify if copy by reference or by value
[tvec, muvec] = UtilMakeTvecMuvec(basic, advanced); % tvec should run from t0+tau to approximately t0+Tmax in small increments

basictmp=basic;
tmpstruct.muvec = muvec(:);

tvec = [advanced.dt 2.^(0:.5:5)]; % FIX: Update this so that only up to TMax - curt samples are considered
fignum = fignum + 1;
figure(fignum); hold off;
figure(fignum+1); hold off;
figure(fignum+2); hold off;
rewmat = [];
lenmat = [];
for tval = advanced.NumGridsForContours:-1:0
    curt = basic.t0 + basic.TMax * tval / advanced.NumGridsForContours;
    basictmp.t0 = curt;
    tmpstruct = DelayOptimalBayesOneStage(basictmp,advanced,tvec,tmpstruct)
    figure(fignum); 
    plot(tmpstruct.muvec,tmpstruct.OptOneShotReward);
    rewmat = [tmpstruct.OptOneShotLength(:), rewmat];
    hold on;
    figure(fignum+1);
    plot(tmpstruct.muvec,tmpstruct.OptOneShotLength);
    lenmat = [tmpstruct.OptOneShotLength(:), lenmat];
    hold on;
    figure(fignum+2);
    plot(tmpstruct.muvec,tmpstruct.OptOneShotLength > 0);
    pause;
end
figure
contour(rewmat)
figure
contour(lenmat)

[tvec, muvec] = UtilMakeTvecMuvec(basic, advanced); % tvec should run from t0+tau to approximately t0+Tmax in small increments
length(tvec)/advanced.PlotCheck
size(mat.Puppermat)
size(mat.muvec)


%%%%%%%%% TESTS PDE vs KG approach %%%%%%%%%%%%%%%
if ~exist('fignum','var'), fignum = 20; end;
[basic, advanced] = SetBaseParameters();
advanced.NumGridsForContours = 60; advanced.MinGridPerStdev = 20; % these are smaller than usual as we iterate multiple graphs
advanced.NumPointsQuadrature = 80;
advanced.saveplot = true;           % set to true to save plots in files, false to not save files automatically
advanced.dirstring = 'PDE';
basic.online = false;
[basic, advanced, rval, msgs] = DelayInputValidator( basic, advanced );
if ~rval, msgs, end;
fieldname = 'DOPDE';       % pick the name of the field whose values are to be changed
basicflag = false;           % set to true if the field is in the 'basic' structure, false if it is in the 'advanced' structure
fieldvec = [true false];    % create an array of values to rotate through for that field 
subtitle = 'PDE versus KG*-type computation';    % give a short subtitle name for the figures
fmodifier = fieldname;      % used to diferentiate file name if it is used for saving plots
[fignum, basicvec, advancedvec, legendvec, matvec] = TestDelayIterate(fignum,basic,advanced,basicflag,fieldname,fieldvec,subtitle,fmodifier);
[fignum] = UtilExperimentVectorPlot( fignum, basicvec, advancedvec, legendvec, matvec, subtitle,fmodifier);
