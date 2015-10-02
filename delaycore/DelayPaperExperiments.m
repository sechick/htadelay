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
% (*) Generate expected value of information for various routines, for both
% single runs, and for TestDelayIterate experiments.
% (*) Test simulation results: why do some differ from the PDE results?
% (*) Learn more about nonmonotonic upper bound with large tau and
% discounting (in some cases). DelayCOntinExpectation calls? TerminalReward
% and Terminal regret test functions to visualize? Testing combinations of
% online / offline, fixedP, etc. Other debugging as arises.
% (*) Do functions for certain frequentist statistics (be it PDE or be it
% via simulation).
% (*) Code for various other trials we might be interested in looking at.
% (*) Do parking lot stuff to improve API/GUI after we know better what is 
% functionality we can do, and what end user usage cases might look like.
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

% An example with 'stent' parameters
if ~exist('fignum','var'), fignum = 20; end;
[basic, advanced] = SetStents();
advanced.DOPLOT = 0;
[basic, advanced, rval, msgs] = DelayInputValidator( basic, advanced );
if ~rval, msgs, end;
[~, mat] = DelayCurvesRecur(basic, advanced);
[mat] = DelayStageOne(basic, advanced, mat ); 
[ fignum, mat ] = DelaySimOverview( fignum, basic, advanced, mat );
fignum = DelaySimOutput( fignum, basic, advanced, mat );     % requires DelaySimOverview to have been called

fignum = DelayDoContours(fignum, basic, advanced, mat);        % Generate a bunch of contour plots
fignum = UtilExperimentVectorPlot( fignum, basic, advanced, legend, mat, '', 'tst' )

fignum = DelayDoThreeTestPlots(fignum, basic, advanced);

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
advanced.NumGridsForContours = 50; advanced.MinGridPerStdev = 30; % these are smaller than usual as we iterate multiple graphs
advanced.NumPointsQuadrature = 100;
advanced.saveplot = true;           % set to true to save plots in files, false to not save files automatically
advanced.dirstring = 'tau';
basic.online = false;
[basic, advanced, rval, msgs] = DelayInputValidator( basic, advanced );
if ~rval, msgs, end;
fieldname = 'tau';       % pick the name of the field whose values are to be changed
basicflag = true;           % set to true if the field is in the 'basic' structure, false if it is in the 'advanced' structure
%fieldvec = [0 20 40 80 250 500 1000];    % create an array of values to rotate through for that field 
fieldvec = [0 20 40 250 1000];    % create an array of values to rotate through for that field 
subtitle = 'Vary delay \tau';    % give a short subtitle name for the figures
fmodifier = fieldname;      % used to diferentiate file name if it is used for saving plots
[fignum, basicvec, advancedvec, legendvec, matvec] = TestDelayIterate(fignum,basic,advanced,basicflag,fieldname,fieldvec,subtitle,fmodifier);
[fignum] = UtilExperimentVectorPlot( fignum, basicvec, advancedvec, legendvec, matvec, subtitle,fmodifier);

%' for fun, try varying tau above by trying with and without the "regret" penalty, which can be done'
%' by setting "RegretPenalty" above to 0 to have no regret in objective function, and to 1.0 to have
%' a penalty in the regret. 
%' interestingly, the upper boundary gets "higher" when the regret is included: that is, it is important'
%' to sample for longer if there is an additional penalty for regret over a potentially incorrect decision'

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
