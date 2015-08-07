function [ basicout, advancedout, rval, errorlist ] = DelayInputValidator( basic, advanced )
%DELAYINPUTVALIDATOR checks the basic and advanced data structures for
%validity, and corrects values if needed, for usage in the Delay Differential 
% equations for the delay paper of Chick, Forster, Pertile (in alpha order). 
%It returns two main parameters, one with 'basic' parameters, and which 
%should be used instead of the basic and advanced parameters which are
%input to this routine. That is, use the basic and advanced OUT parameters
%for passing into DelayCurvesRecur (stage II of the backwards recursion, 
%which is a backwards recursion for a continuous time stopping problem),
%and then DelayStageOne (stage I of the backwards recursion, which is a
%one-stage lookahead story).
%
% INPUTS:
%   basic, advanced: two structures with parameters for the sequential
%       sampling problem and its computation. Typically this should be created
%       by calling the constructor method, DelayInputConstructor.
% OUTPUTS:
%   basicout, advancedout: updated copies of the structures which have been
%       passed as inputs, but with values 'tweaked' in order to make them
%       amenable for the computations. Some values are complicated
%       functions of the other parameters, and so it is important to use
%       the 'out' version of the parameters when calling the recursions.
%   rval: true if all is ok, false if there were some errors worth noting.
%   errorlist: text string which contains messages which explains some of
%       the errors worth noting, and how DelayInputValidator went about
%       fixing them.
%
% NOTE: This routine actually creates new fields in the data structures,
%   including: NumT, PlotCheck, muvec
%
% Designed for Delay Sequential Trials paper of Chick, Forster, Pertile
% (alpha order).
%
% (c) 2014, S Chick, P Pertile, M Forster
% Created: 14 April 2014
% Last touched: 14 April 2014
% 
% The first two batches of parameters are most critical in terms of defining the
% problem structure of the clinical trial.
routinename='DelayInputValidator';
errorlist = '';

MINBIGX = 5;        % mimimum number of dmu per sigma for numerical analysis

[defaultbasic, defaultadvanced] = DelayInputConstructor();  % find the default values used in the construtor method.

%currentDate = datestr(now,'mmmdd');
%myStruct.(currentDate) = [1,2,3]
rval = 1;

if basic.sigma <= 0
    basic.sigma = defaultbasic.sigma;  % standard deviation of sample differences (between people in two trial arms)
    errorlist = sprintf('%s\n%s: sampled std dev was negative: sigma reset to %f\n',errorlist,routinename,basic.sigma);
    rval = 0;
end
if basic.c < 0
    basic.c = defaultbasic.c;  % cost per sample
    errorlist = sprintf('%s\n%s: negative sampling cost, c reset to %f\n',errorlist,routinename,basic.c);
    rval = 0;
end
if (basic.theta > 1) || (basic.theta <= 0.0)
    tmpval = basic.theta;
    basic.theta = defaultbasic.theta;  % discount factor
    errorlist = sprintf('%s\n%s: invalid discount factor (%f), theta reset to %f\n',errorlist,routinename,tmpval,basic.theta);
    rval = 0;
end
if basic.PPatients <= 0
    tmpval = basic.PPatients;
    basic.PPatients = defaultbasic.PPatients; % number of patients covered post-trial
    errorlist = sprintf('%s\n%s: PPatients should be positive, reset to %f\n',errorlist,routinename,tmpval);
    rval = 0;
end
if basic.ICost < 0
    basic.ICost = defaultbasic.ICost;  % fixed cost if adoption of new technology is made
    errorlist = sprintf('%s\n%s: invalid fixed cost ICost, ICost reset to %f\n',errorlist,routinename,basic.ICost);
    rval = 0;
end
if basic.tau < 0
    basic.tau = defaultbasic.tau;  % delay (in number of patients treated) until data is observed
    errorlist = sprintf('%s\n%s: tau must be nonnegative, and was reset to %f\n',errorlist,routinename, basic.tau);
    rval = 0;
end
if basic.TMax <= basic.tau
    basic.TMax = 2 * basic.tau;  % maximum number of samples (patient pairs) to be allowed during trial
    errorlist = sprintf('%s\n%s: TMax too small relative to tau, and was reset to %f\n',errorlist,routinename,basic.TMax);
    rval = 0;
end

%basic.online = true;    % default: online learning, meaning results of patients in trial are counted in expected reward

if basic.t0 <= 1.0
    basic.t0 = defaultbasic.t0;  % effective number of samples in prior distribution for unknown mean
    errorlist = sprintf('%s\n%s: t0 should exceed 1.0, reset to %f\n',errorlist,routinename,basic.t0);
    rval = 0;
end
if (basic.mu0 > basic.mumax) || (basic.mu0 < basic.mumin)
    basic.mu0 = (basic.mumax + basic.mumin) / 2;
    errorlist = sprintf('%s\n%s: warning, potentially invalid value of mu0, reset to %f\n',errorlist,routinename,basic.mu0);
end

if advanced.UnkVariance     % if variance is unknown, then validate scale parameter for unknown variance
    if advanced.UnkVarianceShape == -1.0  % if we don't use the three parameter form for the unknown variance, 
        if basic.t0 < 4     % xi_0 <- (t0-1)/2, and var unknown mean has 2xi_0 = (t0-1) dof, and want var to exist
            rval = 0;
            errorlist = sprintf('%s\n%s: error, basic.t0 is %f but should be at least 4 if variance is unknown\n',errorlist,routinename,basic.t0);
        end
    elseif advanced.UnkVarianceShape < 2
        rval = 0;
        errorlist = sprintf('%s\n%s: error, advanced.UnkVarianceShape is %f but should be at least 2 if variance is unknown\n',errorlist,routinename,advanced.UnkVarianceShape);        
    end
    % NOTE: if variance is unknown, 
    if advanced.DOPDE
        errorlist = sprintf('%s\n%s: warning, if sampling variance is UNKnown, normally one should set advanced.DOPDE to false\n',errorlist,routinename,advanced.UnkVarianceShape);        
    end
%    warning('%s\n%s: error, advanced.UnkVarianceShape is not supported for dynamic program solution, but is for simulation tests in sense of using plug-in estimators for variance\n',errorlist,routinename);        
else
    if ~advanced.DOPDE
        errorlist = sprintf('%s\n%s: warning, if sampling variance is Known, normally one should set advanced.DOPDE to true\n',errorlist,routinename,advanced.UnkVarianceShape);        
    end
end

if length(advanced.DistributionType) > 0    % if a sampling distribution type is given, check to see if it is implemented in current directory, and check to see if the distribution instance is of the correct object type
% Following commented code identifies whether a given object is implemented
% or not in the 'current directory'.  removed, so that it does not require
% the sampling distribution object to be implemented in the current working
% directory - can be elsewhere on path
%    dirData = dir();      %# Get the data for the current directory
%    dirIndex = [dirData.isdir];  %# Find the index for directories
%    fileList = {dirData(dirIndex).name};  % Get list of directories
%    definedRVs={fileList{strncmp(fileList, '@Dist', 4)}};  % set of RVs available should include those in directories starting with '@Dist'
%    if ~sum(strcmp(definedRVs, strcat('@',func2str(advanced.DistributionType))))
%        errorlist = sprintf('%s\n%s: warning, advanced.DistributionType %s not found in current directory\n',errorlist,routinename,strcat('@',func2str(advanced.DistributionType)));        
%    end
    if ~strcmp(class(advanced.Distribution),func2str(advanced.DistributionType))
        errorlist = sprintf('%s\n%s: warning, advanced.DistributionType is %s but instance advanced.Distribution is of type %s\n',...
            errorlist,routinename,strcat('@',func2str(advanced.DistributionType)), class(advanced.Distribution));        
    end
end

if advanced.verbose == true
    newval = 20;     % in fact, to get a reasonable contour plot this should be at least 20 or 40
    if advanced.NumGridsForContours <= newval
%        errorlist = sprintf('%s\n%s: NumGridsForContours too small, was reset to %f\n',errorlist,routinename,newval)
        errorlist = sprintf('%s\n%s: NumGridsForContours might be too small, suggest at least %f\n',errorlist,routinename,newval);
%        basic.NumGridsForContours = newval;  % standard deviation of sample differences (between people in two trial arms)
    end
end

if advanced.StageOneChecks > 0 % for optimal stage 1 sampling budget, set to 0 or less to check at all integers, use positive
                    % number, such as 50, to specify that one should check at s = tau/50, 2 tau/50, 3 tau/50, ..., tau
    if ceil(advanced.StageOneChecks) ~= floor(advanced.StageOneChecks)  % insure it is an integer if it is positive
        advanced.StageOneChecks = -1;
        errorlist = sprintf('%s\n%s: StageOneChecks should be negative (to check all integer stage 1 sample sizes) or a positive integer number of values to check, was reset to %f\n',errorlist,routinename,advanced.StageOneChecks);
    end
end

%    advanced.fixedP = true; % set to true if expected reward on stopping is for P * expected reward per patient, false if patients not tested due to early stopping can also benefit from better alterantive
%    advanced.nochangeallowed = false;   % default to false, so that one waits til all data in (tau samples more), before selecting. If true, then choice must be made before tau outstanding samples are observed
%    advanced.verbose = true; % print out summary information during the run
%    advanced.DOPLOT = false; % print out plots which are useful for debugging or seeing the dynamics of the backward recursions
%    advanced.RegretPenalty = 0;   % 0 for no regret penalty, positive for regret penalty
%    advanced.smoothed = true;   % smoothen the stopping boundaries if true, otherwise do not smoothen them.
if (advanced.MAXPU <= .1) || (advanced.MAXPU > 0.49)
    advanced.MAXPU = defaultadvanced.MAXPU;
    errorlist = sprintf('%s\n%s: MAXPU was reset to %f\n',errorlist,routinename,advanced.MAXPU);
end

if isnan(advanced.RegretPenalty)
   advanced.RegretPenalty = defaultadvanced.RegretPenalty;
   errorlist = sprintf('%s\n%s: NaN no longer supported for RegretPenalty: reset to %f\n',errorlist,routinename, advanced.RegretPenalty);
end
if advanced.RegretPenalty < 0
   advanced.RegretPenalty = max(0,defaultadvanced.RegretPenalty);
   errorlist = sprintf('%s\n%s: RegretPenalty needs to be nonnegative: reset to %f\n',errorlist,routinename, advanced.RegretPenalty);
end

%    advanced.titlestring = 'Title of trial for purposes of plot.';  % name used in plots
%    advanced.filestring = 'fName';      % root/base of name of files if any files are to be saved. make this an empty string if no files are to be saved
%    advanced.fontname = 'Helvetica';    % name of font to use for plots
%    advanced.bigfontsize = 16;          % size of font for larger items in plot
%    advanced.smallfontsize = 14;        % size of font for smaller items in plot
if (advanced.fracheight < 0.1) || (advanced.fracheight >= 1.0)
    advanced.fracheight = defaultadvanced.fracheight;  % if 0, then pure expected reward. if >0 then reward is penalized by expected regret of stopping at that time
    errorlist = sprintf('%s\n%s: fracheight should be fraction of screen size for plot, was reset to 0.8\n',errorlist,routinename,advanced.fracheight);
end

    % Monte Carlo data
%    advanced.simNumReps = 400;           % 0 for no simulations, >2 for running sample paths of the trials, and simulation estimates of power curves
%    advanced.CRN = true;                % true to use CRN for noise across Bayes and frequentist estimations, across MAT/experiments
%    advanced.CRNAcrossExperiment = 0;   % seed to use if CRN is desired across the MAT experiment
%    advanced.CRNAcrossBayesMu = true;   % true to use CRN acruss mu values for a given MAT in TestDelayIterate experiments 
%    advanced.PLOTSIMS = true;
%    advanced.keepAllOutput = false;      % true to keep all output from all replications, false if only means and std of various items are to be kept
%    advanced.PDEWithSimPlots = true;     % true if simulation plots should be accompanied by PDE results where possible
%    advanced.NumPointsQuadrature = 81; % SHOULD BE ODD INTEGER: number of points for quadrature for computing expected regret once tau samples come in
if advanced.simNumReps < 0
    advanced.simNumReps = 0;
    errorlist = sprintf('%s\n%s: simnumreps must be at least 0 (0 to run no reps)\n',errorlist,routinename);
end
if advanced.NumPointsQuadrature < 1
    advanced.NumPointsQuadrature = defaultadvanced.NumPointsQuadrature;
    errorlist = sprintf('%s\n%s: NumPointsQuadrature must be at least 1 (reset to %f)\n',errorlist,routinename),advanced.NumPointsQuadrature;
end

%advanced.z = 0.0; % test statsitic for frequentist tests of whether new treatment exceeds the old. 
% should default to 0 if only means are to be used, set to 1.96 for example, if stronger evidence of new better than old is needed

% next few lines figure out how big dt and dw should be, to have pu = pd =
% MAXPU when we get back to start of time horizon.  
% For this, see page 2 of notes from Steve on 22 Mar 09 from Arlotto Chick
% Gans project on leraning curves.
MINTzero = basic.t0;
if advanced.UnkVariance
    if advanced.UnkVarianceShape == -1
        tmpxi = (basic.t0) / 2;         % use if there is a noninformative prior for mu, sigma^2 of sampling distribution
%        tmpxi = (basic.t0 - 1) / 2;         % use if there is a noninformative prior for mu, sigma^2 of sampling distribution
    else
        tmpxi = (advanced.UnkVarianceShape - 1) / 2;    % use if the prior is otherwise specified
    end
	fudgefactor = sqrt( (2 * tmpxi) / (2 * tmpxi - 2) );
else
    fudgefactor = 1.0;
end
if advanced.MinGridPerStdev >= 1
    bigX = max(advanced.MinGridPerStdev, sqrt(1.01*2*advanced.MAXPU*MINTzero) );   % how many dmu per sigma do we need?
    advanced.dmu = basic.sigma / bigX ;
    advanced.dt = ( 2*advanced.MAXPU*MINTzero^2/bigX^2/fudgefactor ) / ( 1 - (2*advanced.MAXPU*MINTzero/bigX^2/fudgefactor) );
elseif advanced.MinGridPerStdev < 0    % try to set dt to the negative of MinGridPerStdev
    advanced.dt=-advanced.MinGridPerStdev;
    bigX = sqrt( 2*advanced.MAXPU*MINTzero^2 / advanced.dt );   % how many dmu per sigma do we need?
%    advanced.dmu = basic.sigma / bigX ;
    advanced.dmu = fudgefactor * basic.sigma * sqrt( advanced.dt / (2 * advanced.MAXPU * MINTzero * (MINTzero + advanced.dt)) );
    bigX = basic.sigma / advanced.dmu;
else    % try dt = 1
    advanced.dt=1;
%    bigX = sqrt(2*advanced.MAXPU*MINTzero^2 / advanced.dt);   % how many dmu per sigma do we need?
%    advanced.dmu = basic.sigma / bigX ;
    advanced.dmu = fudgefactor * sqrt( basic.sigma^2 * advanced.dt / (2 * advanced.MAXPU * MINTzero * (MINTzero + advanced.dt)) );
    bigX = basic.sigma / advanced.dmu;
end

if bigX < MINBIGX
    bigX=MINBIGX;
    advanced.dmu = basic.sigma / bigX ;
    advanced.dt = ( 2*advanced.MAXPU*MINTzero^2/bigX^2/fudgefactor ) / ( 1 - (2*advanced.MAXPU*MINTzero/bigX^2/fudgefactor) );
    errorlist = sprintf('%s\n%s: might need smaller time step or more delta-mu per standard deviation for accuracy, dt reset to %f\n',errorlist,routinename,advanced.dt);
end

% These are some new parameters which might not have been set up in by
% DelayInputConstructor. Indeed, they are best set by DelayInputValidator
% rather than by the constructor, as they may depend in a complicated way
% on the other parameters. Be sure to call UtilMakeTvecMuvec in conjuction
% with these issues.
% Places where UnkVariance is handled differently than known variance
% UtilMakeTvecMuvec()
% DelayCurvesRecur()
% DelaySimComputer()
if advanced.UnkVariance
%FIX: Can set dt to 1 for the unknown variance case, but then the PCS and expected number of samples
% might not be computed correctely. Would need fixing / numerical stability checks in DelayCurvesRecur
%   advanced.dt = 1;   % SEC: Note, this is due to KG* formulation rather than PDE for unknown variance case.
end
tlen = 1+max(1,ceil( (basic.TMax - basic.tau) / advanced.dt ));
advanced.PlotCheck = max(2,floor(tlen/advanced.NumGridsForContours));  % do plot every 'PlotCheck' steps in time

% find grid of mu values over full hi-low range, and with grid point at mu0
% this grid allows us to return a value function estimation at precisely
% the value of mu0 without having to interpolate 
advanced.mushift=mod(basic.mu0,advanced.dmu);
%advanced.mushift=mod(basic.ICost / basic.PPatients,advanced.dmu);
basic.mumax = advanced.dmu * ceil( (advanced.mushift+basic.mumax) / advanced.dmu );
basic.mumin = advanced.dmu * floor( (-advanced.mushift+basic.mumin) / advanced.dmu );

% The following code reconstructs the time grid and mu grid in full. Better
% to call it where it is needed: UtilMakeTvecMuvec
%advanced.tvec= (basic.t0 + basic.tau) + advanced.dt * (0:(advanced.tlen-1));
%advanced.muvec = advanced.mushift+(basic.mumin:advanced.dmu:basic.mumax) ;        % handle integer number of grid points, cover range of mu needed

if ~((advanced.simNumReps == 0) || (advanced.simNumReps > 2))
    advanced.simNumReps = defaultadvanced.simNumReps;  % effective number of samples in prior distribution for unknown mean
    errorlist = sprintf('%s\n%s: simNumReps should be 0 or greater than 2, reset to %f\n',errorlist,routinename, advanced.simNumReps);
    rval = 0;
end

if isfield(advanced,'simFreqDeltaVec')
    if advanced.simNumReps == 0
        advanced.simNumReps = defaultadvanced.simNumReps;  % effective number of samples in prior distribution for unknown mean
        errorlist = sprintf('%s\n%s: simNumReps should begreater than 2 if simFreqDeltaVec non-empty, simNumReps reset to %f\n',errorlist,routinename,advanced.simPowerSteps);
        rval = 0;
    end
else
    numinsimFreqDeltaVec=200;  % can change the number of points in freqdeltavec. the 11/3 in the next line is to get \pm 11/3 standard errors in width for vector's range
    advanced.simFreqDeltaVec = (basic.sigma/sqrt(basic.t0)) * (11/3) * (-ceil(numinsimFreqDeltaVec/2):ceil(numinsimFreqDeltaVec/2)) / ceil(numinsimFreqDeltaVec/2);
end
        
%rval = isempty(errorlist);
%if advanced.verbose && (~rval)
if ~isempty(errorlist)
    warning(errorlist)
end
if rval
    basicout = basic;
    advancedout = advanced;
else
    basicout = [];
    advancedout = [];
end

end

