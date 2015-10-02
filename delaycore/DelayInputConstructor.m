function [ basic, advanced ] = DelayInputConstructor( basicarray, advancedarray  )
%DELAYINPUTCONSTRUCTOR creates two data structures which are used as inputs
%to the Delay Differential equations for the delay paper of Chick, Forster,
%Pertile (in alpha order). It returns two main parameters, one with 'basic'
%parameters, each of which should be reset/checked for adapting the code to
%the specific problem at hand (which assumes that two arms are involved in
%the trial, and that the trial may have delays in results coming, the
%results are real-valued with normal distributions (or at least
%approximately so), the prior distribution of the mean is assumed to be a
%normal distribution, there is assumed to be a maximum number of patients
%which can be tested following which a specified number of patients are put
%on the new therapy.
%
% The 'advanced' parameter structure has additional numerical computation
% parameters, parameters for displaying and/or saving figures, etc.
%
% The default fields for the arrays will be filled in. There is no need to
% call the Constructor with basicarray or advancedarry. However, if those
% parameters are passed, they are assumed to be a cell array of even
% length, with the names of the parameters and the default values assigned
% to it.
%       % basicarray = { 'c', 2.1, 'theta', 0.999 }, for example
%
%DELAYINPUTVALIDATOR can be used to validate parameter choices which are
%passed to it, to check values which may have been set by the user, and in
%any event SHOULD BE CALLED anyway after this routine is used to construct
%the input parameters and they are tweaked by the end user, in order to
%complete the dmu, dt and other values which are useful in the main
%computations (NumT, tvec, tlen, PlotCheck, muvec)
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
%
% The 'usual' way to set advanced.MinGridPerStdev is to choose it to be the
% number of 'dmu' steps per standard deviation (e.g. 20 for course grid, 40
% or bigger for a finer grid). Bigger means higher precision (but slower
% run times). 
%
% NOTE: if variance is unknown, it is best to use a larger dt (e.g. set
% reduce advanced.MinGridPerStdev, OR set advanced.MinGridPerStdev to 
% -(the desired dt), for example set it to -0.5 to get a dt of 0.5 and to 
% set dmu appropriately).
    if nargin < 2
        advancedarray = {};
    end
    if nargin < 1
        basicarray = {};
    end

    basic.c = 0;            % marginal cost per sampled pair
    basic.theta = 0.999;    % discount factor per sampled pair. Typically user should set at least one of c > 0 or theta somewhat less than 1. otherwise, there is no cost of cotinuing to sample (for the offline learning case)
    basic.PPatients = 10000;% number of patients to be covered by a contract if an adoption or rejection decision is made
    basic.ICost = 0.0;      % fixed investment cost.
    basic.tau = 40;         % default delay of test results is 0. Otherwise, is number of patients which can be treated until data comes in
    basic.TMax = 1000;      % maximum number of patient pairs allowed in the sequential variation of the trial
    basic.numpairs = 700;   % number of patient pairs to be tested for simulations of the 'as puplished' trial (e.g fixed sample size)
    basic.online = true;    % default: online learning, meaning results of patients in trial are counted in expected reward

% If there is a sampling distribution which is normally distributed, and
% which has unknown mean and known variance, and the PDE calculations can
% rely on knowing the variance, then use the following:
    basic.mu0 = 0;          % mean of prior distribution for uknown mean reward
    basic.t0 = 5;           % effective number of samples in prior distribution for unknown mean reward
    basic.sigma = 100;      % standard deviation of sample differences (between people in two trial arms)
    advanced.UnkVariance = false;   % false if variance is known, true if unknown (default to known variance)
    advanced.DistributionType = []; % set to empty string unless object used for sampling distribution
    advanced.Distribution = []; % 

% If there is a sampling distribution which is normally distributed, and
% which has unknown mean and known variance, and the PDE calculations can
% NOT rely on knowing the variance, then use the following:
%    basic.mu0 = 0;          % mean of prior distribution for uknown mean reward
%    basic.t0 = 5;           % effective number of samples in prior distribution for unknown mean reward
%    basic.sigma = 100;      % standard deviation of sample differences (between people in two trial arms)
%    advanced.UnkVariance = true;   % identify if the stopping boundary will be adjusted to account for sampling variance or not
%    advanced.DistributionType = @DistNormalMu; % identify the type of sampling distribution object
%    advanced.Distribution = advanced.DistributionType(basic.mu0, basic.t0, basic.sigma); % create an instance of the object with the right hyperparameters

% If there is a sampling distribution which is normally distributed, and
% which has unknown mean and UNknown variance, and the PDE calculations can
% NOT rely on knowing the variance, then use the following:
%    basic.mu0 = 0;          % mean of prior distribution for uknown mean reward
%    basic.t0 = 5;           % effective number of samples in prior distribution for unknown mean reward
%    basic.sigma = 100;      % standard deviation of sample differences (between people in two trial arms)
%    xi0 = 20;      % shape parameter for unknown variance
%    advanced.UnkVariance = true;   % identify if the stopping boundary will be adjusted to account for sampling variance or not
%    advanced.DistributionType = @DistNormalMuSig; % identify the type of sampling distribution object
%    advanced.Distribution = advanced.DistributionType(basic.mu0, basic.t0, basic.sigma, xi0); % create an instance of the object with the right hyperparameters
    
    advanced.UnkVarianceShape = -1.0; % Ignored for the moment. Might later be used for fudge factor for boundary
        % Set UnkVarianceShape = -1 to have a 3-parameter conjugate prior
        %   for the unknown mean and variance, as in Chick and Inoue (Oper Res 2001), eq. (2) and (3)
        % Set UnkVarianceShape = shape parameter of inverted gamma distribution for
        %   unknown variance, as in \xi_{i,0} of Chick and Frazier (Man Sci 2012), eq. (20).
        %   In this case, \chi_{i,0} is set so that basic.sigma =
        %   \chi_{i,0} / (\xi_{i,0} - 1), meaning that basic.sigma gives the
        %   a priori mean value of the unknown variance

    advanced.mushift = 0;
    advanced.NumGridsForContours = 50;  % number of time points to include for time axis in contour plots
    advanced.MinGridPerStdev = 30;    % minimum number of delta-mus per standard deviation
    advanced.dt = NaN;       % delta t for the recursion: set to NaN if it should be computed from dmu
    advanced.fixedP = true; % set to true if expected reward on stopping is for P * expected reward per patient, false if patients not tested due to early stopping can also benefit from better alterantive
    advanced.nochangeallowed = false;   % default to false, so that one waits til all data in (tau samples more), before selecting. If true, then choice must be made before tau outstanding samples are observed
    advanced.verbose = true; % print out summary information during the run
    advanced.DOPLOT = true ; % print out plots which are useful for debugging or seeing the dynamics of the backward recursions
    advanced.RegretPenalty = 0; % positive for regret, or 0 for standard calculation (no regret penalty).
    advanced.smoothed = true;   % smoothen the stopping boundaries if true, otherwise do not smoothen them.
    advanced.MAXPU = 0.475; % maximum probability of going up (or down) in the trinomial tree which is to be constructed.
    advanced.StageOneChecks = -1; % for optimal stage 1 sampling budget, set to 0 or less to check at all integers, use positive
    advanced.DOPDE = true;  % use the PDE mechanism to compute stage II if true (the default), or a KG-type otherwise (DelayDriver uses a [dt, 1, 1.414, 2,...] step look ahead).
%    advanced.KGNoPCS = false; % don't compute the PCS and ENumSamps, ONLY
%    active if DOPDE is false (so that KG* computation is active)
    % number, such as 50, to specify that one should check at s = tau/50, 2 tau/50, 3 tau/50, ..., tau
    
    advanced.titlestring = 'Title of trial for purposes of plot.';  % name used in plots
    advanced.filestring = 'fName';      % root/base of name of files if any files are to be saved. make this an empty string if no files are to be saved
    advanced.dirstring = 'dirName';    % root/base of name of directory for output files (figures, data) related to this analysis
    advanced.graphicextension='eps'; % can also use 'pdf', etc. as in saveas() function for figures.
    advanced.saveplot=true ;        % set to 'true' if you want to have plots saved, the default is false, to not save plots
    advanced.fontname = 'Helvetica';    % name of font to use for plots
    advanced.bigfontsize = 16;          % size of font for larger items in plot
    advanced.smallfontsize = 14;        % size of font for smaller items in plot
    advanced.fracheight = 0.95;          % fraction of size of screen for standardizing the plot size.

    % Monte Carlo data
    advanced.simNumReps = 200 ;           % 0 for no simulations, > 2 for running sample paths of the trials, and simulation estimates of power curves
% example usage of this parameter - set in DelayInputValidator if not set otherwise.
    advanced.CRN = true;                % true to use CRN for noise across Bayes and frequentist estimations, across MAT
    advanced.CRNAcrossExperiment = 37;  % seed for use with common random numbers (if CRN is true), should be type int32.
    advanced.CRNAcrossBayesMu = true;   % true to use CRN across mu in TestDelayIterate experiments too
    advanced.PLOTSIMS = true;
    advanced.keepAllOutput = true;      % true to keep all output from all replications, false if only means and std of various items are to be kept
    advanced.PDEWithSimPlots = true;     % true if simulation plots should be accompanied by PDE results where possible
    advanced.DoRegretIntegral = false;  % TRUE if regret should be done to high accuracy (integral) or FALSE if quadrature approximation is to be used
    advanced.NumPointsQuadrature = 80; % If advanced.DoRegretIntegral is false, then INTEGER number of points for quadrature approximation for regret.

    advanced.z = 0.0;       % used as a test statistic for frequentist routines: critical value for statistical inference test (e.g. 1.96), set to 0 if only care about mean 
% Now that basic parameters are set up, check to see if values have been
% passed by end user. Note that the parameters further down which are
% computed especially by parameters above need to be handled specially.
% Same for SOME parameters which are set by the validator routine, e.g.
% allow end user to specify simFreqDeltaVec,but not tlen, PlotCheck,
% mushift as those are selected for internal consistency purposes.
    [rval, basic, advanced] = DelayInputModifier(basic, advanced, basicarray, advancedarray);
    
% COMPUTED PARAMETERS: If more are put here, please fix the code to insure
% that they don't overwrite parameter values which may have been entered by
% the end user.
if ~rval
    basic = [];
    advanced = [];
else
    if advanced.UnkVariance
        if advanced.UnkVarianceShape == -1
            tmpxi = (basic.t0) / 2;         % use if there is a noninformative prior for mu, sigma^2 of sampling distribution
%            tmpxi = (basic.t0 - 1) / 2;         % use if there is a noninformative prior for mu, sigma^2 of sampling distribution
        else
            tmpxi = (advanced.UnkVarianceShape - 1) / 2;    % use if the prior is otherwise specified
        end
        fudgefactor = sqrt( (2 * tmpxi) / (2 * tmpxi - 2) );
    else
        fudgefactor = 1.0;
    end
	if ~isfield(advanced,'dmu') % the increment used for the grid for the posterior mean
        advanced.dmu = fudgefactor * basic.sigma / advanced.MinGridPerStdev;  
	end
	if ~isfield(basic,'mumax')  % the upper value to be included in the grid for the posterior mean
        basic.mumax = 15*basic.sigma/sqrt(basic.t0);           
	end
    if ~isfield(basic,'mumin')  % the lower value to be included in the grid for the posterior mean
        basic.mumin =  -basic.mumax;       
    end
end


end

