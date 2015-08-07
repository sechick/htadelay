function [basic, advanced] = SetBaseParameters()
% SETBASEPARAMETERS gives and example of setting the input parameters for
% the delay sequential trials project. If you are interested in modeling a
% new trial, make a copy of this file (or another Set*. file) and adapt the
% parameters. Please modify parameters after the initial call to
% DelayInputConstructor. The fields which may be changed in the basic and 
% advanced data structures can be found in the DelayInputValidator file.

% Please do not modify that file or the DelayInputValidator file (or do so 
% at your own risk). 
%
% Designed for Delay Sequential Trials paper of Chick, Forster, Pertile
% (alpha order).
%
% (c) 2014, S Chick
% Created: 14 April 2014
% Last touched: 14 April 2014
% 
% The first two batches of parameters are most critical in terms of defining the
% problem structure of the clinical trial.
[basic, advanced] = DelayInputConstructor();

%%%%%% Problem dependent parameters %%%%%%
% sampling model for Bayesian inference: normal distribution with known
% sigma and unknown mean
basic.sigma = 376.46;    % sample standard deviation
basic.mu0 = 0; %18.2;    % prior: mean for unknown mean
basic.t0 = 10;
%sigma0 = basic.sigma/sqrt(basic.t0);   % prior: std dev for unknown mean

% Patient delay and pool
basic.TMax = 2000;    %total number of patients in the population for testing
basic.c = 1;  % variable cost per patient of doing trial

delayinyears = 1/12;     % one month is 1/12.
annualdiscountrate = .045;  % real valued number, with value like 0.05 for 5% discount rate
patientsperyear = 1000; % maximum number of patients put into trial per year
basic.tau = patientsperyear * delayinyears;

basic.online = false;  % set to true if per patient reward is gained, false for no 
        % per patient reward for those in trial, but only reward is gained
        % upon completion of trial (and thus all reward is embedded in the
        % reward function, except for sampling and discounting costs
advanced.fixedP = true;  % set to true if there are exactly PPatients to follow once a treatment or alternative is selected
        % set to false if the patients which are not included in trial due
        % to early stopping are also included in the treatment arm, that
        % is, if fixedP is false, then there are P + (TMAX - T) patients
        % which are included, where T is the stopping time (number of
        % patients) included in the trial (and T \in 0, 1, ..., TMAX).
advanced.nochangeallowed = false;    % set to false (the better default) meaning that it is ok
        % to change what you think is best after the added tau data points
        % come in. set this to true if you want to ignore the added tau
        % pending data points to make a decision.

% Implementation stuff
basic.PPatients = 1*10e6;  %total number of patients to be treated following adoption
basic.ICost = 1*10e6;       % fixed cost of investment

% computed data
basic.theta = exp( - log(1+annualdiscountrate) / patientsperyear ); %  per patient discount rate

%%%%%% Numerical computation parameters %%%%%%

basic.mumax = +1200;    % maximum posterior mean per patient; pick big enough to be above upper stopping boundary
basic.mumin = -1200;    % minimum posterior mean per patient; pick small enought to be below lower stopping boundary

%advanced.MinGridPerStdev = -0.004;    % if this is negative, code tries to use its negative as the delta-t
advanced.MinGridPerStdev = 80;    % minimum number of grid points per stdev of output, or negative to use dt=1

% For plotting
advanced.smoothed = true;   % true or false if curve smoothing is desired.
advanced.DOPLOT = false; % put true or false depending if you do or do not want dynamic plots
advanced.NumGridsForContours = 80; % number of time values for generating contour plots (also num plots in animation to use during recursion over t)

advanced.bigfontsize=16;
advanced.smallfontsize=14;
advanced.fontname='Helvetica';

advanced.titlestring='AZ-ER versus AC';    % for use in figures / display
advanced.filestring='JRSSA';     % for use in saving files - should be a string with only characters, no blanks etc.
advanced.dirstring='JRSSA';     % for use in saving files - should be a string with only characters, no blanks etc.

end