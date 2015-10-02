function [basic, advanced] = SetTestExample()
% SETTESTEXAMPLE gives and example of setting the input parameters for
% the delay sequential trials project. 
%
% Designed for Delay Sequential Trials paper of Chick, Forster, Pertile
% (alpha order).
%
% (c) 2014, S Chick
% Created: 17 April 2014
% Last touched: 17 April 2014
% 
% The first two batches of parameters are most critical in terms of defining the
% problem structure of the clinical trial.
[basic, advanced] = DelayInputConstructor();

%%%%%% Problem dependent parameters %%%%%%
% 95% CI (-1795, 1960)
% mean: 82.46
% numpatients: 60 in one arm 66 in the other: so, lets use 60.
% delayinyears = 1;
% length of study : 3 years
% pick TMax bigger.
% discount: was not discounted, 

CIlow=-1795; CIhi = 1960;
Conflevel = 0.95;
numpairs = 60;
stder = (CIhi - CIlow) / norminv( 1 - (1-Conflevel)/2)
postmean = 82.46;
studyyears = 3;
delayinyears = 1     % one month is 1/12.

basic.sigma = stder * sqrt(numpairs);
basic.mu0=82.46;
basic.t0 = 5
%sigma0 = basic.sigma/sqrt(basic.t0);   % prior: std dev for unknown mean % no need to assign, will be computed

% Patient delay and pool
basic.TMax = numpairs;
annualdiscountrate = .03;  % real valued number, with value like 0.05 for 5% discount rate
patientsperyear = numpairs / (studyyears - delayinyears)  % maximum number of patients put into trial per year
basic.c = 10  % variable cost per patient of doing trial
basic.tau = delayinyears * patientsperyear;

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
basic.PPatients = 100000  %total number of patients to be treated following adoption
basic.ICost = 1000     % fixed cost of investment

% computed data
basic.theta = exp( - log(1+annualdiscountrate) / patientsperyear ) %  per patient discount rate

%%%%%% Numerical computation parameters %%%%%%

basic.mumax =  15000;    % maximum posterior mean per patient; pick big enough to be above upper stopping boundary
basic.mumin = -40000;    % minimum posterior mean per patient; pick small enought to be below lower stopping boundary

%MinGridPerStdev = -0.004;    % if this is negative, code tries to use its negative as the delta-t
advanced.MinGridPerStdev = 100;    % minimum number of grid points per stdev of output, or negative to use dt=1

% For plotting
advanced.smoothed = true;   % true or false if curve smoothing is desired.
advanced.DOPLOT = false % put true or false depending if you do or do not want dynamic plots
advanced.NumGridsForContours = 100; % number of time values for generating contour plots (also num plots in animation to use during recursion over t)

advanced.bigfontsize=16;
advanced.smallfontsize=14;
advanced.fontname='Helvetica';

%titlestring='Total Hip Arthroplasty versus Resurfacing Arthroplasty';
advanced.titlestring='Test example';
advanced.filestring='Test';

end