function [basic, advanced] = SetHipArthro() 
% SETHIPARTHRO gives and example of setting the input parameters for
% the delay sequential trials project. 
%
%%%% General parameters to set. For delay clinical trial paper.
% EXAMPLE 1: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1a.	Cost-effectiveness of total hip arthroplasty versus resurfacing arthroplasty.	Edlin et al.	2012	BMJ open	2		CEA	Incremental C, QALY, INMB
% 1b.	Total hip arthrplasty versus resurfacing. 	Costa et al.	2012	BMJ	344		(Clinical trial results for 1a)	
% UNIT: POUNDS UK
%
% Designed for Delay Sequential Trials paper of Chick, Forster, Pertile
% (alpha order).
%
% Created: 14 April 2014
% Last touched: 14 April 2014
% 

%%%%%% Problem dependent parameters %%%%%%
% 95% CI (-1795, 1960)
% mean: 82.46
% numpatients: 60 in one arm 66 in the other: so, lets use 60.
% delayinyears = 1;
% length of study : 3 years
% pick TMax bigger.
% discount: was not discounted, 
[basic, advanced] = DelayInputConstructor();

% For info: ICER rate for study was 20 000 UK pounds per QALY.
% QALYs were assumed normally distributed, conditional on trial arm and
% EQ-5D-3L value. icrement measured at 12 months.
% total costs over 12 months assumed to be lognormal.
% 

%CIlow=-1795; 
CIhi = 1960;
Conflevel = 0.95;
numpairs = 60;
postmean = 82.46;
stder = (CIhi - postmean) / norminv( 1 - (1-Conflevel)/2);
studyyears = 3;
delayinyears = 1;     % one month is 1/12.

basic.numpairs = numpairs;
basic.sigma = stder * sqrt(numpairs);
basic.mu0=postmean;
basic.t0 = 2;
%sigma0 = basic.sigma/sqrt(basic.t0);   % prior: std dev for unknown mean

% Patient delay and pool
basic.TMax = 300;
annualdiscountrate = .03;  % real valued number, with value like 0.05 for 5% discount rate, paper used 1.9% for deflation, otherwise no discounting
patientsperyear = basic.numpairs / (studyyears - delayinyears);  % maximum number of patients put into trial per year
basic.c = 200;  % variable cost per patient of doing trial
basic.tau = delayinyears * patientsperyear;
advanced.z  = 0; % critical value for classical statistical inference test (e.g. 1.96)


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
basic.PPatients = 70000 * 0.06 * 10  % 70K is number of hip arthroplasty ops per year in England and Wales, of which 6% are hip resurfacings, and we assume a 10 year horizon for contract
basic.ICost = 0     % fixed cost of investment: this is ad hoc choice

% computed data
basic.theta = exp( - log(1+annualdiscountrate) / patientsperyear ) ; %  per patient discount rate

%%%%%% Numerical computation parameters %%%%%%

basic.mumax =  3e4;    % maximum posterior mean per patient; pick big enough to be above upper stopping boundary
basic.mumin = -3e4;    % minimum posterior mean per patient; pick small enought to be below lower stopping boundary

%MinGridPerStdev = -0.004;    % if this is negative, code tries to use its negative as the delta-t
advanced.MinGridPerStdev = 120;    % minimum number of grid points per stdev of output, or negative to use dt=1

% For plotting
advanced.smoothed = true;   % true or false if curve smoothing is desired.
advanced.DOPLOT = false; % put true or false depending if you do or do not want dynamic plots
advanced.NumGridsForContours = 100; % number of time values for generating contour plots (also num plots in animation to use during recursion over t)

advanced.bigfontsize=16;
advanced.smallfontsize=14;
advanced.fontname='Helvetica';

%titlestring='Total Hip arthroplasty versus Resurfacing arthroplasty';
advanced.titlestring='Total Hip Arthroplasty';
advanced.filestring='Hip';
advanced.dirstring = 'hip';    % root/base of name of directory for output files (figures, data) related to this analysis
    
