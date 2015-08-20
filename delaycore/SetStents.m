function [basic, advanced] = SetStents()
% SETSTENTS gives and example of setting the input parameters for
% the delay sequential trials project. Makes parameters so that in
% principle they match the problem structure of the example used in the
% code developed by Paolo.
%
% Designed for Delay Sequential Trials paper of Chick, Forster, Pertile
% (alpha order).
%
% (c) 2014, S Chick, P Pertile
% Created: 17 April 2014
% Last touched: 17 April 2014
% 
% The first two batches of parameters are most critical in terms of defining the
% problem structure of the clinical trial.
[basic, advanced] = DelayInputConstructor();

basic.PPatients = 2000000;
studyyears = 7/12 + 1; % recruitment lasted 7 months, follow-up is one year
delayinyears = 1;     % one month is 1/12.
%CIhi = ;
%Conflevel = ;
basic.numpairs = 529;
postmean = 252;
%stder = (CIhi - postmean) / norminv( 1 - (1-Conflevel)/2); 
%basic.sigma = stder * sqrt(numpairs);
%basic.mu0=postmean;
basic.t0 = 20;
basic.sigma = 17538; % I'LL TRY TO SEE WHETHER THIS IS CONSISTENT WITH WHAT ONE GETS BY USING THE SAME APPROACH AS IN THE OTHER FILE
annualdiscountrate = .01;  % real valued number, with value like 0.05 for 5% discount rate
basic.ICost = 0;
patientsperyear = basic.numpairs / (studyyears - delayinyears);  % maximum number of patients put into trial per year
basic.tau = round(delayinyears * patientsperyear);
basic.TMax = 2000;
basic.c = 200;
%sigma0 = basic.sigma/sqrt(basic.t0); % no need to assign, will be computed
%in routines from sigma and t0 as needed.
basic.mumax = 50 * 500000 / (basic.t0 + basic.TMax); % 500000 was the max s_Y in Martin-Paolo code for the same application
basic.mumin = - 2 * basic.mumax; % -2
basic.mu0=0;
advanced.z  = 0; % critical value for classical statistical inference test (e.g. 1.96)
basic.theta = exp( - log(1+annualdiscountrate) / patientsperyear ); %  per patient discount rate
%%%%%%%%%
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

%%%%%% Numerical computation parameters %%%%%%

advanced.MinGridPerStdev = 30;    % minimum number of grid points per stdev of output, or negative to use dt=1
% advanced.MinGridPerStdev = 200; 

% For plotting
advanced.smoothed = true;   % true or false if curve smoothing is desired.
advanced.DOPLOT = false; % put true or false depending if you do or do not want dynamic plots
advanced.NumGridsForContours = 80; % number of time values for generating contour plots (also num plots in animation to use during recursion over t)

advanced.bigfontsize=16;
advanced.smallfontsize=14;
advanced.fontname='Helvetica';

advanced.titlestring='Stent Example from Paolo';
advanced.filestring='stent';
advanced.dirstring='stent';

end