function[ basic, advanced] = SetStents_sensit()
% SetSinus_SECTION5 gives an example of setting the input parameters for
% the delay sequential trials project. This set of parameters are used for
% the sinusitis application in Section 5 of the paper
%
% Designed for Delay Sequential Trials paper of Chick, Forster, Pertile
% (alpha order).
%
% (c) 2014, S Chick, P Pertile
% Created: 17 April 2014
% Last touched: 24 March 2015
% 
% The first two batches of parameters are most critical in terms of defining the
% problem structure of the clinical trial.
[basic, advanced] = DelayInputConstructor();
basic.PPatients = 100000*20;      
%basic.accrualperiod =  7 / 12 ;             
basic.delayinyears = 1 ;     
% studyyears = accrualperiod + delayinyears ; % MF 27/03 no longer reqd?
basic.numpairs = 529 ;          % number of pairwise allocations reported by study 
%postmean = 252;
%stder = (CIhi - postmean) / norminv( 1 - (1-Conflevel)/2); 
%basic.sigma = stder * sqrt(numpairs);
%basic.mu0=postmean;
basic.t0 = 20 ;                     % Assumption 
basic.sigma = 17538 ;               % Derived from Cohen et al. (2004). See also Table 1 of Pertile et al. (2010)
basic.annualdiscountrate = 0.01 ;         % 0.03 assumed 
basic.ICost = 0 ;                   % Assumption 
basic.TMax = 2000 ;              % Assumption 
basic.c = 200;                  % Assumption
basic.mu0 = 0;
basic.z  = 0 ; % critical value for classical statistical inference test (e.g. 1.96)
basic.mumax = 10 * 50 * 250000 / (basic.t0 + basic.TMax) ; % 500000 was the max s_Y in Martin-Paolo code for the same application
% basic.mumax = 4000 ; 
basic.mumin = -2 * basic.mumax ;
basic.patientsperyear = 907  ;  % maximum number of patients put into trial per year

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

advanced.MinGridPerStdev = 50;    % minimum number of grid points per stdev of output, or negative to use dt=1
% advanced.MinGridPerStdev = -0.01;    % if this is negative, code tries to use its negative as the delta-t
% advanced.MinGridPerStdev = 200; 

% For plotting
advanced.smoothed = true;   % true or false if curve smoothing is desired.
advanced.DOPLOT = false ; % put true or false depending if you do or do not want dynamic plots
advanced.NumGridsForContours = 100; % number of time values for generating contour plots (also num plots in animation to use during recursion over t)
advanced.bigfontsize=16;
advanced.smallfontsize=14;
advanced.fontname='Helvetica';
advanced.smallfontsize=14;

%advanced.titlestring='Stent Example from Paolo';
advanced.titlestring='Drug Eluting Stents versus Bare Metal Stents';
advanced.filestring='stent';
advanced.dirstring='stent';

end

%PARAMETER CAHNGES TO PRODUCE FIGURES STARTING FROM N.3 (replicating each figure requires to go back to the baseline value of the parameter, i.e. baseline value, after each experiment)

% Figure 3

% basic.c = 1;

% Figure 4

%basic.c = 500;

% Figure 5

%annualdiscountrate = 0;

% Figure 6

%annualdiscountrate = 0.03;

% Figure 7

%basic.t0 = 5;

% Figure 8

%basic.t0 = 80;

% Figure 9

%basic.online = true;






