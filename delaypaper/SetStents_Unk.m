function[ basic, advanced] = SetStents_Unk()
% SetStents_Unk gives an example of setting the input parameters for
% the delay sequential trials project. This set of parameters are used for
% the sinusitis application in Section 5 of the paper
%
% Designed for Delay Sequential Trials paper of Chick, Forster, Pertile
% (alpha order).
%
% (c) 2016, S Chick, 2 april 2016
%
% The first two batches of parameters are most critical in terms of defining the
% problem structure of the clinical trial.

% use same parameters as for the stent example with known variance...
[basic, advanced] = SetStents_SECTION5(); 
% except for resetting a few parameters ...
accrualperiod =  7 / 12 ;             
patientsperyear = basic.numpairs / accrualperiod  ;  % maximum number of patients put into trial per year

% do for undiscounted
annualdiscountrate = 0.01 ;         % Assumption
basic.theta = exp( - log(1+annualdiscountrate) / patientsperyear ); %  per patient discount rate
% and reset for unknown variance with a particular setting for shape
% parameter for unknown variance.
advanced.UnkVariance = true;
advanced.UnkVarianceShape = 2*basic.t0;   % if = basic.t0/2, then uninformative prioradvanced.UnkVarianceShape = 2*basic.t0;  %P: original:2*basic.t0 % if = basic.t0/2, then uninformative prior
%basic.c = 800 ;                  % Assumption
% basic.sigma0 = basic.sigma/sqrt(basic.t0) ; no need to assign, will be

% redo a few parameters for computation and for output
advanced.MinGridPerStdev = 35;    % minimum number of grid points per stdev of output, or negative to use dt=1
% advanced.MinGridPerStdev = 200;
advanced.NumGridsForContours = 60; % number of time values for generating contour plots (also num plots in animation to use during recursion over t)

advanced.titlestring='DES versus BMS with Unknown Variance';
advanced.filestring='stentunk';
advanced.dirstring='stentunk';

end
