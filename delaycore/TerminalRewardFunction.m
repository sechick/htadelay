function [rval] = TerminalRewardFunction(valuevec,discountfactor,predvarsample,nochangeallowed)
% (c) 2014 Chick, Forster, Pertile
% This file is for input to matlab, and does calculations
% to support the 'selecting the best system' paper in the
% bayesian environment.
% Edited 2014 to include improvement in computation value from Mills ratio
%
% Code is 'as is' and no guarantees for correctness.
%
% INPUTS:
% valuevec : vector of means of posterior distribution
% discountfactor: discountfactor to apply to reward after the patients are
% seen
% predvarsample: a scalar with the predictive variance of the posterior mean
% nochangeallowed: optional, SHOULD DEFAULT TO FALSE
%
% OUTPUTS:
% rval : vector of outputs
%
%
if nargin < 4
    nochangeallowed = false;
end
if nargin < 3
  'TerminalRewardFunction: not enough arguments'
end

if nochangeallowed % if no change is allowed, then it is equivalent to saying that there is no variance in posterior decision, or pred variance is 0
    predvarsample = 0;
end
if predvarsample > 0
    zvec = - valuevec / sqrt(predvarsample);
    rval = discountfactor^(1-nochangeallowed) * sqrt(predvarsample) * PsiNorm(zvec);
elseif predvarsample == 0       % this line indicates that if no change is allowed, the decision can be implemented without waiting for the delay of tau patients
    rval = (discountfactor^(1-nochangeallowed))*max(valuevec,0);
else
    rval=0*valuevec;
    'TerminalRewardFunction: called with negative predvarsample'    
end
