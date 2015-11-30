function [rval] = TerminalRewardFunctionUnk(valuevec,discountfactor,predvarsample,nochangeallowed, dof)
% (c) 2014 Chick, Forster, Pertile
% This file is for input to matlab, and does calculations
% to support the 'selecting the best system' paper in the
% bayesian environment.
% Assumes sampling variance is unknown, so t-distribution loss function is used.
%
% Code is 'as is' and no guarantees for correctness.
%
% INPUTS:
% valuevec : vector of means of posterior distribution
% discountfactor: discountfactor to apply to reward after the patients are
% seen
% predvarsample: a scalar with the predictive variance of the posterior mean
% nochangeallowed: optional, SHOULD DEFAULT TO FALSE
% dof: degrees of freedom
% OUTPUTS:
% rval : vector of outputs
%
%
if nargin < 5
  warning('TerminalRewardFunctionUnk: not enough arguments');
end

if nochangeallowed % if no change is allowed, then it is equivalent to saying that there is no variance in posterior decision, or pred variance is 0
    predvarsample = 0;
end
if predvarsample > 0
    zvec = abs(valuevec) / sqrt(predvarsample);
    rval = sqrt(predvarsample) * PsiNormUV(zvec,dof) + valuevec .* (valuevec>0);
    rval = discountfactor^(1-nochangeallowed) * rval;
elseif predvarsample == 0       % this line indicates that if no change is allowed, the decision can be implemented without waiting for the delay of tau patients
    rval = (discountfactor^(1-nochangeallowed))*max(valuevec,0);
else
    rval=0*valuevec;
    warning('TerminalRewardFunction: called with negative predvarsample');  
end
