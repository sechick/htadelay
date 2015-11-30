function [rval] = TerminalProbPickNewUnk(valuevec,predvarsample,nochangeallowed,dof)
% (c) 2014 Chick, Forster, Pertile
% This file is for input to matlab, and does calculations
% to support the 'selecting the best system' paper in the
% bayesian environment.
% Call if variance is unknown
%
% Code is 'as is' and no guarantees for correctness.
%
% INPUTS:
% valuevec : vector of means of posterior distribution
% predvarsample: a scalar with the predictive variance of the posterior mean
% nochangeallowed: optional, SHOULD DEFAULT TO FALSE
% dof: degrees of freedom for relevant t distribution
%
% OUTPUTS:
% rval : vector of outputs
%
if nargin < 4
  rval=0;
  warning('TerminalProbPickNew: not enough arguments');
end

if nochangeallowed % if no change is allowed, then it is equivalent to saying that there is no variance in posterior decision, or pred variance is 0
    predvarsample = 0;
end
if predvarsample > 0
    zvec = valuevec / sqrt(predvarsample);
    rval = tcdf(zvec,dof);
elseif predvarsample == 0       % this line indicates that if no change is allowed, the decision can be implemented without waiting for the delay of tau patients
    rval = (valuevec >= 0);
else
    rval = 0;
    warning('TerminalProbPickNewUnk: called with negative predvarsample');
end

end
