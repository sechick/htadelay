function [rval] = TerminalRewardFunctionVecUnk(valuevec,discountvector,predvarvec,nochangeallowed, dof)
% (c) 2014 Chick, Forster, Pertile
% This file is for input to matlab, and does calculations
% to support the 'selecting the best system' paper in the
% bayesian environment.
% Assumes sampling variance is unknown, so t-distribution loss function is used.
%
% does not yet include mills ratio improvement.
%
% Code is 'as is' and no guarantees for correctness.
%
% INPUTS:
% valuevec : ROW VECTOR of means of posterior distribution
% discountfactor: VECTOR of discount factors to apply to reward after the patients are
% seen (should be same length 
% predvarvec: a VECTOR with the predictive variance of the posterior mean
% nochangeallowed: optional, SHOULD DEFAULT TO FALSE
% dof: degrees of freedom
%
% OUTPUTS:
% rval : MATRIX of rewards, with number of rows equal to size of vector of
% means, and number of columns equal to the number of discount factors.
%
%
if nargin < 5
  'TerminalRewardFunctionVecUnk: not enough arguments'
end

if nochangeallowed % if no change is allowed, then it is equivalent to saying that there is no variance in posterior decision, or pred variance is 0
    rval = max(valuevec' * discountvector,0);
else
    zvec = - valuevec' * (1./sqrt(predvarvec));
    rval = PsiNormUV(zvec,dof); % pick up valid values first
    rval(isnan(rval))=0;rval(isinf(rval))=0;
    rval = rval * diag(discountvector .* sqrt(predvarvec));
    rval = rval + valuevec' * (predvarvec==0);
end
