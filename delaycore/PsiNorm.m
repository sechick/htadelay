function [rval] = PsiNorm(zval,optionalzcdf)
% (c) 2004 Stephen E. Chick, all rights reserved
% This file is for input to matlab, and does calculations
% to support the 'selecting the best system' paper in the
% bayesian environment.
% Edited 2014 to include improvement in computation value from Mills ratio
%
% Code is 'as is' and no guarantees for correctness.
%
% INPUTS:
% zval : test statistic values (vector)
% optionalzcdf : optional argument, if passed, it is assumed to be a
% precomputed version of normcdf(zval). Typically it will not be passed.
%
% OUTPUTS:
% rval : vector of outputs
%
if nargin < 2
    needcdf = 1;
else
    needcdf = 0;
end

ZVALLIM = 10;
zpdf = normpdf(zval);

zbig=zval>ZVALLIM;
zsmall=zval<=ZVALLIM;
rbig = zbig .* (zpdf ./ (zval.*zval + 1)); % uses Mills ratio for improving stability when z is big
if needcdf
    rsmall = zsmall .* (zpdf - zval .* normcdf( - zval));
else
    rsmall = zsmall .* (zpdf - zval .* (1-optionalzcdf));
end
rval = rbig + rsmall;

% following naive version is mathematically correct but computationally unstable for large |zval|
%rval = normpdf(zval) - zval .* normcdf( - zval);
%rval = max(rval,0);
