function [rval] = PsiNormUV(zval,dof,optionaltcdf)
% (c) 2004 Stephen E. Chick, all rights reserved
% This file is for input to matlab, and does calculations
% to support the 'selecting the best system' paper in the
% bayesian environment.
%
% Assumes unknown sampling variance in loss function.
%
% Code is 'as is' and no guarantees for correctness.
%
% INPUTS:
% zval : test statistic values (vector)
% dof : degrees of freedom.
% optionalzcdf : optional argument, if passed, it is assumed to be a
% precomputed version of tcdf(zval). Typically it will not be passed.

% OUTPUTS:
% rval : vector of outputs
%\frac{\nu + s^2}{\nu-1} \phi_\nu(s) - s (1- \Phi_\nu(s) )
%
% For BIG values of zval, one might try to use formula 
% based on Mills Ratio and Soms (1976) as cited in Soms (Jasa, 1980, vol
% 75, number 370, page 438, equation 2.1), or as in the 2001 Machine
% Learning article about student distribution mill's ratio approximations.
%
% NOTE: THIS APPROX IS TO BE VALIDATED.
%
    if nargin < 4
        needcdf = 1;
    else
        needcdf = 0;
    end

%    ZVALLIM = 100;
%    NOMILLSRATIO = true;   % for now, the Mills ratio stuff for student distribution is not debugged, so don't use it.

    tpdfvals = tpdf(zval, dof);
%    zbigpos=zval>ZVALLIM;
%    zbigneg=zval<-ZVALLIM;
%    zmed=abs(zval)<=ZVALLIM;
    
%    if NOMILLSRATIO
%        zbigpos = 0*zbigpos;
%        zbigneg = zbigpos;
%        zmed=ones(size(zmed));
%    end

    %rbig = zbig .* tpdfvals .* (1./zval - dof ./ ((dof+2) * zval.^3)) ./ zval; % uses Mills ratio for improving stability when z is big
    %rbig = zbig .* tpdfvals .* (1./zval.^2 - dof ./ ((dof+2) * zval.^4)) ; % uses Mills ratio for improving stability when z is big
 %   rbigpos = zbigpos .* tpdfvals .* ((dof + zval.^2) ./ (dof - 1) - zval .* sqrt(1+ zval.^2/dof) * (1/2 + 1/sqrt(dof))) ; % uses Mills ratio for improving stability when z is big
 %   rbigpos(isnan(rbigpos))=0;    % handle case of divide by 0 in preceding line.
 %   rbigneg = zbigneg .* (-zval + tpdfvals .* ((dof + zval.^2) ./ (dof - 1) - zval .* sqrt(1+ zval.^2/dof) * (1/2 + 1/sqrt(dof)))) ; % uses Mills ratio for improving stability when z is big
 %   rbigneg(isnan(rbigneg))=0;    % handle case of divide by 0 in preceding line.
    if needcdf
%        rmed = zmed .* ((dof + zval.^2) .* tpdfvals ./ (dof - 1) - zval .* tcdf( - zval, dof));
        rval = (dof + zval.^2) .* tpdfvals ./ (dof - 1) - zval .* tcdf( - zval, dof);
    else
%        rmed = zmed .* (tpdfvals - zval .* (1-optionaltcdf));
        rval = tpdfvals - zval .* (1-optionaltcdf);
    end
%    rval = rbigpos + zbigneg + rmed;

end