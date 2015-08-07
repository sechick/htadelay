function [wval] = CGApproxBoundW(sval)
% (c) 2009 Stephen E. Chick, all rights reserved
% This file is for input to matlab, and does calculations
% to support the chick gans paper
%
% Code is 'as is' and no guarantees for correctness.
%
% INPUTS:
% sval : s value for standardized problem (can be a vector)
%
% OUTPUTS:
% wval : approximation to boundary

    wval = sval;
    wval(sval <= (1/7)) = sval(sval<= (1/7)) / sqrt(2);
% replacing formula for middle patch with a slightly different version
% which matches the endpoints better (but might call for slightly earlier
% stoppping). If the formula is changed further, please also make the
% change in CGApproxBoundWinv.m
%    wval(sval>(1/7) & sval<=100) = exp( -0.02645*log( sval(sval>(1/7) & sval<=100) ).^2 + ...
%        0.89106* log( sval(sval>(1/7) & sval<=100) ) - 0.4873 );
    wval(sval>(1/7) & sval<=100) = exp( -0.0275*log( sval(sval>(1/7) & sval<=100) ).^2 + ...
        0.8797 * log( sval(sval>(1/7) & sval<=100) ) - 0.5024 );
    wval(sval > 100) = sqrt(sval(sval > 100)) .* sqrt( 2* log(sval(sval > 100)) - ...
        log(log(sval(sval > 100))) - log(16*pi));

end    