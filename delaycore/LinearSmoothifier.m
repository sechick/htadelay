function [ tout, valout ] = LinearSmoothifier( tvec, valvec, tighten )
%LinearSmoothifier: smoothens the stopping boundary from a grid based
%setting to a linearly smoothed, by taking the convex hull below the points
% INPUTS:
%   tvec: times at which the vector should be smoothed
%   valvec: vector containing the discretized stopping values
% OUTPUTS:
%   tout, valout: the smoothed time values and observation values
%
% This handles the case of jagged boundary and tries to 'smoothen' it (in
% case of numerical stabilties which sometimes arise with the upper
% boundary for the case of unknown variance. Essentially does a weighted
% average of nearby points.

% 1) Error checking and setup
tlen = length(tvec);
vlen = length(valvec);
if tlen ~= vlen
    'LinearSmootherLower: time and value vector not of same dimension'
    tout = tvec;
    valout = valvec;
    return;
end

% 2) Find the vector of 'jumps', and the times they jump.
tlen = length(tvec);
vlen = length(valvec);

tout = tvec(tlen);
valout = valvec(tlen);
for i=(tlen-1):-1:1
    if valvec(i) < valout(1)
        tout = [tvec(i), tout];
        valout = [valvec(i), valout];
    end
end
if tout(1) ~= tvec(1)
    % if the first time step is not at a jump, then need to do something a
    % bit special: here, we extrapolate linearly from preceding steps
    i=1;
    if length(tout) > 2
        ratio = (valout(2) - valout(1)) / (tout(2) - tout(1));
        newval = valout(2) - ratio * (tout(2) - tvec(1));
        tout = [tvec(i), tout];
        valout = [newval, valout];
    else
        'LinearSmootherLower: potentially odd smoothing at t0'
        tout = [tvec(i), tout];
        valout = [valvec(i), valout];
    end
end
    
% SECOND: Go back and find the linear interpolations if smoothening is
% needed
% cleanup: doublecheck if last entries are the same
toutlen = length(tout);
if length(tout) > 2
    if (tout(toutlen-1)==tout(toutlen))
        tout = tout(1:toutlen-1);
        valout(toutlen-1) = min(valout(toutlen-1:toutlen));
        valout = valout(1:toutlen-1);
    end
end
    
if ~tighten
    valout = interp1(tout,valout,tvec,'linear');  %
    tout = tvec;
end

end
