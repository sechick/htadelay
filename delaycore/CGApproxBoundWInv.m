function [sval] = CGApproxBoundWInv(mval)
% (c) 2009 Stephen E. Chick, all rights reserved
% This file is for input to matlab, and does calculations
% to support the chick gans paper
%
% Code is 'as is' and no guarantees for correctness.
%
% INPUTS:
% mval : set of scaled values to invert through boundary, assuming c=0, delta > 0
%
% OUTPUTS:
% sval : set of s values which correspond to inverting boundary: sval =
% b-wiggle-inverse(mval)
%

kink1 = 1/7; % these thresholds are picked as the 'kinks' in CGApproxBoundW function
kink2 = 100;
breakpnt1 = CGApproxBoundW(kink1);    
breakpnt2 = CGApproxBoundW(kink2);

sval = mval;       % allocate space
sval(mval <= 0) = Inf;     % nonpositive values have infinite inversse

sval(mval <= breakpnt1) = sqrt(2) * mval(mval <= breakpnt1); % first patch can do analytically

if sum((mval > breakpnt1) & (mval <= breakpnt2))   % if there are points in the middle patch...
    loval = 0.99*kink1; % pick values just outside of band for middle patch
    hival = 1.01*kink2; % just in case there is a bit of roundouff at the kinks
    numsteps = 600;
    myvec = exp( log(loval) + (log(hival)-log(loval))*(0:numsteps)/numsteps);  % space things out on log scale
    % in the following formula, an updated version of the boundary is given
    % (relatvie to the CG paper) which does a better job of matching the
    % boundaries at the 'kinks', but it does call for slightly earlier
    % stopping times in this middle band (relative to stopping times in the
    % CG paper).
    bndsvec = exp( -0.0275*log( myvec ).^2 + 0.8797* log( myvec ) - 0.5024 ); % needs to match formula in CFApproxBoundW
    sval((mval > breakpnt1) & (mval <= breakpnt2))=...
        interp1(bndsvec,myvec,mval((mval > breakpnt1) & (mval <= breakpnt2)));
end

if sum(mval > breakpnt2)   % if there are points in the third patch...
    % find an upper bound for the boundary so that we are sure to get a valid
    % interpolation: that is, so that kink2 and kink3 give the support of
    % values of s for all possible values of mval in the third patch....
    goldenratio = ( 1 + sqrt(5) ) / 2;
    kink3 = goldenratio*kink2;
    deltanumsteps = 200;
    numsteps = deltanumsteps;
    while max(mval) > CGApproxBoundW(kink3)
        kink3 = kink3 * goldenratio;
        numsteps = numsteps + deltanumsteps;
    end
    
    loval = 0.99*kink2; % pick values just below the third patch, and as big as the biggest value of the boundary needed for inversion
    hival = kink3; % just in case there is a bit of roundouff at the kinks
    myvec = exp( log(loval) + (log(hival) - log(loval))*(0:numsteps)/numsteps );
    myvec(length(myvec)) = kink3;
    bndsvec = sqrt(myvec) .* sqrt( 2* log(myvec) - log(log(myvec)) - log(16*pi));
    sval(mval > breakpnt2)=interp1(bndsvec,myvec,mval(mval > breakpnt2));
end

end     