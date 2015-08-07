function [sval] = UtilInterpSVec(muvec,mu0,bestsvec,tau,TMax) %P have added Tmax
% For a vector muvec of mu values and a vector of samples to take in one stage
% sampling (bestsvec), interpolate to find the optimal number of samples to
% take in stage 1 if the mean is mu0. Assumes delay of tau in selecting
% samples. If tau is set to something negative, then the effect of tau is
% not used in interpolating (meaning, set tau to be -1 for example if
% interpolating to obtain the optimal one-stage sample size).
%
% Returns sval as an integer (rounded up if there was a real valued sample
% size)
%
% NB: Assumes muvec runs up from low values to high values.
%
% if mu0 is at or below lowest value in muvec, or above highest value, return 0 samples as optiomal

if tau < 0
    IsOneStageExperiment = true;
else
    IsOneStageExperiment = false;
end

%    stage2true = (bestsvec == (tau+1)); % get booleans to see if one progresses to stage 2 or not
%    nonzerotrue = (bestsvec > 0);       % get booleans to see if one should take at least to stage 1
%    stage1not2 = nonzerotrue & ~stage2true;    % booleans if one should do stage 1 but not stage 2

    if mu0 < muvec(1)
        sval = bestsvec(1);   
    elseif mu0 > max(muvec)
        sval = bestsvec(length(bestsvec));
    else % ok have a valid mu0 in the resonable range.
        % first, check if mu0 is exactly one of the entries in muvec. if so,
        % can return without interpolation
        [mu0inmuvec, binx] = max(muvec == mu0);
        if mu0inmuvec  % found exact match, no interpolation needed
            sval = bestsvec(binx);
        elseif IsOneStageExperiment % if we should ignore the effect of tau, just go ahead an interpolate
            %% PAOLO CHANGES START
            [~, binx] = max(muvec > mu0);    % find index of first element greater than mu0
            lowx = binx-1;          % and index of first element less than mu0
            if bestsvec(lowx) == bestsvec(binx) % if the optimal svec values are equal above and below
                sval = bestsvec(lowx);          % then can set sval equal to that value
            elseif min(bestsvec(lowx),bestsvec(binx)) == 0  % if either at 0 then set to 0
                sval = 0;       % NB: Could set to 'max' here to insure more sampling for stage I... to be tested.
            elseif (max(bestsvec(lowx),bestsvec(binx)) <= TMax) % if both between 0 and TMax, interpolate
                sval = interp1(muvec,bestsvec,mu0,'linear');
            elseif min(bestsvec(lowx),bestsvec(binx)) > TMax % if both greater than TMax, interpolate
                sval = interp1(muvec,bestsvec,mu0,'linear');
            else   % if one less than TMAX and the other greater than TMax, pick the MINIMUM (%P). this prevents interpolating when one is less than TMax and the other is greater than TMax
                sval = min(bestsvec(lowx),bestsvec(binx));
            end
            %% PAOLO CHANGES END
           % sval = interp1(muvec,bestsvec,mu0,'linear'); %P this was the
           % original line
        else    % did not find an exact match: mu0 between grid points. check if interpolation is needed, accounting for effect of 'tau' gap
            [~, binx] = max(muvec > mu0);    % find index of first element greater than mu0
            lowx = binx-1;          % and index of first element less than mu0
            if bestsvec(lowx) == bestsvec(binx) % if the optimal svec values are equal above and below
                sval = bestsvec(lowx);          % then can set sval equal to that value
            elseif min(bestsvec(lowx),bestsvec(binx)) == 0  % if either at 0 then set to 0
                sval = 0;       % NB: Could set to 'max' here to insure more sampling for stage I... to be tested.
            elseif (max(bestsvec(lowx),bestsvec(binx)) <= tau) % if both between 0 and tau, interpolate
                sval = interp1(muvec,bestsvec,mu0,'linear');
            elseif min(bestsvec(lowx),bestsvec(binx)) > tau % if both greater than tau, interpolate
                sval = interp1(muvec,bestsvec,mu0,'linear');
            else   % if one less than tau and the other greater than tau, pick the MINIMUM (%P). this prevents interpolating when one is less than tau and the other is greater than tau
                sval = min(bestsvec(lowx),bestsvec(binx));
            end
        end
    end

    sval = ceil(sval);  % round up the number of samples

end