function [thresvec] = UtilGetABCDVector(bestsvec,tau)
% Try to compute points A, B, C, D, if they are defined, where A, B, C, D
% are all indices into the vector muvec, and refer to points where muvec(.)
% takes on a value with a particular property. A and B are defined if and
% only if there is value of mu in muvec for which one samples at least
% once. C and D are only defined if and only if there is a value of mu in
% muvec for which one continues to stage II sampling.
%   A: minimum value of mu in muvec which is above the zone where it is
%   optimal to sample at least once, and for which 0 samples are optimal
%   B: maximum value of mu in muvec which is below the zone where it is
%   optimal to sample at least once, and for which 0 samples are optimal
%   C: minimum value of mu in muvec which is above the zone where it is
%   optimal to sample at least tau times, and for which one samples more than
%   0 times (or, if this is the same point as A, then the mu just below A
%   in muvec.
%   D: similar rules as for C.
%
% If these points are not defined, then the value of the assigned index is
% -1.
%
% These points are stored in the vector mat.Threshpoint (length 4).

    thresvec = -1*ones(4,1); % by default, the points are not defined unless shown otherwise
    len = length(bestsvec);
    if sum(bestsvec>0) > 0 % A and B might also be defined
        if bestsvec(len) == 0   % need to check if the upper bound is actually in the grid of muvec: if not, should run with higher values of mu in muvec
            [testval, testindx] = max(fliplr(bestsvec) > 0);
            thresvec(1) = len - testindx + 2;
        end
        if bestsvec(1) == 0   % need to check if the lower bound is actually in the grid of muvec: if not, should run with lower values of mu in muvec
            [testval, testindx] = max(bestsvec > 0);
            thresvec(2) = testindx-1;
        end
        if sum(bestsvec>tau) > 0 % C and D might be defined
            if bestsvec(len) <= tau   % need to check if the upper bound is actually in the grid of muvec: if not, should run with higher values of mu in muvec
                [testval, testindx] = max(fliplr(bestsvec) > tau);
                thresvec(3) = len - testindx + 2;
                if thresvec(1)==thresvec(3)
                    thresvec(3) = thresvec(1) - 1;
                end
            end
            if bestsvec(1) <= tau   % need to check if the lower bound is actually in the grid of muvec: if not, should run with lower values of mu in muvec
                [testval, testindx] = max(bestsvec > tau);
                thresvec(4) = testindx-1;
                if thresvec(2)==thresvec(4)
                    thresvec(4) = thresvec(2) + 1;
                end
            end
        end
    end
end
