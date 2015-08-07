function [rval] = DelayContinExpectation(first,last,pu,pd,psame,myvec)
% (c) 2014 Chick
%
% For Chick, Forster, Pertile paper on 'delay'.
%
% Edited 2014 to include improvement in computation value from Mills ratio
%
% Code is 'as is' and no guarantees for correctness.
%
% INPUTS:
% first,last: integers with first and last indices of elements of myvec
% which are considered to be in the continuation set
% pu, pd, psame: probability of going 'up', 'down' or 'stay same' in
% trinomial tree in backward recursion
% myvec : the vector upon which backward recursion is to be operated on
%
% OUTPUTS:
% rval : vector of outputs which is the expected reward mixing of the 'up'
% and 'down' parts
%
% NOTES: Assumes values with 'up' correspond to higher indices in the
% vector 'myvec'
%

    if (first <= 1) || (last >= length(myvec))
        warning('DelayContinExpectation: assumes first and last entries of stopped are true. error in calling routine?');
    end
    if first > last
        warning('DelayContinExpectation: assumes someplace is stopped');
    end
    
    try
        rval = pu*myvec((first+1):(last+1)) + psame*myvec(first:last) + pd * myvec((first-1):(last-1));
    catch
        [first, last, length(myvec)]
        rval = myvec;
    end

end