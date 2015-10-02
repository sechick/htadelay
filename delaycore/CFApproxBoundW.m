function [wval] = CFApproxBoundW(sval)
% (c) 2009 Stephen E. Chick, all rights reserved
% This file is for input to matlab, and does calculations
% to support the chick frazier paper
%
% Code is 'as is' and no guarantees for correctness.
%
% Updated 2015 Aug 14, to bring in Chernoff Sequential test III result for
% boundary for very large sval.
%
% INPUTS:
% sval : s value for standardized problem (can be a vector)
%
% OUTPUTS:
% wval : approximation to boundary

    wval = sval;
    wval(sval <= 1) = 0.233*sval(sval<=1).^2;
    wval(sval>1 & sval<=3) = 0.00537*sval(sval>1 & sval<=3).^4 - 0.06906*sval(sval>1 & sval<=3).^3 + ...
        0.3167*sval(sval>1 & sval<=3).^2 - 0.02326*sval(sval>1 & sval<=3);
    wval(sval>3 & sval<=40) = 0.705*sqrt(sval(sval>3 & sval<=40)).*log(sval(sval>3 & sval<=40));
    wval(sval > 40) = sqrt( sval(sval > 40) .* (3*log(sval(sval > 40))) - log(8*3.14159) - 2./log(sval(sval > 40)) - 170);

end     