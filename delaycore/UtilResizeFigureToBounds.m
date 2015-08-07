function tmp = UtilResizeFigureToBounds( mat )
% UtilResizeFigureToBounds: resets the x-axis and y-axis bounds, assuming
% that time values are on the x-axis, and posterior mean values are on the
% y-axis, so that the x-axis displays values from 0 to the maximum time
% value in the structure 'mat' from the Delay differential equation
% routines, and the y-axis contains a range that is just a bit wider than
% the range from the minimum of the lower bound to the maximum of the
% upper bound.
%
% If mat is a vector of such structures, as in a comparison of several
% analysis with different input parameters, then the bounds are big enough
% to account for the largest such range.

len = length(mat);

GAPPARAM = 16;  % parameter to control how much above the top boundary and how much below the bottom boundary
% one should keep in plot: bigger GAPPARAM means smaller amount on top and bottom. 

tmp=axis(); 
i=1;
tmp(1) = min(0, min(mat(i).tvec));
tmp(2) = max(0, max(mat(i).tvec));
tst = max(mat(i).bndupper) - min(mat(i).bndlower);
tmp(3)=min(mat(i).bndlower)-tst/GAPPARAM; 
tmp(4)=max(mat(i).bndupper)+tst/GAPPARAM; 
for i=2:len
    tmp(1) = min(tmp(1), min(mat(i).tvec));
    tmp(2) = max(tmp(2), max(mat(i).tvec));
    tst = max(mat(i).bndupper) - min(mat(i).bndlower);
    tmp(3)=min(tmp(3), min(mat(i).bndlower)-tst/GAPPARAM); 
    tmp(4)=max(tmp(4), max(mat(i).bndupper)+tst/GAPPARAM); 
end
    
axis(tmp);

end

