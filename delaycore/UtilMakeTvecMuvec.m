function [tvec, muvec] =  UtilMakeTvecMuvec( basic, advanced )
%UtilMakeTvecMuvec returns the time grid and mu grid, given the problem structures defined in the parameters basic, advanced.
%
%
%

% For use in Delay sequential sampling project.
tlen = 1+max(1,ceil( (basic.TMax - basic.tau) / advanced.dt ));
tvec= (basic.t0 + basic.tau) + advanced.dt * (0:(tlen-1));
muvec = advanced.mushift+advanced.dmu *(round(basic.mumin/advanced.dmu):round(basic.mumax/advanced.dmu));        % handle integer number of grid points, cover range of mu needed

end
