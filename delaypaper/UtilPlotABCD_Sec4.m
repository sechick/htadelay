function tmp = UtilPlotABCD( basic, advanced, mat, ishoriz, ist0, lohivec )
% UtilPlotABCD: Plots the points A, B, C, D on the current figure, assuming
% those points are defined. Those points are computed by Stage I
% calculations (thus, the Stage I and Stage II calculations are presumed to
% have been run in order to compute mat from basic and advanced); the
% points A, B, C, D are indices into the vector mat.Threshpoint.
% The parameter ishoriz is set to be true if the means are on the x-axis,
% and false if means are on the y-axis. ishoriz defaults to true.
% lohivec is a length-2 vector with the low and y values to display for the
% lines at points A,B,C,D. it defaults to [0, tau].
% The parameter ist0 shifts A,B,C,D to t=t0 if true (default is t=0)
%
% This routine is intended to be called when the active figure has the
% number of samples as a function of the prior mean. The parameter ishoriz
% should be true (the default) if the means are on the x-axis and the
% number of samples is on the y-axis. If ishoriz is false, then it is
% presumed that the number of samples is on the x-axis and the values of
% the mean is on the y-axis.
%
% Project with Chick, Forster, Pertile (alpha order). 

% Original version: 19th December 2014.
% Updated: 11 Jan 2015 (SC)
% Updated: 24 Feb 2015 (MF): added functionality for ist0, plot
%   enhancements
% Updated: 11 Jan 2015 (SC): fixed default valus for ishoriz and ist0

    tmp = 1;
    % check parameter values
    if nargin < 4
        ishoriz = true;
    end
    if nargin < 5
        ist0 = false;
    end
    if nargin < 6
        lohivec = [0, basic.tau+1];
    end
    
    % MF 25/03: the following commands drop the relevant letters A, B, C, D in
    % cases where stage 1 does not exist. Code permits stage 1 not to exist
    % for AC but to exist for BD, to exist for BD but not for AC, not to
    % exist for either AC or BD, and to exist for both AC and BD. 
    if ( mat.Threshpoint( 1 ) == mat.Threshpoint( 3 ) + 1 && mat.Threshpoint( 2 ) == mat.Threshpoint( 4 ) - 1 ) 
        cellvec = { 'A', 'B', ' ', ' ' } ; 
    elseif ( mat.Threshpoint( 1 ) == mat.Threshpoint( 3 ) + 1 && mat.Threshpoint( 2 ) ~= mat.Threshpoint( 4 ) - 1 )
        cellvec = { ' ', 'B', ' ', 'D' };
    elseif ( mat.Threshpoint( 1 ) ~= mat.Threshpoint( 3 ) + 1 && mat.Threshpoint( 2 ) == mat.Threshpoint( 4 ) - 1 )
        cellvec = { 'A',' ', 'C', ' D' } ;
    else
        cellvec = { 'A', 'B', 'C', 'D' } ; 
    end

    for i=1:4 % for points A, B, C, D
        if mat.Threshpoint(i) > 0   % if there is a valid point (the points A, B, C, D might not exist: if not then threshpoint(i) will be negative.
            threshers = mat.Threshpoint(i);
            if i == 1 % 'shift in' points A and B on graph to account for the rounding to 0 in UtilInterpSVec() when values of mu between grid points have one of those grid points having a 0 optiml sample size for stage I sampling
                threshers = threshers -1;  % If UtilInterpSVec() had nonzero interpolation there, then these tweaks should be 0 instead of -1 and +1 respectively
            elseif i == 2
                threshers = threshers + 1;
            end

            if ishoriz   % if means are on x-axis, plot with means on x-axis
                if i > 2
                    vertalign = 'bottom';
                else
                    vertalign = 'bottom'; %'baseline';
                end
                if (i == 2) || (i == 3)
                    horizalign = 'right';
                else
                    horizalign = 'left';
                end
                scatter(mat.muvec(threshers) , lohivec(1), '+k');
                text(mat.muvec(threshers-1), lohivec(1) , cellvec{i},...
                    'HorizontalAlignment', horizalign, 'VerticalAlignment', vertalign,...
                    'FontSize',advanced.smallfontsize,'FontName',advanced.fontname);
                line( [ mat.muvec(threshers), mat.muvec(threshers)  ], [ lohivec(1), lohivec(2) ], ...
                    'Linestyle', ':' , 'Color', 'k'  ) ; % 
            else         % if means are on y-axis, plot means on y-axis...
                horizalign = 'left';
                if (i == 2) || (i == 3)
                    vertalign = 'top';
                else
                    vertalign = 'bottom';
                end
                scatter(lohivec(1) + ist0 * basic.t0, mat.muvec(threshers), '+k');
                text(lohivec(1) + ist0 * basic.t0+15, mat.muvec(threshers), cellvec{i},...
                    'HorizontalAlignment', horizalign, 'VerticalAlignment', vertalign,...
                    'FontSize',10,'FontName',advanced.fontname); % Fig 2-3
                   % 'FontSize',advanced.smallfontsize,'FontName',advanced.fontname);
                line( [ lohivec(1) + ist0 * basic.t0, lohivec(2) + ist0 * basic.t0 ], [ mat.muvec(threshers), mat.muvec(threshers)  ], ...
                    'Linestyle', ':', 'Color', 'k' ) ; %
            end
        end
    end
end