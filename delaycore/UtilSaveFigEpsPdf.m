function [ output_args ] = UtilSaveFigEpsPdf( h, dirname, fname, varargin )
%UtilSaveFigEpsPdf: saves content of figure number 'h' as a fig file, an eps
%   file, and pdf file, using the base name 'fname', in directory 'dirname', 
%   and with print flags 'varargin'. The argument 'varargin' is optional, 
%   and may actually have a variable number of arguments.

    if ~isempty(dirname)      
        tmp = dir(dirname);
        if ~length(tmp) %~isdir(dirname)
            rval = mkdir(dirname);
        end
        fullpath = sprintf('%s/',dirname);
    else
        fullpath = '';
    end
    fullpath = [fullpath fname];

    nVarargs = length(varargin);

    % Save EPS in regular format: have not yet fully figured out how to get
    % the bounding box 'tight' for EPS format
    % This saves EPS in middle of full size page ... lots of space around
    % it
   
    
    if nVarargs < 1
        print(fullpath,'-depsc');
    else
        print(fullpath,'-depsc',varargin{1:nVarargs});
    end

    % the following seems to work for getting the PDF and FIG to be tighter
    % for the bounding box
% get the current axes
ax = get(h, 'CurrentAxes');

% make it tight
ti = get(ax,'TightInset');
set(ax,'Position',[ti(1) ti(2) 1-ti(3)-ti(1) 1-ti(4)-ti(2)]);

% adjust the papersize
set(ax,'units','centimeters');
pos = get(ax,'Position');
ti = get(ax,'TightInset');
set(h, 'PaperUnits','centimeters');
set(h, 'PaperSize', [pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
set(h, 'PaperPositionMode', 'manual');
set(h, 'PaperPosition',[0 0  pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
    
	saveas(h,fullpath) ; 
% SEC:
%    saveas(h,fullpath,'eps2'); % this saves EPS, in lower left corner, but box is too big
% PERHAPS use GS tools to 'clip' the EPS file automatically?
    if nVarargs < 1
        print(fullpath,'-dpdf');
    else
        print(fullpath,'-dpdf',varargin{1:nVarargs});
    end
    
end

