function [ figout, basicvec, advancedvec, legendvec, matvec ] = TestDelayIterate(fignum,basic,advanced,basicflag,fieldname,fieldvec,subtitle,filemodifier)
%TestDelayIterate Generates a bunch of graphs based on changing certain
%values of a single parameter. It takes the setup in basic and advanced. if
%basicflag is true, then it creates a vector of (basic,advanced)
%structures, changes the fieldname 'field' across them accoring to values
%found in 'fieldvec' (using basic structure if basicflag is true, and
%advanced structure if not). 'subtitle' is used to help with expanding the
%title field in the plot when displayed, 'filemodifier' is used to
%differentiate the name of the file if it is to be saved.
%
% Designed for Delay Sequential Trials paper of Chick, Forster, Pertile
% (alpha order).
%
% (c) 2014, S Chick
% Created: 17 April 2014
% Last touched: 17 April 2014 

[basicvec, advancedvec, legendvec, veclen] = UtilExperimentVectorCreate( basic, advanced, basicflag, fieldname, fieldvec);

% set up some file names in case files to differentiate output files if some are to be saved
for i=1:veclen
    advancedvec(i).filestring = [advancedvec(i).filestring filemodifier num2str(i)];
end

% run the analysis for each of the structures
for i=1:veclen
    if advancedvec(i).verbose 
        sprintf('running analysis %d: %s',i,advancedvec(i).filestring)
    end
    [rval, mat] = DelayCurvesRecur(basicvec(i), advancedvec(i));        % Do stage II of DP
    [mat] = DelayStageOne(basicvec(i), advancedvec(i), mat );           % Do stage I of DP
    tmp = advancedvec(i).DOPLOT; 
    advancedvec(i).subtitle = sprintf('%s / %s', subtitle, legendvec{i});
    [ fignum, mat ] = DelaySimOverview( fignum, basicvec(i), advancedvec(i), mat );  % Run simulations (if replications requested)
    if i==1 matvec = mat; else matvec = [matvec mat]; end
end

% do the plots based on the structures
for i=1:veclen
    advancedvec(i).subtitle = sprintf('%s / %s', subtitle, legendvec{i});
    fignum = DelayDoContours(fignum, basicvec(i), advancedvec(i), matvec(i));        % Generate a bunch of contour plots
end
%[fignum] = UtilExperimentVectorPlot( fignum, basicvec, advancedvec, legendvec, matvec, subtitle,filemodifier);

figout = fignum;

end

