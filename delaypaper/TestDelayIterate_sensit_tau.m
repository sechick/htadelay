function [ figout, basicvec, advancedvec, legendvec, matvec ] = TestDelayIterate_sensit_tau(fignum,basic,advanced,basicflag,fieldname,fieldvec,subtitle,filemodifier)
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

[basicvec, advancedvec, legendvec, veclen] = UtilExperimentVectorCreate_sensit_tau( basic, advanced, basicflag, fieldname, fieldvec);

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
    %%%%
   if i==1
        advancedvec(i).plot_upper = mat(i).muvec(mat(i).Threshpoint(1)) + min(mat(i).muvec(mat(i).Threshpoint(1))/4,- mat(i).muvec(mat(i).Threshpoint(2))/4) ; % set upper limit of plot region ( + /  %P
        advancedvec(i).plot_lower =  mat(i).muvec(mat(i).Threshpoint(2)) - min(-mat(i).muvec(mat(i).Threshpoint(2))/4,mat(i).muvec(mat(i).Threshpoint(1))/4)  ; % set lower limit of plot region%P add comments
        advancedvec(i).simFreqDeltaVec = (advancedvec(i).plot_lower:(advancedvec(i).plot_upper-advancedvec(i).plot_lower)/advanced.numinsimFreqDeltaVec:advancedvec(i).plot_upper);% P : EXPERIMENT
   else
    advancedvec(i).simFreqDeltaVec =  advancedvec(1).simFreqDeltaVec; 
   end     
    
    %%%%%
    [ fignum, mat ] = DelaySimOverview( fignum, basicvec(i), advancedvec(i), mat );  % Run simulations (if replications requested)
    if i==1 matvec = mat; else matvec = [matvec mat]; end
end

% do the plots based on the structures %P: i've tried to remove next loop
% for i=1:veclen
%     advancedvec(i).subtitle = sprintf('%s / %s', subtitle, legendvec{i});
%     fignum = DelayDoContours(fignum, basicvec(i), advancedvec(i), matvec(i));        % Generate a bunch of contour plots
% end
%[fignum] = UtilExperimentVectorPlot( fignum, basicvec, advancedvec, legendvec, matvec, subtitle,filemodifier);

figout = fignum;

end

