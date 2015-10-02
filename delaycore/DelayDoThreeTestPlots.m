function [ figout, mata, matb, matc, matd, mate, matf ] = DelayDoThreeTestPlots( fignum, basic, advanced )
%DELAYDOTHREETESTPLOTS 
%%%%
%%%% (c) 2014, Chick, Pertile, Forster
%%%%
%%%% For paper on optimal bayesian clinical trial stopping when there are
%%%% delays in information which arrives. This module can be called to
%%%% generate the following three plots automatically:
%%%%        1. Comparison of stopping boundaries with off line learning
%%%%        versus on line learning (outputs mata and matb)
%%%%        2. Comparison of stopping boundaries with a fixed number of
%%%%        patients to be treated when the decision to adopt or not is
%%%%        made, e.g. fixedP is true, versus what happens when the
%%%%        patients initially allocated for the trial (but not used in the
%%%%        trial due to early stopping) are given the 'better' treatment
%%%%        (outputs matc and matd)
%%%%        3. Comparison of stopping boundaries when the adoption decision
%%%%        can be changed tau time units after stopping, to allow for all
%%%%        data to be used in the stopping decision, versus a requirement
%%%%        that only patient data at the time of the stopping decision can
%%%%        be used to make the adoption decision ('no change'). The latter
%%%%        ignores the data from the final 'tau' patients whose data has
%%%%        not been finalized at the time of stopping the trial (outputs matf and matg)
%%%%
%%%% REQUIREMENTS:
%%%%    1. This is to be run as a function, with fignum the figure number
%%%%    of the prior figure to have been plotted, and with basic and
%%%%    advanced being validated input parameters for DelayCurvesRecur and
%%%%    DelayStageOne.
%%%%    2. It assumes that all clinical trial and numerical analysis
%%%%    parameters have been set, as in one of the 'Set*.m' files.
%%%%    3. It presumes that the variable 'fignum' has been set to the
%%%%    lowest number for the plots to be generated by this script.
%%%%
%%%% Last edits:
%   14 April 2014: Revamped from script into a functoin, and enabled this
%       to use the new parameter interface format (basic, advanced).
%   26 March 2014: Created, based on earlier version of DelayDriver.m.
%       Made modular so this script can be called from DelayDriver.m and
%       generate the graphs. Assumes that the base values of the parameters
%       for the clinical trial model have all been set up using the naming
%       conventions of the files Set*.m. Also assumes fignum variable has
%       been set up.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% COMPARE ON LINE LEARNING AND OFF LINE LEARNING BOUNDARY             %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%[basic, advanced] = SetBaseParameters();
%[basic, advanced, rval, msgs] = DelayInputValidator( basic, advanced );

% The setup
basic1 = basic; advanced1 = advanced; basic2 = basic; advanced2 = advanced;
basic1.online=false; basic2.online=true;
legend1='offline'; 
legend2='online';
advanced1.titlestring = sprintf('%s (%s)',advanced.titlestring,legend1);
advanced2.titlestring = sprintf('%s (%s)',advanced.titlestring,legend2);
subtitle='Online and offline learning';

if advanced.saveplot % If files to be saved, write information about the version of code used.
    UtilSaveVersionFile( advanced.dirstring );
end

% The standard machinery, part 1
[~, mat1] = DelayCurvesRecur(basic1,advanced1);
[mat1] = DelayStageOne(basic1,advanced1,mat1); 
fignum = fignum+1; figure(fignum); hold off;
hold off
plot(mat1.tvec,mat1.bndupper,'-o',mat1.tvec,mat1.bndlower,'-o');
hold on
if sum(mat1.bestsvec<=basic.tau)>0 
    plot(basic.t0+mat1.bestsvec(mat1.bestsvec<=basic.tau),mat1.muvec(mat1.bestsvec<=basic1.tau),'o');
end  % plot stopping region at time t0
[~, mat2] = DelayCurvesRecur(basic2,advanced2);
[mat2] = DelayStageOne(basic2,advanced2,mat2); 
figure(fignum);
plot(mat2.tvec,mat2.bndupper,'-+',mat2.tvec,mat2.bndlower,'-+');
if sum(mat2.bestsvec<=basic.tau)>0 
    plot(basic.t0+mat2.bestsvec(mat2.bestsvec<=basic.tau),mat2.muvec(mat2.bestsvec<=basic2.tau),'.');
end  % plot stopping region at time t0
legend(legend1,legend1,legend1,legend2,legend2,legend2);
title(sprintf('%s\n%s',advanced.titlestring,subtitle),'FontSize',advanced.bigfontsize,'FontName',advanced.fontname);
xlabel('t_0 + number of patients started','FontSize',advanced.smallfontsize,'FontName',advanced.fontname);
ylabel('posterior mean per patient','FontSize',advanced.smallfontsize,'FontName',advanced.fontname);
UtilResizeFigureToBounds([mat1 mat2]);       % reset so widest boundary is accounted for
UtilStdizeFigure(fignum,advanced);
if advanced.saveplot
    texttodifferentiate = 'Do3OnOffa';
    UtilSaveFigFile(fignum, advanced.dirstring, advanced.filestring, texttodifferentiate, advanced.graphicextension);
end

% The standard machinery, part 2
% compare the optimal number of samples for stage I
fignum=fignum+1; figure(fignum); hold off;
plot(mat1.muvec,mat1.bestsvec,'o',mat2.muvec,mat2.bestsvec,'.');
title(sprintf('%s\n%s',advanced.titlestring,subtitle),'FontSize',advanced.bigfontsize,'FontName',advanced.fontname)
xlabel('t_0 + number of patients started','FontSize',advanced.smallfontsize,'FontName',advanced.fontname);
ylabel('optimal num samples, stage I','FontSize',advanced.smallfontsize,'FontName',advanced.fontname);
legend(legend1, legend2);
UtilStdizeFigure(fignum,advanced);
if advanced.saveplot
    texttodifferentiate = 'Do3OnOffb';
    UtilSaveFigFile(fignum, advanced.dirstring, advanced.filestring, texttodifferentiate, advanced.graphicextension);
end

% The standard machinery, part 3
%[ fignum ] = DelayDoContours( fignum, basic1, advanced1, mat1 );
%[ fignum ] = DelayDoContours( fignum, basic2, advanced2, mat2 );

% custom assignments
mata = mat1; matb = mat2;
%advanceda = advanced1; advancedb = advanced2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPARE FIXED P VERSUS CAN TREAT OTHER PATIENTS NOT DONE IN TRIAL WITH  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The setup
basic1 = basic; advanced1 = advanced; basic2 = basic; advanced2 = advanced;
advanced1.fixedP=true; advanced2.fixedP=false;
legend1='fixed P'; 
legend2='early switch';
advanced1.titlestring = sprintf('%s (%s)',advanced.titlestring,legend1);
advanced2.titlestring = sprintf('%s (%s)',advanced.titlestring,legend2);
subtitle=sprintf('Fixed P v. Early stop allows switch (online = %d)',basic.online);

% The standard machinery, part 1
[~, mat1] = DelayCurvesRecur(basic1,advanced1);
[mat1] = DelayStageOne(basic1,advanced1,mat1); 
fignum = fignum+1; figure(fignum); hold off;
hold off
plot(mat1.tvec,mat1.bndupper,'-o',mat1.tvec,mat1.bndlower,'-o');
hold on
if sum(mat1.bestsvec<=basic.tau)>0 
    plot(basic.t0+mat1.bestsvec(mat1.bestsvec<=basic.tau),mat1.muvec(mat1.bestsvec<=basic1.tau),'o');
end  % plot stopping region at time t0
[~, mat2] = DelayCurvesRecur(basic2,advanced2);
[mat2] = DelayStageOne(basic2,advanced2,mat2); 
figure(fignum);
plot(mat2.tvec,mat2.bndupper,'-+',mat2.tvec,mat2.bndlower,'-+');
if sum(mat2.bestsvec<=basic.tau)>0 
    plot(basic.t0+mat2.bestsvec(mat2.bestsvec<=basic.tau),mat2.muvec(mat2.bestsvec<=basic2.tau),'.');
end  % plot stopping region at time t0
legend(legend1,legend1,legend1,legend2,legend2,legend2);
title(sprintf('%s\n%s',advanced.titlestring,subtitle),'FontSize',advanced.bigfontsize,'FontName',advanced.fontname);
xlabel('t_0 + number of patients started','FontSize',advanced.smallfontsize,'FontName',advanced.fontname);
ylabel('posterior mean per patient','FontSize',advanced.smallfontsize,'FontName',advanced.fontname);
UtilResizeFigureToBounds([mat1 mat2]);       % reset so widest boundary is accounted for
UtilStdizeFigure(fignum,advanced);
if advanced.saveplot
    texttodifferentiate = 'Do3fixedPa';
    UtilSaveFigFile(fignum, advanced.dirstring, advanced.filestring, texttodifferentiate, advanced.graphicextension);
end

% The standard machinery, part 2
% compare the optimal number of samples for stage I
fignum=fignum+1; figure(fignum); hold off;
plot(mat1.muvec,mat1.bestsvec,'o',mat2.muvec,mat2.bestsvec,'.');
title(sprintf('%s\n%s',advanced.titlestring,subtitle),'FontSize',advanced.bigfontsize,'FontName',advanced.fontname)
xlabel('posterior mean per patient','FontSize',advanced.smallfontsize,'FontName',advanced.fontname);
ylabel('optimal num samples in stage I','FontSize',advanced.smallfontsize,'FontName',advanced.fontname);
legend(legend1, legend2);
UtilStdizeFigure(fignum,advanced);
if advanced.saveplot
    texttodifferentiate = 'Do3fixedPb';
    UtilSaveFigFile(fignum, advanced.dirstring, advanced.filestring, texttodifferentiate, advanced.graphicextension);
end

% The standard machinery, part 3
%[ fignum ] = DelayDoContours( fignum, basic1, advanced1, mat1 );
%[ fignum ] = DelayDoContours( fignum, basic2, advanced2, mat2 );

% custom assignments
matc = mat1; matd = mat2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPARE base case (change ok after remaining data arrives) vs no change %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The setup
basic1 = basic; advanced1 = advanced; basic2 = basic; advanced2 = advanced;
advanced1.nochangetest=true; advanced2.nochangetest=false;
legend1='pick after data arrives'; 
legend2='pick immediately on stop';
advanced1.titlestring = sprintf('%s (%s)',advanced.titlestring,legend1);
advanced2.titlestring = sprintf('%s (%s)',advanced.titlestring,legend2);
subtitle=sprintf('Pick best after new data arrives or immediately upon stopping (online = %d)',basic.online);

% The standard machinery, part 1
[~, mat1] = DelayCurvesRecur(basic1,advanced1);
[mat1] = DelayStageOne(basic1,advanced1,mat1); 
fignum = fignum+1; figure(fignum); hold off;
plot(mat1.tvec,mat1.bndupper,'-o',mat1.tvec,mat1.bndlower,'-o');
hold on
if sum(mat1.bestsvec<=basic.tau)>0 
    plot(basic.t0+mat1.bestsvec(mat1.bestsvec<=basic.tau),mat1.muvec(mat1.bestsvec<=basic1.tau),'o');
end  % plot stopping region at time t0
[~, mat2] = DelayCurvesRecur(basic2,advanced2);
[mat2] = DelayStageOne(basic2,advanced2,mat2); 
figure(fignum);
plot(mat2.tvec,mat2.bndupper,'-+',mat2.tvec,mat2.bndlower,'-+');
if sum(mat2.bestsvec<=basic.tau)>0 
    plot(basic.t0+mat2.bestsvec(mat2.bestsvec<=basic.tau),mat2.muvec(mat2.bestsvec<=basic2.tau),'.');
end  % plot stopping region at time t0
legend(legend1,legend1,legend1,legend2,legend2,legend2);
title(sprintf('%s\n%s',advanced.titlestring,subtitle),'FontSize',advanced.bigfontsize,'FontName',advanced.fontname);
xlabel('t_0 + number of patients started','FontSize',advanced.smallfontsize,'FontName',advanced.fontname);
ylabel('posterior mean per patient','FontSize',advanced.smallfontsize,'FontName',advanced.fontname);
UtilResizeFigureToBounds([mat1 mat2]);       % reset so widest boundary is accounted for
UtilStdizeFigure(fignum,advanced);
if advanced.saveplot
    texttodifferentiate = 'Do3changeOKa';
    UtilSaveFigFile(fignum, advanced.dirstring, advanced.filestring, texttodifferentiate, advanced.graphicextension);
end

% The standard machinery, part 2
% compare the optimal number of samples for stage I
fignum=fignum+1; figure(fignum); hold off;
plot(mat1.muvec,mat1.bestsvec,'o',mat2.muvec,mat2.bestsvec,'.');
title(sprintf('%s\n%s',advanced.titlestring,subtitle),'FontSize',advanced.bigfontsize,'FontName',advanced.fontname)
xlabel('posterior mean per patient','FontSize',advanced.smallfontsize,'FontName',advanced.fontname);
ylabel('optimal num samples in stage I','FontSize',advanced.smallfontsize,'FontName',advanced.fontname);
legend(legend1, legend2);
UtilStdizeFigure(fignum,advanced);
if advanced.saveplot
    texttodifferentiate = 'Do3changeOKb';
    UtilSaveFigFile(fignum, advanced.dirstring, advanced.filestring, texttodifferentiate, advanced.graphicextension);
end

% The standard machinery, part 3
%[ fignum ] = DelayDoContours( fignum, basic1, advanced1, mat1 );
%[ fignum ] = DelayDoContours( fignum, basic2, advanced2, mat2 );

% custom assignments
mate = mat1; matf = mat2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figout = fignum;
