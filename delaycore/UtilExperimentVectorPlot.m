function [ figout ] = UtilExperimentVectorPlot( fignum, basicvec, advancedvec, legendvec, matvec, subtitle, fnamediff )
%UTILEXPERIMENTVECTORPLOT 
%
% Given a vector of basic and advanced parameters, and a vector of the
% output analysis in matvec, generate some plots which use the legend
% vector. Useful when trying to systematically vary a given field over a
% range of values in fieldvec.
%%%%
%%%% REQUIREMENTS:
%%%%    1. Assumes input vectors have been validated and the stage II and
%%%%    stage I analyses have been completed, etc.
%%%%    2. See DelayDriver for an example usage case.
%%%% INPUTS:
%%%%    basicvec, advancedvec, legendvec, matvec: created by above routines
%%%%    subtitle: used in plot title in second line to describe plot
%%%%    fnamediff: used in file name of plot to help differentiate name of
%%%%    file, if a file is to be written with the plot
%%%%
%%%% Last edits:
%   17 April 2014: Written.
%   18 April 2014: added fnamediff
%
%   For 'delay' sequential trial paper of Paolo Martin and Steve
%
% coded by steve.

veclen = length(matvec);
mycode = {'-o', '-+', '-x', '-s', '-d', '-^', '-v'} ;
mycode2 ={ 'o',  '.',  'x',  's',  'd',  '^',  'v'} ;

if veclen <= 1
    'UtilExperimentVectorPlot.m called with fewer than 2 experiments in vector. Normally should be called with at least 2, based on output from TestDelayIterate'
%    warning('UtilExperimentVectorPlot.m called with fewer than 2 experiments in vector. Normally should be called with at least 2, based on output from TestDelayIterate');
end
% The standard machinery, part 1a
% plot the upper and lower stopping boundaries
fignum = fignum+1; figure(fignum); hold off;
hold off
for i=1:veclen
    plot(matvec(i).tvec,matvec(i).bndupper,mycode{1+mod(i-1,length(mycode))});
    hold on
end
if length(legendvec) > 1 
    legend(legendvec,'Location','SouthEast'); 
end;
for i=1:veclen
    plot(matvec(i).tvec,matvec(i).bndlower,mycode{1+mod(i-1,length(mycode))});
    if sum(matvec(i).bestsvec<=basicvec(i).tau)>0, plot(basicvec(i).t0+matvec(i).bestsvec(matvec(i).bestsvec<=basicvec(i).tau),matvec(i).muvec(matvec(i).bestsvec<=basicvec(i).tau),mycode2{1+mod(i-1,length(mycode2))});end  % plot stopping region at time t0
    if advancedvec(i).saveplot
        UtilSaveVersionFile( advancedvec(i).dirstring );
    end
end
%legend(legend1,legend1,legend1,legend2,legend2,legend2);
i=1;
title(sprintf('%s\n%s',advancedvec(i).titlestring,subtitle),'FontSize',advancedvec(i).bigfontsize,'FontName',advancedvec(i).fontname);
xlabel('t_0 + number of patients started','FontSize',advancedvec(i).smallfontsize,'FontName',advancedvec(i).fontname);
ylabel('posterior mean per patient','FontSize',advancedvec(i).smallfontsize,'FontName',advancedvec(i).fontname);
UtilResizeFigureToBounds( matvec );
UtilStdizeFigure(fignum,advancedvec(i));
if advancedvec(i).saveplot
    texttodifferentiate = sprintf('ExpVec1%s',fnamediff);
    UtilSaveFigFile(fignum, advancedvec(i).dirstring, advancedvec(i).filestring, texttodifferentiate, advancedvec(i).graphicextension);
end

% The standard machinery, part 1b
% plot the upper and lower stopping boundaries
allsametau = 1; keeptau=basicvec(1).tau;
for i=2:veclen          % check if all taus for all experiments are the same
    allsametau = allsametau && (keeptau == basicvec(i).tau);
end
% if not all tau are the same, then plot an extra graph, which displays the number
% of samples seen up to the time of stopping
fignum = fignum+1; figure(fignum); hold off;
hold off
for i=1:veclen
    plot(matvec(i).tvec-basicvec(i).tau,matvec(i).bndupper,mycode{1+mod(i-1,length(mycode))});
    hold on
end
if length(legendvec) > 1 
    legend(legendvec,'Location','SouthEast'); 
end;
for i=1:veclen
    plot(matvec(i).tvec-basicvec(i).tau,matvec(i).bndlower,mycode{1+mod(i-1,length(mycode))});
%    if sum(matvec(i).bestsvec<=basicvec(i).tau)>0, plot(basicvec(i).t0+matvec(i).bestsvec(matvec(i).bestsvec<=basicvec(i).tau),matvec(i).muvec(matvec(i).bestsvec<=basicvec(i).tau),mycode2{1+mod(i-1,length(mycode2))});end  % plot stopping region at time t0
end
%legend(legend1,legend1,legend1,legend2,legend2,legend2);
i=1;
title(sprintf('%s\n%s',advancedvec(i).titlestring,subtitle),'FontSize',advancedvec(i).bigfontsize,'FontName',advancedvec(i).fontname);
xlabel('t0 + number of samples seen upon stopping','FontSize',advancedvec(i).smallfontsize,'FontName',advancedvec(i).fontname);
ylabel('posterior mean per patient','FontSize',advancedvec(i).smallfontsize,'FontName',advancedvec(i).fontname);
UtilResizeFigureToBounds( matvec );
UtilStdizeFigure(fignum,advancedvec(i));
if advancedvec(i).saveplot
    texttodifferentiate = sprintf('ExpVec1b%s',fnamediff);
    UtilSaveFigFile(fignum, advancedvec(i).dirstring, advancedvec(i).filestring, texttodifferentiate, advancedvec(i).graphicextension);
end

% The standard machinery, part 2
% compare the optimal number of samples for stage I
fignum=fignum+1; figure(fignum); hold off;
hold off
for i=1:veclen
    plot(matvec(i).muvec,matvec(i).bestsvec,mycode2{1+mod(i-1,length(mycode2))});
    hold on
end
i=1;
title(sprintf('%s\n%s',advancedvec(i).titlestring,subtitle),'FontSize',advancedvec(i).bigfontsize,'FontName',advancedvec(i).fontname)
xlabel('t_0 + number of patients started','FontSize',advancedvec(i).smallfontsize,'FontName',advancedvec(i).fontname);
ylabel('optimal num samples, stage I','FontSize',advancedvec(i).smallfontsize,'FontName',advancedvec(i).fontname);
if length(legendvec) > 1 
    legend(legendvec,'Location','SouthEast'); 
end;
UtilStdizeFigure(fignum,advancedvec(i));
if advancedvec(i).saveplot
    texttodifferentiate = sprintf('ExpVec2%s',fnamediff);
    UtilSaveFigFile(fignum, advancedvec(i).dirstring, advancedvec(i).filestring, texttodifferentiate, advancedvec(i).graphicextension);
end

%%%%%% Mean number of samples to be observed
fignum = fignum+1; figure(fignum);  hold off;
hold off;
for i=1:veclen
    plot(matvec(i).muvec,matvec(i).optENumSamps,mycode{1+mod(i-1,length(mycode))});
    hold on;
end
if length(legendvec) > 1 
    legend(legendvec,'Location','SouthEast'); 
end;
for i=1:veclen
    plot(matvec(i).muvec,basicvec(i).TMax,mycode2{1+mod(i-1,length(mycode2))},[min(matvec(i).bndlower) max(matvec(i).bndupper)],[basicvec(i).tau basicvec(i).tau],mycode2{1+mod(i-1,length(mycode2))});
end
i=1;
%    UtilResizeFigureToBounds(mat);
UtilStdizeFigure(fignum,advancedvec(i));
title(sprintf('%s\nE[Num. samples until stopping], %s',advancedvec(i).titlestring,subtitle),'FontSize',advancedvec(i).bigfontsize,'FontName',advancedvec(i).fontname);
xlabel('Prior mean','FontSize',advancedvec(i).smallfontsize,'FontName',advancedvec(i).fontname);
ylabel('E[num samples]','FontSize',advancedvec(i).smallfontsize,'FontName',advancedvec(i).fontname); 
if advancedvec(i).saveplot
    texttodifferentiate = 'MultiENSamp';
    UtilSaveFigFile(fignum, advancedvec(i).dirstring, advancedvec(i).filestring, texttodifferentiate, advancedvec(i).graphicextension);
end

%%%%%%%%%%%%%%%% Plot expected regret for sequential
fignum = fignum+1; figure(fignum); hold off;
hold off;
for i=1:veclen
    SeqRegret = matvec(i).RewardPIT0I(:) - matvec(i).B0hat(:);
    plot(matvec(i).muvec,SeqRegret,mycode{1+mod(i-1,length(mycode))});
    hold on;
end
if length(legendvec) > 1 
    legend(legendvec,'Location','SouthEast'); 
end;
i=1;
%    UtilResizeFigureToBounds(mat);
UtilStdizeFigure(fignum,advancedvec(i));
title(sprintf('%s\nE[Regret] for Sequential, %s',advancedvec(i).titlestring,subtitle),'FontSize',advancedvec(i).bigfontsize,'FontName',advancedvec(i).fontname);
xlabel('Prior mean','FontSize',advancedvec(i).smallfontsize,'FontName',advancedvec(i).fontname);
ylabel('E[Regret | sampling scheme]','FontSize',advancedvec(i).smallfontsize,'FontName',advancedvec(i).fontname); 
if advancedvec(i).saveplot
    texttodifferentiate = 'MultiRegSeq';
    UtilSaveFigFile(fignum, advancedvec(i).dirstring, advancedvec(i).filestring, texttodifferentiate, advancedvec(i).graphicextension);
end

%%%%%%%%%%%%%%%% Plot expected regret for sequential and one shot runs
fignum = fignum+1; figure(fignum); hold off;
hold off;
for i=1:veclen
    OneshotRegret = matvec(i).RewardPIT0I(:) - matvec(i).BayesOneshotExpectedReward(:);
    plot(matvec(i).muvec,OneshotRegret,mycode2{1+mod(i-1,length(mycode2))},[min(matvec(i).bndlower) max(matvec(i).bndupper)],[0.0 0.0],mycode2{1+mod(i-1,length(mycode2))});
    hold on;
end
if length(legendvec) > 1 
    legend(legendvec,'Location','SouthEast'); 
end;
i=1;
%    UtilResizeFigureToBounds(mat);
UtilStdizeFigure(fignum,advancedvec(i));
title(sprintf('%s\nE[Regret] for One Shot, %s',advancedvec(i).titlestring,subtitle),'FontSize',advancedvec(i).bigfontsize,'FontName',advancedvec(i).fontname);
xlabel('Prior mean','FontSize',advancedvec(i).smallfontsize,'FontName',advancedvec(i).fontname);
ylabel('E[Regret | sampling scheme]','FontSize',advancedvec(i).smallfontsize,'FontName',advancedvec(i).fontname); 
if advancedvec(i).saveplot
    texttodifferentiate = 'MultiRegOneStep';
    UtilSaveFigFile(fignum, advancedvec(i).dirstring, advancedvec(i).filestring, texttodifferentiate, advancedvec(i).graphicextension);
end

%%%%%% Bayes Posterior prob correct selection
fignum = fignum+1; figure(fignum);  hold off;
hold off;
for i=1:veclen
    plot(matvec(i).muvec,matvec(i).Bayespcsvec,mycode{1+mod(i-1,length(mycode))});
    hold on;
end
if length(legendvec) > 1 
    legend(legendvec,'Location','SouthEast'); 
end;
for i=1:veclen
    plot(matvec(i).muvec,matvec(i).OneshotPCS,mycode2{1+mod(i-1,length(mycode2))},[min(matvec(i).bndlower) max(matvec(i).bndupper)],[0.0 0.0],mycode2{1+mod(i-1,length(mycode2))});
end
i=1;
%    UtilResizeFigureToBounds(mat);
UtilStdizeFigure(fignum,advancedvec(i));
title(sprintf('%s\nE[PCS] for Sequential and One Shot, %s',advancedvec(i).titlestring,subtitle),'FontSize',advancedvec(i).bigfontsize,'FontName',advancedvec(i).fontname);
xlabel('Prior mean','FontSize',advancedvec(i).smallfontsize,'FontName',advancedvec(i).fontname);
ylabel('E[PCS | sampling scheme]','FontSize',advancedvec(i).smallfontsize,'FontName',advancedvec(i).fontname); 
if advancedvec(i).saveplot
    texttodifferentiate = 'MultiEPCSSamp';
    UtilSaveFigFile(fignum, advancedvec(i).dirstring, advancedvec(i).filestring, texttodifferentiate, advancedvec(i).graphicextension);
end

%%%%%% PDE computation for B0hat
fignum = fignum+1; figure(fignum); hold off;
hold off;
for i=1:veclen
    plot(matvec(i).muvec,matvec(i).B0hat,mycode{1+mod(i-1,length(mycode))});
    hold on;
end
if length(legendvec) > 1 
    legend(legendvec,'Location','SouthEast'); 
end;
for i=1:veclen
    plot([min(matvec(i).bndlower) max(matvec(i).bndupper)],[0.0 0.0],mycode2{1+mod(i-1,length(mycode2))});
end
i=1;
%    UtilResizeFigureToBounds(mat);
UtilStdizeFigure(fignum,advancedvec(i));
title(sprintf('%s\nB0hat = E[Value | at start of stage I], %s',advancedvec(i).titlestring,subtitle),'FontSize',advancedvec(i).bigfontsize,'FontName',advancedvec(i).fontname);
xlabel('Prior mean','FontSize',advancedvec(i).smallfontsize,'FontName',advancedvec(i).fontname);
ylabel('E[reward | both stages]','FontSize',advancedvec(i).smallfontsize,'FontName',advancedvec(i).fontname); 
if advancedvec(i).saveplot
    texttodifferentiate = 'B0hat';
    UtilSaveFigFile(fignum, advancedvec(i).dirstring, advancedvec(i).filestring, texttodifferentiate, advancedvec(i).graphicextension);
end

%%%%%% Bayes Posterior B0vec (stage two only)
fignum = fignum+1; figure(fignum);  hold off;
hold off;
for i=1:veclen
    plot(matvec(i).muvec,matvec(i).B0vec(:),mycode{1+mod(i-1,length(mycode))});
    hold on;
end
if length(legendvec) > 1 
    legend(legendvec,'Location','SouthEast'); 
end;
for i=1:veclen
    plot([min(matvec(i).bndlower) max(matvec(i).bndupper)],[0.0 0.0],mycode2{1+mod(i-1,length(mycode2))});
end
i=1;
%    UtilResizeFigureToBounds(mat);
UtilStdizeFigure(fignum,advancedvec(i));
title(sprintf('%s\nB0, = E[Value, at start of stage II], %s',advancedvec(i).titlestring,subtitle),'FontSize',advancedvec(i).bigfontsize,'FontName',advancedvec(i).fontname);
xlabel('Prior mean','FontSize',advancedvec(i).smallfontsize,'FontName',advancedvec(i).fontname);
ylabel('E[reward | start stage II]','FontSize',advancedvec(i).smallfontsize,'FontName',advancedvec(i).fontname); 
if advancedvec(i).saveplot
    texttodifferentiate = 'B0vec';
    UtilSaveFigFile(fignum, advancedvec(i).dirstring, advancedvec(i).filestring, texttodifferentiate, advancedvec(i).graphicextension);
end

%%%%% PUT SOME GRAPHS OF SOME SIMULATIONS!
fields=fieldnames(matvec(1).outBayes);
codePDE = '-.'; %codePDEFr = ':';
for i=1:numel(fields)
    a = strfind( fields{i}, 'Ave' );
    if ~isempty(a)          % if 'Ave' is in name, then do some plots
        basename = fields{i}(1:(a-1));
        fignum=fignum+1; figure(fignum); hold off;
        for j=1:veclen
            avevec = matvec(j).outBayes.(fields{i});
            plot(advancedvec(j).simFreqDeltaVec,avevec,mycode{1+mod(j-1,length(mycode))},'LineWidth',2);
            hold on
        end
        if length(legendvec) > 1 
            legend(legendvec,'Location','SouthEast'); 
        end;
        for j=1:veclen
            avevec = matvec(j).outBayes.(fields{i});
            sevec = matvec(j).outBayes.(strcat(basename,'SE'));
            avevecfreq = matvec(j).outFreq.(fields{i});
            sevecfreq = matvec(j).outFreq.(strcat(basename,'SE'));
            % plot frequentist result
            plot(advancedvec(j).simFreqDeltaVec,avevecfreq,mycode2{1+mod(j-1,length(mycode2))},'LineWidth',2);
            % now plot error bounds for bayes and frequentist results
            plot(advancedvec(j).simFreqDeltaVec,avevec+sevec,'--',advancedvec(j).simFreqDeltaVec,avevec-sevec,'--','LineWidth',0.5);
            plot(advancedvec(j).simFreqDeltaVec,avevecfreq+sevecfreq,'-.',advancedvec(j).simFreqDeltaVec,avevecfreq-sevecfreq,'-.','LineWidth',0.5);
            if advancedvec(j).PDEWithSimPlots
                switch basename     % add some theoretical results if desired
                  case 'ENumSeen'   %% mean number of samples as function of prior mean (works).
                    plot(matvec(j).muvec(:),matvec(j).optENumSamps(:),codePDE);
        %                  case 'DecisionMean'      % has slope of less than 45 degrees. Regret?
        %                  case 'DecisionMeanAll' % has slope of 45 degrees
        %                  case 'RealizedMean'
        %                  case 'RealizedReward'
        %                  case 'RealizedRewardAll'
                  case 'ExpectedReward'        % FIX: NOTE, PDE solution is 'higher' than the simulaiton solution. Discounting? Or because simulation is suboptimal to PDE when mean is high?
                    plot(matvec(j).muvec(:),matvec(j).B0hat(:),codePDE);
        %                plot(mat.muvec(:),mat.BayesOneshotExpectedReward(:),codePDEFr);
                  case 'PrHitTop'   % Prob (top stopping boundary is hit before bottom boundary)(works).
                    plot(matvec(j).muvec(:),matvec(j).Puppermat(:,1),codePDE);
                  case 'PrNewSelected'          % Prob (new alternative is selected) (works).
                    plot(matvec(j).muvec(:),matvec(j).Pnewmat(:,1),codePDE);
        %                  case 'PrNewSelectedAll'
                  case 'PrTrueBestSelected'     % PCS, given sequential procedure 
                    plot(matvec(j).muvec(:),matvec(j).Bayespcsvec(:),codePDE);
    %                    plot(matvec(j).muvec(:),matvec(j).OneshotPCS(:),codePDEFr);
                end % switch
            end % PDEWithSimPlots
        end % for
        
        % do some formatting and save if necessary
        UtilStdizeFigure( fignum, advancedvec(j) );
        title(sprintf('%s\n%s (sample mean +/- std err)',advancedvec(j).titlestring,basename),'FontSize',advancedvec(j).bigfontsize,'FontName',advancedvec(j).fontname);
        xlabel('Mean reward (prior or true)','FontSize',advancedvec(j).smallfontsize,'FontName',advancedvec(j).fontname); 
        ylabel(basename,'FontSize',advancedvec(j).smallfontsize,'FontName',advancedvec(j).fontname);
        if strncmp('Pr',fields{i},2)   % if field is for a probability, set y-axis to be 0-1
            plot([min(matvec(j).bndlower) max(matvec(j).bndupper)],[1 1],mycode2{1+mod(j-1,length(mycode2))});
            if strncmp('PrTrueBest',fields{i},10)
                tmp = axis; tmp(3)=0.4; tmp(4)=1.01; axis(tmp);
            else
                tmp=axis; tmp(3)=0.0; tmp(4)=1.01; axis(tmp);
            end
        else
            plot([min(matvec(j).bndlower) max(matvec(j).bndupper)],(basicvec(j).t0+basicvec(j).tau)*[0 0],mycode2{1+mod(j-1,length(mycode2))});
        end
        %axis('square');
        if advancedvec(j).saveplot
            texttodifferentiate = sprintf('Sim%sN%d',basename,i);
            UtilSaveFigFile(fignum, advancedvec(j).dirstring, advancedvec(j).filestring, texttodifferentiate, advancedvec(j).graphicextension);
        end
    end %~isempty(a)
end % for i = 1:numfields

%NB: Don't plot the contours - that can be done separately, otherwise it
%creates a lot of graphs.
% The standard machinery, part 3
%for i = 1:veclen
%    [ fignum ] = DelayDoContours( fignum, basicvec(i), advancedvec(i), matvec(i) );
%end

figout = fignum;
