function [ figout ] = DelayPlotBayesInfo( fignum, basic, advanced, mat )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    if isfield(advanced,'subtitle')
        subtitle = sprintf(' (%s)', advanced.subtitle);
    else
        subtitle = '';
    end

    if advanced.saveplot % If files to be saved, write information about the version of code used.
        UtilSaveVersionFile( advanced.dirstring );
    end
    
    %%%%%% Mean number of samples to be observed
    fignum = fignum+1; figure(fignum); hold off;
    tauplusstageII = basic.tau + mat.ENumSamps;
    tauplusstageII(mat.ENumSamps == 0) = 0;
%    plot(mat.muvec,mat.optENumSamps,'-',mat.muvec,basic.TMax,'-.',mat.muvec,tauplusstageII,'.',[min(mat.bndlower) max(mat.bndupper)],[basic.tau basic.tau],'x')
    plot(mat.muvec,mat.optENumSamps,'-',mat.muvec,tauplusstageII,'.',[min(mat.bndlower) max(mat.bndupper)],[basic.tau basic.tau],'x')
%    UtilResizeFigureToBounds(mat);
    UtilStdizeFigure(fignum,advanced);
    title(sprintf('%s\nMean number of samples until stopping, %s',advanced.titlestring,subtitle),'FontSize',advanced.bigfontsize,'FontName',advanced.fontname);
    xlabel('Prior mean','FontSize',advanced.smallfontsize,'FontName',advanced.fontname);
    ylabel('E[num samples]','FontSize',advanced.smallfontsize,'FontName',advanced.fontname); 
%    legend('Sequential','All Tmax Samples','tau + Stage II samples');
    legend('Sequential (all stages)','tau + Stage II samples');
    if advanced.saveplot
        texttodifferentiate = 'ENSamp';
        UtilSaveFigFile(fignum, advanced.dirstring, advanced.filestring, texttodifferentiate, advanced.graphicextension);
    end
    
    %%%%%% reward plot
    fignum = fignum+1; figure(fignum); hold off;
    plot(mat.muvec,mat.RewardPIT0I,'-',mat.muvec,mat.BayesOneshotExpectedReward,'-.',mat.muvec,mat.B0hat','.',[min(mat.bndlower) max(mat.bndupper)],[basic.tau basic.tau],'x')
%    UtilResizeFigureToBounds(mat);
    UtilStdizeFigure(fignum,advanced);
    title(sprintf('%s\nExpected Reward, %s',advanced.titlestring,subtitle),'FontSize',advanced.bigfontsize,'FontName',advanced.fontname);
    xlabel('Prior mean','FontSize',advanced.smallfontsize,'FontName',advanced.fontname);
    ylabel('E[Reward | sampling scheme]','FontSize',advanced.smallfontsize,'FontName',advanced.fontname); 
    legend('with Perfect Info at t=0','All Tmax Samples','Sequential');
    if advanced.saveplot
        texttodifferentiate = 'EEOCSamp';
        UtilSaveFigFile(fignum, advanced.dirstring, advanced.filestring, texttodifferentiate, advanced.graphicextension);
    end

        %%%%%% regret plot
    fignum = fignum+1; figure(fignum); hold off;
    OneshotRegret = mat.RewardPIT0I(:) - mat.BayesOneshotExpectedReward(:);
    SeqRegret = mat.RewardPIT0I(:) - mat.B0hat(:);
    plot(mat.muvec(:),SeqRegret,'-',mat.muvec(:),OneshotRegret,'-.',[min(mat.bndlower) max(mat.bndupper)]',[basic.tau basic.tau]','x')
%    UtilResizeFigureToBounds(mat);
    UtilStdizeFigure(fignum,advanced);
    title(sprintf('%s\nBayes Posterior EOC, %s',advanced.titlestring,subtitle),'FontSize',advanced.bigfontsize,'FontName',advanced.fontname);
    xlabel('Prior mean','FontSize',advanced.smallfontsize,'FontName',advanced.fontname);
    ylabel('E[Regret | sampling scheme]','FontSize',advanced.smallfontsize,'FontName',advanced.fontname); 
    legend('Sequential','All Tmax Samples','Bound of Contin Set');
    if advanced.saveplot
        texttodifferentiate = 'EEOCSamp';
        UtilSaveFigFile(fignum, advanced.dirstring, advanced.filestring, texttodifferentiate, advanced.graphicextension);
    end

    %%%%%% Bayes Posterior prob correct selection
    fignum = fignum+1; figure(fignum); hold off;
    plot(mat.muvec,mat.Bayespcsvec,'-',mat.muvec,mat.OneshotPCS,'-.',[min(mat.bndlower) max(mat.bndupper)],[1 1],'x')
%    UtilResizeFigureToBounds(mat);
    UtilStdizeFigure(fignum,advanced);
    hold on;
    title(sprintf('%s\nBayes Posterior PCS, %s',advanced.titlestring,subtitle),'FontSize',advanced.bigfontsize,'FontName',advanced.fontname);
    xlabel('Prior mean','FontSize',advanced.smallfontsize,'FontName',advanced.fontname);
    ylabel('E[PCS | sampling scheme]','FontSize',advanced.smallfontsize,'FontName',advanced.fontname); 
    legend('Sequential','All Tmax Samples','Bound of Contin Set');
    if advanced.saveplot
        texttodifferentiate = 'EPCSSamp';
        UtilSaveFigFile(fignum, advanced.dirstring, advanced.filestring, texttodifferentiate, advanced.graphicextension);
    end

    
% mat.Puppermat = Puppermat;  % Probability of stopping due to hitting/exceeding upper boundary
% mat.Pnewmat = Pnewmat;      % Probability that if one stopped at a given (t, mu), that new technology is adopted
% mat.ENumSamps = ENumSampsin;
% mat.Regretcontmat = Regretcontmat;  % Expected regret of incorrect decision given that there is an opttion to continue sampling
% mat.Regretnowmat = Regretnowmat;  % Expected regret of incorrect decision if stopping at given (t,mu)
% mat.Bayespcsvec = pcsin;

figout = fignum;
    
end



