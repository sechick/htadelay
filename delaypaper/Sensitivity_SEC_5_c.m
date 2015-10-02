%% MACRO TO BE CALLED FROM SEC-Test-Experiments
%
% Generates a bunch of plots given vector of analysis created in the above
% file. Useful for plotting graphs when tau is varied, e.g. for a
% presentation which displays how anaylsis might be sensitive to the delay.
%
% SEC: 7 april 2015

 
veclen = length(matvec);
mycode = {'-', '--', '-.', '-s', '-d', '-^', '-v'} ;
mycode2 ={ 'o',  'x',  '^',  's',  'd',  '.',  'v'} ;
%mycode = {'-o', '-+', '-x', '-s', '-d', '-^', '-v'} ;
%mycode2 ={ 'o',  '.',  'x',  's',  'd',  '^',  'v'} ;
mycolor ={ 'k',  'k',  'k',  'k',  'k',  'k',  'k'} ;
linewid = 2;
legendvec = {'cost of sampling: 200','cost of sampling: 5000'} ; %


%%%%%%%%%%%%%%%% Plot bayes versus trial
fignum = fignum+1; figure(fignum); hold off;
topmax = 0;
for i=1:veclen
    advancedvec(i).dirstring = 'Figure';
    advancedvec(i).graphicextension = 'eps';
    DiffTrial = matvec(i).outBayes.RealizedRewardAve(:) - matvec(i).outBayes.RealizedRewardTrialAve(:);
    topmax = max(topmax,max(DiffTrial));
    DiffAll = matvec(i).outBayes.RealizedRewardAve(:) - matvec(i).outBayes.RealizedRewardAllAve(:);
    topmax = max(topmax,max(DiffAll));
    DiffOneShot = matvec(i).outBayes.RealizedRewardAve(:) - matvec(i).outBayes.RealizedRewardFixedAve(:);
    topmax = max(topmax,max(DiffOneShot));
    plot(advancedvec(i).simFreqDeltaVec,DiffTrial, strcat(mycode{1+mod(i-1,length(mycode))},mycolor{1+mod(i-1,length(mycolor))} ) ,'LineWidth',linewid);
%    plot(advancedvec(i).simFreqDeltaVec,DiffTrial,strcat(mycode{1+mod(i-1,length(mycode))},mycolor{1+mod(i-1,length(mycolor))} ),'LineWidth',linewid);
%    plot(matvec(i).muvec,DiffOneShot,strcat(mycode{1+mod(i-1,length(mycode))},mycolor{1+mod(i-1,length(mycolor))} ),'LineWidth',linewid);
    hold on;
end
if length(legendvec) > 1 
         legend(legendvec,'Location','NorthWest'); 
end;
i=1;
%    UtilResizeFigureToBounds(mat);
tmp=axis();
tmp(1)=MULOW;
tmp(2)=MUHIGH;
tmp(3)= -VALDIFMAX/10;
tmp(4)=50000000;
axis(tmp);
UtilStdizeFigure(fignum,advancedvec(i));set(gca,'fontsize',advancedvec(i).bigfontsize);
%title(sprintf('%s\nSequential Bayes vs. Trial, %s',advancedvec(i).titlestring,subtitle),'FontSize',advancedvec(i).bigfontsize,'FontName',advancedvec(i).fontname);
xlabel('Prior mean','FontSize',advancedvec(i).bigfontsize,'FontName',advancedvec(i).fontname);
ylabel('Reward: Optimal Bayes vs. Trial','FontSize',advancedvec(i).bigfontsize,'FontName',advancedvec(i).fontname); 
if advancedvec(i).saveplot
    texttodifferentiate = 'OptVsTrial';
%    UtilSaveFigFile(fignum, advancedvec(i).dirstring, advancedvec(i).filestring, texttodifferentiate, advancedvec(i).graphicextension);
    UtilSaveFigEpsPdf(fignum,advancedvec(i).dirstring,strcat(advancedvec(i).dirstring,texttodifferentiate),'-r600');
end

%%%%%%%%%%%%%%%% Plot expected reward optimal sequential versus Tmax
% fignum = fignum+1; figure(fignum); hold off;
% hold off;
% for i=1:veclen
% %    DiffOneShot = matvec(i).B0hat(:) - matvec(i).OptOneShotReward(:);
%     DiffAll = matvec(i).outBayes.RealizedRewardAve(:) - matvec(i).outBayes.RealizedRewardAllAve(:);
%     plot(advancedvec(i).simFreqDeltaVec,DiffAll,strcat(mycode{1+mod(i-1,length(mycode))},mycolor{1+mod(i-1,length(mycolor))} ),'LineWidth',linewid);
% %    plot(matvec(i).muvec,DiffOneShot,strcat(mycode{1+mod(i-1,length(mycode))},mycolor{1+mod(i-1,length(mycolor))} ),'LineWidth',linewid);
%     hold on;
% end
% if length(legendvec) > 1 
%     legend(legendvec,'Location','best','FontSize',12); 
% end;
% i=1;
% %    UtilResizeFigureToBounds(mat);
% tmp=axis();
% tmp(1)=MULOW;
% tmp(2)=MUHIGH;
% tmp(3)= -VALDIFMAX/10;
% tmp(4)=VALDIFMAX;
% axis(tmp);
% UtilStdizeFigure(fignum,advancedvec(i));set(gca,'fontsize',advancedvec(i).bigfontsize);
% %title(sprintf('%s\nSequential Bayes vs. Tmax, %s',advancedvec(i).titlestring,subtitle),'FontSize',advancedvec(i).bigfontsize,'FontName',advancedvec(i).fontname);
% xlabel('Prior mean','FontSize',advancedvec(i).bigfontsize,'FontName',advancedvec(i).fontname);
% ylabel('Reward: Optimal Bayes vs. Tmax','FontSize',advancedvec(i).bigfontsize,'FontName',advancedvec(i).fontname); 
% if advancedvec(i).saveplot
%     texttodifferentiate = 'OptVsTmax';
% %    UtilSaveFigFile(fignum, advancedvec(i).dirstring, advancedvec(i).filestring, texttodifferentiate, advancedvec(i).graphicextension);
%     UtilSaveFigEpsPdf(fignum,advancedvec(i).dirstring,strcat(advancedvec(i).dirstring,texttodifferentiate),'-r600');
% end

%%%%%%%%%%%%%%%% Plot expected reward optimal sequential versus optimal one
%%%%%%%%%%%%%%%% shot
fignum = fignum+1; figure(fignum); hold off;
hold off;
for i=1:veclen
%    DiffOneShot = matvec(i).B0hat(:) - matvec(i).OptOneShotReward(:);
    DiffOneShot = matvec(i).outBayes.RealizedRewardAve(:) - matvec(i).outBayes.RealizedRewardFixedAve(:);
    plot(advancedvec(i).simFreqDeltaVec,DiffOneShot,strcat(mycode{1+mod(i-1,length(mycode))},mycolor{1+mod(i-1,length(mycolor))} ),'LineWidth',linewid);
%    plot(matvec(i).muvec,DiffOneShot,strcat(mycode{1+mod(i-1,length(mycode))},mycolor{1+mod(i-1,length(mycolor))} ),'LineWidth',linewid);
    hold on;
end
if length(legendvec) > 1 
    legend(legendvec,'Location','NorthWest'); 
end;
i=1;
%    UtilResizeFigureToBounds(mat);
tmp=axis();
tmp(1)=MULOW;
tmp(2)=MUHIGH;
if sameheightforallmaxloss 
    tmp(3)= -VALDIFMAX/10;
    tmp(4)=VALDIFMAX;
end
axis(tmp);
UtilStdizeFigure(fignum,advancedvec(i));set(gca,'fontsize',advancedvec(i).bigfontsize);
%title(sprintf('%s\nSequential Bayes vs. One-Shot Bayes, %s',advancedvec(i).titlestring,subtitle),'FontSize',advancedvec(i).bigfontsize,'FontName',advancedvec(i).fontname);
xlabel('Prior mean','FontSize',advancedvec(i).bigfontsize,'FontName',advancedvec(i).fontname);
ylabel('Reward: Optimal Bayes vs. Bayes One-Shot','FontSize',advancedvec(i).bigfontsize,'FontName',advancedvec(i).fontname); 
if advancedvec(i).saveplot
    texttodifferentiate = 'OptVsOneshot';
%    UtilSaveFigFile(fignum, advancedvec(i).dirstring, advancedvec(i).filestring, texttodifferentiate, advancedvec(i).graphicextension);
    UtilSaveFigEpsPdf(fignum,advancedvec(i).dirstring,strcat(advancedvec(i).dirstring,texttodifferentiate),'-r600');
end


%%%%%%%%%%%%%%%% Plot probability of reversal of decision for sequential
fignum = fignum+1; figure(fignum); hold off;
hold off;
for i=1:veclen
    plot(advancedvec(i).simFreqDeltaVec,matvec(i).outBayes.PrReversalAve(:),strcat(mycode{1+mod(i-1,length(mycode))},mycolor{1+mod(i-1,length(mycolor))} ),'LineWidth',linewid);
    hold on;
end
if length(legendvec) > 1 
    legend(legendvec,'Location','NorthWest'); 
end;
i=1;
%    UtilResizeFigureToBounds(mat);
tmp=axis();
tmp(1)=MULOW;
tmp(2)=MUHIGH;
axis(tmp);

UtilStdizeFigure(fignum,advancedvec(i));set(gca,'fontsize',advancedvec(i).bigfontsize);
%title(sprintf('%s\nSequential Bayes, %s',advancedvec(i).titlestring,subtitle),'FontSize',advancedvec(i).bigfontsize,'FontName',advancedvec(i).fontname);
xlabel('Prior mean','FontSize',advancedvec(i).bigfontsize,'FontName',advancedvec(i).fontname);
ylabel('Prob(Reverse decision)','FontSize',advancedvec(i).bigfontsize,'FontName',advancedvec(i).fontname); 
if advancedvec(i).saveplot
    texttodifferentiate = 'Reverse';
%    UtilSaveFigFile(fignum, advancedvec(i).dirstring, advancedvec(i).filestring, texttodifferentiate, advancedvec(i).graphicextension);
    UtilSaveFigEpsPdf(fignum,advancedvec(i).dirstring,strcat(advancedvec(i).dirstring,texttodifferentiate),'-r600');
end


%%%%%%%%%%%%%%%% Plot expected number of samples seen
fignum = fignum+1; figure(fignum); hold off;
hold off;
for i=1:veclen
    plot(advancedvec(i).simFreqDeltaVec,matvec(i).outBayes.ENumSeenAve(:),strcat(mycode{1+mod(i-1,length(mycode))},mycolor{1+mod(i-1,length(mycolor))} ),'LineWidth',linewid);
    hold on;
end
if length(legendvec) > 1 
    legend(legendvec,'Location','NorthWest'); 
end;
for i=1:veclen    % Plot the A, B, C, D
    plot(advancedvec(i).simFreqDeltaVec,matvec(i).outBayes.ENumSeenAve(:),strcat(mycode{1+mod(i-1,length(mycode))},mycolor{1+mod(i-1,length(mycolor))} ),'LineWidth',linewid);
    ist0 = 1 ; % set to true if ABCD to be plotted at t0, false if to be plotted on horizontal axis
    ishoriz = true; 
    if i==1
        UtilPlotABCD( basicvec(i), advancedvec(i), matvec(i), ishoriz, ist0, [0 0] );  
    else
        UtilPlotABCD( basicvec(i), advancedvec(i), matvec(i), ishoriz, ist0, [basicvec(i).tau basicvec(i).tau] );  
    end
end
i=1;
%    UtilResizeFigureToBounds(mat);
tmp=axis();
tmp(1)=MULOW;
tmp(2)=MUHIGH;
tmp(4)= max(tmp(4),basicvec(i).TMax);
axis(tmp);
UtilStdizeFigure(fignum,advancedvec(i));set(gca,'fontsize',advancedvec(i).bigfontsize);
%title(sprintf('%s\nSequential Bayes, E[Num seen], %s',advancedvec(i).titlestring,subtitle),'FontSize',advancedvec(i).bigfontsize,'FontName',advancedvec(i).fontname);
xlabel('prior mean','FontSize',advancedvec(i).bigfontsize,'FontName',advancedvec(i).fontname);
ylabel('Expected sample size','FontSize',advancedvec(i).bigfontsize,'FontName',advancedvec(i).fontname); 
if advancedvec(i).saveplot
    texttodifferentiate = 'NumStarted';
    UtilSaveFigEpsPdf(fignum,advancedvec(i).dirstring,strcat(advancedvec(i).dirstring,texttodifferentiate),'-r600');
%    UtilSaveFigFile(fignum, advancedvec(i).dirstring, advancedvec(i).filestring, texttodifferentiate, advancedvec(i).graphicextension);
end

%%%%%%%%%%%%%%%% Plot stopping boundaries
fignum = fignum+1; figure(fignum); hold off;
hold off;
for i=1:veclen
    plot(matvec(i).tvec,matvec(i).bndupper,strcat(mycode{1+mod(i-1,length(mycode))},mycolor{1+mod(i-1,length(mycolor))} ),'LineWidth',linewid);
    hold on;
    
end
if length(legendvec) > 1 
    legend(legendvec,'Location','NorthEast'); 
end;
for i=1:veclen
    plot(matvec(i).tvec,matvec(i).bndlower,strcat(mycode{1+mod(i-1,length(mycode))},mycolor{1+mod(i-1,length(mycolor))} ),'LineWidth',linewid);
% 	if false
%         UtilPlotABCD( basicvec(i), advancedvec(i), matvec(i), false, true, [0 0] );  
%     else
%         UtilPlotABCD( basicvec(i), advancedvec(i), matvec(i), false, true, [1*basicvec(i).tau 1.5*basicvec(i).tau] );  
%     end
hold on
end
for i=1:veclen
    if abs( matvec(i).Threshpoint(2) - matvec(i).Threshpoint(4) ) > 1
        scatter (matvec(i).bestsvec(matvec(i).Threshpoint(2):matvec(i).Threshpoint(4)) + basic.t0, matvec(i).muvec(matvec(i).Threshpoint(2):matvec(i).Threshpoint(4)),strcat(mycode2{1+mod(i-1,length(mycode))},mycolor{1+mod(i-1,length(mycolor))} ));
        hold on ;
    end
    % plot upper Stage 1 boundary if there exists a Stage 1
    if abs( matvec(i).Threshpoint(3) - matvec(i).Threshpoint(1)) > 1
        scatter (matvec(i).bestsvec(matvec(i).Threshpoint(3):matvec(i).Threshpoint(1)) + basic.t0, matvec(i).muvec(matvec(i).Threshpoint(3):matvec(i).Threshpoint(1)),strcat(mycode2{1+mod(i-1,length(mycode))},mycolor{1+mod(i-1,length(mycolor))} ));
        hold on ;
    end
end

i=1;
%    UtilResizeFigureToBounds(mat);
tmp=axis();
tmp(1)=0;
tmp(3)=MULOW;
tmp(4)=MUHIGH;
plotmaxreps=tmp(2);
axis(tmp);
UtilStdizeFigure(fignum,advancedvec(i));set(gca,'fontsize',advancedvec(i).bigfontsize);
%title(sprintf('%s\nSequential Bayes, stage II boundaries, %s',advancedvec(i).titlestring,subtitle),'FontSize',advancedvec(i).bigfontsize,'FontName',advancedvec(i).fontname);
ylabel('posterior mean','FontSize',advancedvec(i).bigfontsize,'FontName',advancedvec(i).fontname);
xlabel('n_0 + t','FontSize',advancedvec(i).bigfontsize,'FontName',advancedvec(i).fontname); 
if advancedvec(i).saveplot
    texttodifferentiate = 'Bounds';
%    UtilSaveFigFile(fignum, advancedvec(i).dirstring, advancedvec(i).filestring, texttodifferentiate, advancedvec(i).graphicextension);
    UtilSaveFigEpsPdf(fignum,advancedvec(i).dirstring,strcat(advancedvec(i).dirstring,texttodifferentiate),'-r600');
end

fignum = fignum+1; figure(fignum); hold off;
hold off;
for i=1:veclen
    plot(matvec(i).tvec-basicvec(i).tau,matvec(i).bndupper,strcat(mycode{1+mod(i-1,length(mycode))},mycolor{1+mod(i-1,length(mycolor))} ),'LineWidth',linewid);
    hold on;
end
if length(legendvec) > 1 
    legend(legendvec,'Location','NorthEast'); 
end;
for i=1:veclen
    plot(matvec(i).tvec-basicvec(i).tau,matvec(i).bndlower,strcat(mycode{1+mod(i-1,length(mycode))},mycolor{1+mod(i-1,length(mycolor))} ),'LineWidth',linewid);
end
i=1;
%    UtilResizeFigureToBounds(mat);
tmp=axis();
tmp(1)=0;
tmp(3)=MULOW;
tmp(4)=MUHIGH;
plotmaxreps=tmp(2);
axis(tmp);
UtilStdizeFigure(fignum,advancedvec(i));set(gca,'fontsize',advancedvec(i).bigfontsize);
title(sprintf('%s\nSequential Bayes, stage II boundaries, %s',advancedvec(i).titlestring,subtitle),'FontSize',advancedvec(i).bigfontsize,'FontName',advancedvec(i).fontname);
ylabel('posterior mean','FontSize',advancedvec(i).bigfontsize,'FontName',advancedvec(i).fontname);
xlabel('n_0 + number of patients SEEN','FontSize',advancedvec(i).bigfontsize,'FontName',advancedvec(i).fontname); 
if advancedvec(i).saveplot
    texttodifferentiate = 'BoundsLessTau';
%    UtilSaveFigFile(fignum, advancedvec(i).dirstring, advancedvec(i).filestring, texttodifferentiate, advancedvec(i).graphicextension);
    UtilSaveFigEpsPdf(fignum,advancedvec(i).dirstring,strcat(advancedvec(i).dirstring,texttodifferentiate),'-r600');
end
