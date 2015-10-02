function [ figout ] = DelayDoContours( fignum, basic, advanced, mat )
%DELAYDOCONTOURS generates several contour plots, starting with figure
%'fignum', and using the basic and advanced input parameters set up by the
%DelayInput*.m functions, which describe a sequential stopping problem with
%delayedr response. The mat structure should have the data as computed by
%the DelayCurvesRecur followed by the DelayStageOne routines.
%
% figout is set to the next figure which would appear in a sequence of figures to follow
%
% Designed for Delay Sequential Trials paper of Chick, Forster, Pertile (alpha order).
%
% (c) 2014, S Chick
% Created: 14 April 2014
% Last touched: 14 April 2014
%
    tvec = mat.tvec;
    bndupper = mat.bndupper;
    bndlower = mat.bndlower;

    if isfield(advanced,'subtitle')
        subtitle = sprintf(' (%s)', advanced.subtitle);
    else
        subtitle = '';
    end
    
    probcontours = [.01 .1:.2:0.9]; % level sets for probabilities which are displayed
    numcontours = 5;   % number of level sets to display for other values
    boundlinecode = '-.';
    
    if advanced.saveplot % If files to be saved, write information about the version of code used.
        UtilSaveVersionFile( advanced.dirstring );
    end
    
    %%%%%% Contour: Probability best alternative eventually selected 
    fignum = fignum+1; figure(fignum); 
    probcontours = [.5:.1:.9 .95 .99]; % level sets for probabilities which are displayed
    [~, h]=contour(mat.tmat,mat.muvec,mat.PCSmat,probcontours);
    UtilResizeFigureToBounds(mat);
    UtilStdizeFigure(fignum,advanced);
    set(h,'ShowText','on')
    hold on;
    plot(tvec,bndupper,boundlinecode,tvec,bndlower,boundlinecode);
    title(sprintf('%s\nProb(eventually select best) %s',advanced.titlestring,subtitle),'FontSize',advanced.bigfontsize,'FontName',advanced.fontname);
    xlabel('t_0 + num patients started','FontSize',advanced.smallfontsize,'FontName',advanced.fontname); 
    ylabel('mean reward','FontSize',advanced.smallfontsize,'FontName',advanced.fontname);
    %axis('square');
    if advanced.saveplot
        texttodifferentiate = 'PrPCS';
        UtilSaveFigFile(fignum, advanced.dirstring, advanced.filestring, texttodifferentiate, advanced.graphicextension);
    end

    probcontours = [.01 .1:.2:0.9 0.99]; % level sets for probabilities which are displayed
    %%%%%% Contour: Probability upper boundary eventually selected 
    fignum = fignum+1; figure(fignum); 
    [~, h]=contour(mat.tmat,mat.muvec,mat.Puppermat,probcontours);
    UtilResizeFigureToBounds(mat);
    UtilStdizeFigure(fignum,advanced);
    set(h,'ShowText','on')
    hold on;
    plot(tvec,bndupper,boundlinecode,tvec,bndlower,boundlinecode);
    title(sprintf('%s\nProb(exit upper bound) %s',advanced.titlestring,subtitle),'FontSize',advanced.bigfontsize,'FontName',advanced.fontname);
    xlabel('t_0 + num patients started','FontSize',advanced.smallfontsize,'FontName',advanced.fontname); 
    ylabel('mean reward','FontSize',advanced.smallfontsize,'FontName',advanced.fontname);
    %axis('square');
    if advanced.saveplot
        texttodifferentiate = 'PrExitUpper';
        UtilSaveFigFile(fignum, advanced.dirstring, advanced.filestring, texttodifferentiate, advanced.graphicextension);
    end

    %%%%%% Contour: Probability new treatment eventually selected 
    
    fignum = fignum+1; figure(fignum); 
    [~, h]=contour(mat.tmat,mat.muvec,mat.Pnewmat,probcontours);
    UtilResizeFigureToBounds(mat);
    UtilStdizeFigure(fignum,advanced);
    set(h,'ShowText','on')
    hold on;
    plot(tvec,bndupper,boundlinecode,tvec,bndlower,boundlinecode);
    title(sprintf('%s\nProb(pick new treatment) %s',advanced.titlestring,subtitle),'FontSize',advanced.bigfontsize,'FontName',advanced.fontname);
    xlabel('t_0 + num patients started','FontSize',advanced.smallfontsize,'FontName',advanced.fontname); 
    ylabel('mean reward','FontSize',advanced.smallfontsize,'FontName',advanced.fontname);
    if advanced.saveplot
        texttodifferentiate = 'PrSelectNew';
        UtilSaveFigFile(fignum, advanced.dirstring, advanced.filestring, texttodifferentiate, advanced.graphicextension);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
% contour plot of value function: for max, look at B0 value at upper
% boundary of stopping boundary
    maxval =  10^ceil(.1+log10(max(mat.B0mat(:,1))));  
    minval = min(min(mat.B0mat));    % needs a bit of defensive programming here
    fignum = fignum+1; figure(fignum); 
%    x = minval:((maxval-minval)/(numcontours+1)):maxval;
%    [~, h]=contour(mat.tmat,mat.muvec,mat.B0mat,x);
    [~, h]=contour(mat.tmat,mat.muvec,mat.B0mat);
    UtilResizeFigureToBounds(mat);
    UtilStdizeFigure(fignum,advanced);
    set(h,'ShowText','on')
    hold on;
    plot(tvec,bndupper,boundlinecode,tvec,bndlower,boundlinecode);
    title(sprintf('%s\nValue function %s',advanced.titlestring,subtitle),'FontSize',advanced.bigfontsize,'FontName',advanced.fontname);
    xlabel('t_0 + num patients started','FontSize',advanced.smallfontsize,'FontName',advanced.fontname); 
    ylabel('mean reward','FontSize',advanced.smallfontsize,'FontName',advanced.fontname);
    if advanced.saveplot
        texttodifferentiate = 'ValFunction';
        UtilSaveFigFile(fignum, advanced.dirstring, advanced.filestring, texttodifferentiate, advanced.graphicextension);
    end
%    fignum = fignum+1; figure(fignum); 
%    [C h]=contour3(mat.tmat,mat.muvec,mat.B0mat,valuecontours);

% contour plot of regret
%    Regretcontmat = mat.RewardPIMatII - mat.RewardMatII;
    Regretcontmat = max(0,mat.RewardPIMatII - mat.B0mat); % the max(0,...) catches the -e-7 type roundoff errors   
    maxval =  10^ceil(.1+log10(max(max(Regretcontmat))));  
    minval = 0;    % needs a bit of defensive programming here
    fignum = fignum+1; figure(fignum); 
%	x = minval:((maxval-minval)/(numcontours+1)):maxval;
%   x = 10.^(maxval - (0:numcontours-1)); x = fix(x./10.^floor(log10(x))).*10.^floor(log10(x));
%    [~, h]=contour(mat.tmat,mat.muvec,Regretcontmat,x);
    [~, h]=contour(mat.tmat,mat.muvec,Regretcontmat);
    UtilResizeFigureToBounds(mat);
    UtilStdizeFigure(fignum,advanced);
    set(h,'ShowText','on')
    hold on;
    plot(tvec,bndupper,boundlinecode,tvec,bndlower,boundlinecode);
    title(sprintf('%s\nE[Regret | Stage II] %s',advanced.titlestring,subtitle),'FontSize',advanced.bigfontsize,'FontName',advanced.fontname);
    xlabel('t_0 + num patients started','FontSize',advanced.smallfontsize,'FontName',advanced.fontname); 
    ylabel('Regret','FontSize',advanced.smallfontsize,'FontName',advanced.fontname);
    if advanced.saveplot
        texttodifferentiate = 'RegreStageII';
        UtilSaveFigFile(fignum, advanced.dirstring, advanced.filestring, texttodifferentiate, advanced.graphicextension);
    end
 %   fignum = fignum+1; figure(fignum); 
 %   [C h]=contour3(mat.tmat,mat.muvec,Regretcontmat,valuecontours);

% contour plot of reward
    maxval = ceil(log10(max(max(mat.RewardMatII))));      % FIX: WARN: code currently assumes bndupper goes positive and bndlower goes negative.
    minval = min(min(mat.RewardMatII)); %0;    % needs a bit of defensive programming here
    fignum = fignum+1; figure(fignum); 
%	x = 10.^(maxval - (0:numcontours-1)/2); x = fix(x./10.^floor(log10(x))).*10.^floor(log10(x));
%	x = minval:((maxval-minval)/(numcontours+1)):maxval;
%    [~, h]=contour(mat.tmat,mat.muvec,mat.RewardMatII,x);
    [~, h]=contour(mat.tmat,mat.muvec,mat.RewardMatII);
    UtilResizeFigureToBounds(mat);
    UtilStdizeFigure(fignum,advanced);
    set(h,'ShowText','on')
    hold on;
    plot(tvec,bndupper,boundlinecode,tvec,bndlower,boundlinecode);
    title(sprintf('%s\nE[Reward | Stage II] %s',advanced.titlestring,subtitle),'FontSize',advanced.bigfontsize,'FontName',advanced.fontname);
    xlabel('t_0 + num patients started','FontSize',advanced.smallfontsize,'FontName',advanced.fontname); 
    ylabel('mean reward','FontSize',advanced.smallfontsize,'FontName',advanced.fontname);
    if advanced.saveplot
        texttodifferentiate = 'RewardStageII';
        UtilSaveFigFile(fignum, advanced.dirstring, advanced.filestring, texttodifferentiate, advanced.graphicextension);
    end
 %   fignum = fignum+1; figure(fignum); 
 %   [C h]=contour3(mat.tmat,mat.muvec,mat.RewardMatII,valuecontours);

    figout = fignum + 1;    % do this to allow calling routine to continue with next plot at higher numbered plots without overwriting plots
end

