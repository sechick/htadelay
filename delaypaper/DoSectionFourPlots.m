function [figout, mat] = DoSectionFourPlots(fignum, basic,advanced,PathReps,ProductionReps,graphicsuffix,colorflag)
% DoSectionFourPlots: This function generates plots for section which
% demonstrates sample paths and basic sensitivity analysis for the delay
% sequential sampling paper.
%
% CALLED BY: DelayExperimentsForPaper 
%
% 21 Mar 2015: Forster, Chick Pertile

% do some error checking of input data
[basic, advanced, rval, msgs] = DelayInputValidator( basic, advanced );
if ~rval, msgs, end; 

if nargin < 7
    colorflag = false;
end

% set up some local data
%dirname = advanced.dirstring; % directory for putting figures
dirname = 'Figure'; % directory for putting figures
linewid = 2;

% Start running the analysis
advanced.PLOTSIMS = true ;
[~, mat] = DelayCurvesRecur(basic, advanced);
[mat] = DelayStageOne(basic, advanced, mat ); 
toneshot = (0:basic.TMax);
[mat]  = DelayOptimalBayesOneStage (basic, advanced, toneshot, mat ); % this function that Steve wrote defines the optimal fixed sample size based on a comparison of EVSI and sampling costs

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% For Figure 1: Generate a few sample paths and show the general
%%%%%% features of the stopping times etc.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

advanced.simNumReps = PathReps;  % set to specific number of replications required for this specific experiment.
advanced.CRN = true;
advanced.CRNAcrossExperiment = 37; 
[ fignum, mat ] = DelaySimOverview( fignum, basic, advanced, mat ) ;

% work the statistics with the value of mu0 which was tested which is
% closest to I/P, the break even point.
[smallestMuInGrid, indexToMuOfInterest]=min(abs(advanced.simFreqDeltaVec - basic.ICost/basic.PPatients));  

% set up a few data structures to be able to plot stopping times
% boundaries, etc. with a few sample paths. Pick the indices for outtime3,
% outtime6, outtime7 and outtime10 to point to the sample paths of interest
tmpvec = mat.simBayesOut(indexToMuOfInterest).ENumSeen; % get number of samples seen by time of decision to stop
tmpvec(tmpvec>=basic.tau) = tmpvec(tmpvec>=basic.tau) - floor(basic.tau);
outtimes = horzcat( tmpvec, ...
    mat.simBayesOut(indexToMuOfInterest).PrHitTop, ...
    mat.simBayesOut(indexToMuOfInterest).PrHitBot, ...
    1 - mat.simBayesOut(indexToMuOfInterest).PrHitTop - mat.simBayesOut(indexToMuOfInterest).PrHitBot, ...
    mat.simBayesOut(indexToMuOfInterest).RealizedMean );
paths = horzcat ( mat.simBayesOut(indexToMuOfInterest).NoStatsBoundTvec,  mat.simBayesOut(indexToMuOfInterest).NoStatsSomePostMean);
outtime3 = outtimes( 3, 1 ) ; 
decidetime3 = outtime3 + floor(basic.tau) ; % Use floors here to keep indices as integers.
outtime6 = outtimes( 6, 1 ) ; 
decidetime6 = outtime6 + floor(basic.tau) ; 
outtime7 = outtimes( 7, 1 ) ; 
decidetime7 = outtime7 + floor(basic.tau) ; 
outtime10 = outtimes( 10, 1 ) ; 
decidetime10 = outtime10 + floor(basic.tau) ; 

% Plotting the optimal stopping rule
fignum = fignum + 1; figure(fignum); % always increment fignum and create a new figure
hold off ; 
tmpfontsize = 10;
if colorflag
    code1='b'; code2=':r';
else 
    code1='k'; code2=':k';
end
plot( paths( 1:outtime6, 1 )', paths( 1 : outtime6, 7 )', 'k' ,'LineWidth',linewid ) ;
hold on ; 
plot( paths( outtime6 + 1:decidetime6, 1 )', paths( outtime6 + 1:decidetime6, 7 )', code2 ,'LineWidth',linewid ) ;
%rect=[1200,-5000,100,100]; %'Position',[220,55,10,10],, 'Fontsize', tmpfontsize
lgd=legend( 'Actual Stage II (observing and sampling)', 'Actual Stage III (observing only)', 'AutoUpdate', 'Off' ) ; 
plot( paths( 1:outtime3, 1 )', paths( 1 : outtime3, 4 )', 'k', paths( outtime3 + 1:decidetime3, 1 )', paths( outtime3 + 1:decidetime3, 4 )', code2 ,'LineWidth',linewid); 
plot( paths( 1:outtime7, 1 )', paths( 1 : outtime7, 8 )', 'k', paths( outtime7 + 1:decidetime7, 1 )', paths( outtime7 + 1:decidetime7, 8 )', code2  ,'LineWidth',linewid); 
plot( paths( 1:outtime10, 1 )', paths( 1 : outtime10, 11 )', 'k', paths( outtime10 + 1:decidetime10, 1 )', paths( outtime10 + 1:decidetime10, 11 )', code2 ,'LineWidth',linewid); 
set(lgd,'FontSize',tmpfontsize);
set(gca,'FontSize',40)
% plot lower Stage 1 boundary
scatter (mat.bestsvec(mat.Threshpoint(2):mat.Threshpoint(4)) + basic.t0, mat.muvec(mat.Threshpoint(2):mat.Threshpoint(4)),code1);
% plot upper Stage 1 boundary
scatter (mat.bestsvec(mat.Threshpoint(3):mat.Threshpoint(1)) + basic.t0, mat.muvec(mat.Threshpoint(3):mat.Threshpoint(1)),code1);
ylabel('Prior / Posterior Mean', 'Fontsize', tmpfontsize);
xlabel('n_0 + t', 'Fontsize', tmpfontsize);
ishoriz = false; % set to true if muvec is on horizontal axis, false if on vertical axis
plot_upper = 6000 ; % set upper limit of plot region ( + / 
plot_lower = -6000  ; % set lower limit of plot region
%advanced.smallfontsize = 11;
ylim( [ plot_lower,plot_upper] ) ;
line( [ basic.t0 + basic.tau, basic.t0 + basic.tau ], [ plot_lower, plot_upper ], 'Linestyle', ':', 'Color', 'k' ) ;
line( [ basic.t0 + basic.TMax, basic.t0 + basic.TMax ], [ plot_lower, plot_upper ], 'Linestyle', ':', 'Color', 'k' ) ;
line( [ basic.t0, basic.t0 ], [ plot_lower, plot_upper ], 'Linestyle', ':', 'Color', 'k' ) ;
line( [ basic.t0 + basic.tau + basic.TMax, basic.t0 + basic.tau + basic.TMax ], [ plot_lower, plot_upper ], 'Linestyle', ':', 'Color', 'k' ) ;
text(  basic.t0,  plot_lower + 300, 'n_0', 'Fontsize', tmpfontsize) ;
text(  basic.t0 + basic.TMax,  plot_lower + 300 , 'n_0 + T_{max}', 'Fontsize', tmpfontsize) ;
text(  basic.t0 + basic.tau,  plot_lower + 300 , 'n_0 + \tau', 'Fontsize', tmpfontsize) ;
text(  basic.t0 + basic.tau + basic.TMax - 420,  plot_lower + 300 , 'n_0 + \tau + T_{max}', 'Fontsize', tmpfontsize) ;
if (length(graphicsuffix) == 0) || colorflag
    tmpfontsize = round(advanced.smallfontsize*.9);
    text(  ( basic.t0 + basic.tau/20 ) , mat.muvec(mat.Threshpoint(1)) + (plot_upper-mat.muvec(mat.Threshpoint(1)))/2 , 'No trial', 'Fontsize', tmpfontsize);
    text(  ( basic.t0 + basic.tau/20 ) , plot_lower + (mat.muvec(mat.Threshpoint(2))-plot_lower)/2 , 'No trial', 'Fontsize', tmpfontsize) ;
    if abs(mat.Threshpoint(1) - mat.Threshpoint(3)) > 1
        text(  ( basic.t0 + basic.tau/20 ) , (mat.muvec(mat.Threshpoint(1)) + mat.muvec(mat.Threshpoint(3)))/2 , 'Fixed trial', 'Fontsize', tmpfontsize) ;
    end
    if abs(mat.Threshpoint(4) - mat.Threshpoint(2)) > 1
        text(  ( basic.t0 + basic.tau/20 ) ,  (mat.muvec(mat.Threshpoint(2)) + mat.muvec(mat.Threshpoint(4)))/2 , 'Fixed trial', 'Fontsize', tmpfontsize) ;
    end
    tmpfontsize = round(advanced.smallfontsize);
    text(  ( basic.t0 + basic.tau/20 ) , (mat.muvec(mat.Threshpoint(3)) + mat.muvec(mat.Threshpoint(4)))/2 , 'Sequential trial recruitment', 'Fontsize', tmpfontsize) ;
    text(  ( basic.t0 + basic.tau/20 ) , plot_upper - (plot_upper-plot_lower)/40, 'Maximum extent of stage I', 'Fontsize', tmpfontsize) ;
    text(  ( basic.t0 + basic.tau*(1+1/20) ),  plot_upper - (plot_upper-plot_lower)/40, 'Maximum extent of stage II', 'Fontsize', tmpfontsize) ;
    text(  ( basic.t0 + basic.tau*(1/20) + basic.TMax ),  plot_upper - (plot_upper-plot_lower)/40, 'Maximum extent of stage III', 'Fontsize', tmpfontsize) ;
end
ist0 = true ; % set to true if ABCD to be plotted at t0, false if to be plotted on horizontal axis
UtilPlotABCD_Sec4( basic, advanced, mat, ishoriz, ist0 );  
plot (mat.tvec, mat.bndupper,'--k',mat.tvec ,mat.bndlower, '--k'  );
text(  paths( decidetime3 ), outtimes( 3, 5) , '*1', 'Fontsize',tmpfontsize) ;
text(  paths( decidetime6 ), outtimes( 6, 5 ) , '*2', 'Fontsize', tmpfontsize) ;
text(  paths( decidetime7 -2 ), outtimes( 7, 5 ) , '*4', 'Fontsize',tmpfontsize) ;
text(  paths( decidetime7 -4 ), paths( decidetime7 , 8 ) - 10 , '4', 'Fontsize', tmpfontsize) ;
text(  paths( decidetime10 ), outtimes( 10, 5) , '*3', 'Fontsize', tmpfontsize) ;
% xlim( [ 0 basic.TMax + basic.t0 + basic.tau ] ); % Baseline
xlim( [ 0 3100 ] ); % Comparator

UtilStdizeFigure_Sec4(fignum,advanced,false);
if colorflag
    title('Optimal Bayes Sequential stopping boundaries','FontSize',18,'FontName','Helvetica');
    ylabel('Posterior mean','FontSize',16,'FontName','Helvetica');
    xlabel('n_{0} + t','FontSize',16,'FontName','Helvetica'); 
    hold on ;
end

UtilSaveFigEpsPdf(fignum,dirname,strcat('paths',graphicsuffix),'-r600');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% For Other figures, use estimates based on a lot of sample
%%%%%% replications. so, the simulations need to be rerun.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~colorflag % only do other figures if color plot is not needed
    advanced.simNumReps = ProductionReps;  % set to specific number of replications required for this specific experiment.
    [ fignum, mat ] = DelaySimOverview( fignum, basic, advanced, mat ) ;
    fignum = DelaySimOutput( fignum, basic, advanced, mat );  

    % calculating the gain in terms of "realized reward" from using our method.
    % MF 08/03: diff_trial not used in Section 4 because we have no reference
    % experiment. SEs are also included to see whether they make a difference
    % when running 15,000 replications. These SE estimates are larger than they
    % really are if CRN is used.
    diff_all = mat.outBayes.RealizedRewardAve - mat.outBayes.RealizedRewardAllAve; % difference between the "realized reward" of our methods and that of taking always a sample size equal to T_max
    %diff_trial = mat.outBayes.RealizedRewardAve - mat.outBayes.RealizedRewardTrialAve; % difference between the "realized reward" of our methods and that of taking always a sample size equal to the sample size in the trial we took the data from
    diff_fixed = mat.outBayes.RealizedRewardAve - mat.outBayes.RealizedRewardFixedAve; % difference between the "realized reward" of our methods and that of taking always an optimal fixed sample size (i.e. ignoring Stage 2) 
    diff_all_approxSEifindep = sqrt(mat.outBayes.RealizedRewardSE.^2 + mat.outBayes.RealizedRewardAllSE.^2);
    % diff_trial_approxSEifindep = sqrt(mat.outBayes.RealizedRewardSE.^2 + mat.outBayes.RealizedRewardTrialSE.^2);
    diff_fixed_approxSEifindep = sqrt(mat.outBayes.RealizedRewardSE.^2 + mat.outBayes.RealizedRewardFixedSE.^2);

    % MF 08/03: plot the differences in the expected rewards (Figure 3 of the manuscript)
    fignum = fignum + 1; figure(fignum); % always increment fignum and create a new figure
    hold off ; 
    zval = 1.96;
    plot_upper = 6000 ; % set upper limit of plot region 
    plot_lower = -6000  ; % set lower limit of plot region
    % MF 04/03: Plotting the gain in terms of realised reward using our method 
    ishoriz = true ; % set to true because muvec is on horizontal axis
    ist0 = false ; % set to true if ABCD to be plotted at t0, false if to be plotted on horizontal axis
    % Plot the differences for fixed and all
    plot(  advanced.simFreqDeltaVec, diff_fixed, 'k', advanced.simFreqDeltaVec, diff_all, '-.k' ,'LineWidth',linewid ) ;  
    xlim( [ plot_lower, plot_upper] ) ;
    ylim( [ -200000, 1000000 ] ) ; % Comparator
    xlabel( 'Prior mean', 'Fontsize', advanced.bigfontsize ) ; 
    ylabel( 'Difference between expected rewards', 'Fontsize', advanced.bigfontsize) ; 
    legend( 'Optimal Bayes One Stage', 'Fixed', 'Location', 'north','AutoUpdate', 'Off' ) ; 
    hold on ; 
    % MF 08/03: Steve's code for SEs: these correctly compute SE estimate when
    % output is run with CRN across experiments.
    if advanced.keepAllOutput   % if all the simulation output was kept, try to produce some CIs for differences based on the use of CRN (if it was used)
        difftrialvecAVE =zeros(length(advanced.simFreqDeltaVec));
        %difftrialvecSTD =difftrialvecAVE;
        diffallvecAVE = difftrialvecAVE;
        diffallvecSTD = difftrialvecAVE;
        difffixedvecAVE =  difftrialvecAVE;
        difffixedvecSTD = difftrialvecAVE;
        for i=1:length(advanced.simFreqDeltaVec)
            diffall = mat.simBayesOut(i).RealizedReward - mat.simBayesOut(i).RealizedRewardAll;
            diffallvecAVE(i) = mean(diffall);
            diffallvecSTD(i) = std(diffall)/sqrt(length(diffall));
            difffixed = mat.simBayesOut(i).RealizedReward - mat.simBayesOut(i).RealizedRewardFixed;
            difffixedvecAVE(i) = mean(difffixed);
            difffixedvecSTD(i) = std(difffixed)/sqrt(length(diffall));
            %difftrial = mat.simBayesOut(i).RealizedReward - mat.simBayesOut(i).RealizedRewardTrial;
            %difftrialvecAVE(i) = mean(difftrial);
            %difftrialvecSTD(i) = std(difftrial)/sqrt(length(diffall));
        end
        %fignum = fignum + 1; figure(fignum);
        plot( advanced.simFreqDeltaVec, diffallvecAVE + zval * diffallvecSTD, ':k',advanced.simFreqDeltaVec, diffallvecAVE - zval * diffallvecSTD, ':k' ,'LineWidth',linewid ) ;
        hold on ; 
        %title(sprintf('Difference in Outcomes between Sequential and All Pairs Trials\n%s',advanced.titlestring),'FontSize',advanced.bigfontsize,'FontName',advanced.fontname)
        %xlabel('posterior mean per patient','FontSize',advanced.smallfontsize,'FontName',advanced.fontname);
        %ylabel('population benefit from sequential trial','FontSize',advanced.smallfontsize,'FontName',advanced.fontname);
        %legend('difference',sprintf('difference +/- %3.2f SE',zval));
        %UtilStdizeFigure( fignum, advanced ); 
        %fignum = fignum + 1; figure(fignum);
        plot ( advanced.simFreqDeltaVec, difffixedvecAVE + zval * difffixedvecSTD, ':k', advanced.simFreqDeltaVec, difffixedvecAVE - zval * difffixedvecSTD, ':k') ;
        hold on ; 
        %title(sprintf('Difference in Outcomes between Sequential and Fixed Samples Trials\n%s',advanced.titlestring),'FontSize',advanced.bigfontsize,'FontName',advanced.fontname)
        %xlabel('posterior mean per patient','FontSize',advanced.smallfontsize,'FontName',advanced.fontname);
        %ylabel('population benefit from sequential trial','FontSize',advanced.smallfontsize,'FontName',advanced.fontname);
        %legend('difference',sprintf('difference +/- %3.2f SE',zval));
        %UtilStdizeFigure( fignum, advanced ); 
        %fignum = fignum + 1; figure(fignum);
        %plot (advanced.simFreqDeltaVec, difftrialvecAVE, '-', ...
        %     advanced.simFreqDeltaVec, difftrialvecAVE + zval * difftrialvecSTD, '-.', ...
        %     advanced.simFreqDeltaVec, difftrialvecAVE - zval * difftrialvecSTD, '-.');
        %title('Difference in outcomes between sequential and trial pairs')
        %title(sprintf('Difference in Outcomes between Sequential and Trial\n%s',advanced.titlestring),'FontSize',advanced.bigfontsize,'FontName',advanced.fontname)
        %xlabel('posterior mean per patient','FontSize',advanced.smallfontsize,'FontName',advanced.fontname);
        %ylabel('population benefit from sequential trial','FontSize',advanced.smallfontsize,'FontName',advanced.fontname);
        %legend('difference',sprintf('difference +/- %3.2f SE',zval));
    end
    lohivec = [ -2000, 1000000 ] ;
    UtilPlotABCD( basic, advanced, mat, ishoriz, ist0, lohivec );  
    UtilStdizeFigure(fignum,advanced);
    UtilSaveFigEpsPdf(fignum,dirname, strcat('diffs_baseline',graphicsuffix),'-r600');

    % SC's commands for plotting the differences in realized rewards (would expect these
    % differences to be always positive over the range of prior mean
    % considered)
    %fignum = fignum + 1; figure(fignum);
    %plot (advanced.simFreqDeltaVec, diff_trial, '-', advanced.simFreqDeltaVec, diff_fixed,'-')
    %xlim([mat.muvec(mat.Threshpoint(2)) mat.muvec(mat.Threshpoint(1))]);
    %legend('trial','fixed')
    %UtilStdizeFigure( fignum, advanced )

    %Plotting the probability of making the correct decision (Figure 5 of the
    %manuscript)

    fignum = fignum + 1; figure(fignum);
    hold off ; 
    ishoriz = true; % set to true if muvec is on horizontal axis, false if on vertical axis
    ist0 = false ;
    %fignum = fignum + 1; figure(fignum);
    plot(advanced.simFreqDeltaVec, mat.outBayes.PrTrueBestSelectedAve, '--k', advanced.simFreqDeltaVec, mat.outBayes.PrTrueBestSelectedAllAve, '-.k', advanced.simFreqDeltaVec, mat.outBayes.PrTrueBestSelectedFixedAve, '-k' ,'LineWidth',linewid ); 
   hold on ; 
     xlabel( 'Prior mean', 'Fontsize', advanced.bigfontsize ) ; 
    ylabel( 'Proportion of correct decisions', 'Fontsize', advanced.bigfontsize ) ; 
    legend( 'Optimal Bayes Sequential',  'Fixed', 'Optimal Bayes One Stage', 'Location', 'north','AutoUpdate', 'Off' ) ;
    xlim( [ plot_lower, plot_upper] ) ;
    lohivec = [ 0.9, 1 ] ;
    UtilPlotABCD( basic, advanced, mat, ishoriz, ist0, lohivec ); 
    UtilStdizeFigure(fignum,advanced);
    UtilSaveFigEpsPdf(fignum,dirname, strcat('pr_select_best',graphicsuffix),'-r600');

    % Plotting the expected number seen (Figure 4 of the manuscript)

    fignum = fignum + 1; figure(fignum);
    ishoriz = true; % set to true if muvec is on horizontal axis, false if on vertical axis
    ist0 = false ;
    hold off ; 
    %Opt_Fixed = dlmread( 'OptFixed.txt' ) ;
    %Opt_Fixed = Opt_Fixed' ; 
    %Opt_Fixed: Using a different vector set ... see plot below
    %ENumSeenAve = dlmread( 'ENumSeenAve.txt' ) ; 
    ENumSeenAve = mat.outBayes.ENumSeenAve ; 
    TMax = ones( 1, length( advanced.simFreqDeltaVec ) ) * basic.TMax ; 
    %plot(  advanced.simFreqDeltaVec, ENumSeenAve', '--k', advanced.simFreqDeltaVec, TMax, '-.k', advanced.simFreqDeltaVec, Opt_Fixed, 'k' ) ;  
    plot(  advanced.simFreqDeltaVec, ENumSeenAve', '--k', advanced.simFreqDeltaVec, TMax, '-.k',  mat.muvec,mat.OptOneShotLength, 'k' ,'LineWidth',linewid ) ;  
    xlabel( 'Prior mean', 'Fontsize', advanced.bigfontsize ) ; 
    ylabel( 'Expected number of samples' , 'Fontsize', advanced.bigfontsize) ; 
    legend( 'Optimal Bayes Sequential',  'Fixed', 'Optimal Bayes One Stage', 'Location', 'north','AutoUpdate', 'Off' ) ;
    hold on ; 
    ylim( [ 0, 2100 ] ) ; 
    xlim( [ plot_lower,plot_upper] ) ;
    lohivec = [ 0, 2100 ] ;
    UtilPlotABCD( basic, advanced, mat, ishoriz, ist0, lohivec ); 
    hold on ; 
    UtilStdizeFigure(fignum,advanced);
    UtilSaveFigEpsPdf(fignum,dirname, strcat('E_num_seen',graphicsuffix),'-r600');

    fignum = fignum + 1; figure(fignum);
    hold off ; 
    ishoriz = true; % set to true if muvec is on horizontal axis, false if on vertical axis
    ist0 = false ;
    %fignum = fignum + 1; figure(fignum);
    plot(advanced.simFreqDeltaVec, mat.outBayes.PrReversalAve, '--k' ,'LineWidth',linewid);
    hold on ; 
    xlim( [ plot_lower, plot_upper] ) ;
    hold on ; 
    xlabel( 'Prior mean', 'Fontsize', advanced.bigfontsize ) ; 
    ylabel( 'Prob(reversing decision with pipeline data)', 'Fontsize', advanced.bigfontsize ) ; 
    %legend( 'Optimal Bayes Sequential',  'Fixed', 'Optimal Bayes One Stage', 'Location', 'north' ) ;
    hold on ; 
    lohivec = [ 0.0, min(1.0, 1.2*max(mat.outBayes.PrReversalAve)) ] ;
    UtilPlotABCD( basic, advanced, mat, ishoriz, ist0, lohivec ); 
    UtilStdizeFigure(fignum,advanced);
    UtilSaveFigEpsPdf(fignum,dirname, strcat('pr_reversal',graphicsuffix),'-r600');

    %UtilStdizeFigure( fignum, advanced )
    %%legend('our','all','trial','fixed')
    %hold off
    %fignum = fignum + 1; figure(fignum);
    %plot( advanced.simFreqDeltaVec, mat.outBayes.RealizedRewardAve, advanced.simFreqDeltaVec, mat.outBayes.RealizedRewardAllAve, advanced.simFreqDeltaVec, mat.outBayes.RealizedRewardTrialAve, advanced.simFreqDeltaVec, mat.outBayes.RealizedRewardFixedAve)
    % xlim([mat.muvec(mat.indlow1) mat.muvec(mat.indup1)]);
    %xlim([mat.muvec(mat.Threshpoint(2)) mat.muvec(mat.Threshpoint(1))]);
    %legend('our','all','trial','fixed')
    %UtilStdizeFigure( fignum, advanced )

    %fignum = fignum + 1; figure(fignum);
    %plot (advanced.simFreqDeltaVec, mat.outBayes.RealizedRewardAve, '-', ...
    %     advanced.simFreqDeltaVec, mat.outBayes.RealizedRewardAve + mat.outBayes.RealizedRewardSE, '-.', ...
    %     advanced.simFreqDeltaVec, mat.outBayes.RealizedRewardTrialAve, '--', ...
    %     advanced.simFreqDeltaVec, mat.outBayes.RealizedRewardTrialAve + mat.outBayes.RealizedRewardTrialSE,'.' , ...
    %     advanced.simFreqDeltaVec, mat.outBayes.RealizedRewardAve - mat.outBayes.RealizedRewardSE, '-.', ...
    %     advanced.simFreqDeltaVec, mat.outBayes.RealizedRewardTrialAve - mat.outBayes.RealizedRewardTrialSE,'.' )
    %xlim([mat.muvec(mat.Threshpoint(2)) mat.muvec(mat.Threshpoint(1))]);
    %title('diff for ours and trial');
    %legend('our',sprintf('our +/- %3.2f SE',zval),'trial',sprintf('trial +/- %3.2f SE',zval));
    %UtilStdizeFigure( fignum, advanced ) ;
end %if colorflag

    figout = fignum;

end