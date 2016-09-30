function [figout, mat] = DoSectionFivePlotsStentsUnk(fignum,basic,advanced,PathReps,ProductionReps,graphicsuffix)

dirname = 'Figure'; % directory for putting figures
linewid = 2;


[basic, advanced, rval, msgs] = DelayInputValidator( basic, advanced );
if ~rval, msgs, end; 
advanced.dt

advanced.PLOTSIMS = false ;

[~, mat] = DelayCurvesRecur(basic, advanced);
[mat] = DelayStageOne(basic, advanced, mat ); 


%% Plotting the optimal stopping rule
fignum = fignum + 1; figure(fignum); % always increment fignum and create a new figure
hold off ; 
%plot( paths( 1:outtime6, 1 )', paths( 1 : outtime6, 7 )', 'k', paths( outtime6 + 1:decidetime6, 1 )', paths( outtime6 + 1:decidetime6, 7 )', ':k' ) ;  
%legend( 'Actual Stage II (observing and sampling)', 'Actual Stage III (observing only)','Location','South') ; 
%hold on ; 
% plot Stage 1 boundary if there exists a Stage 1
plot (mat.tvec, mat.bndupper,'--k',mat.tvec ,mat.bndlower, '--k' ,'LineWidth',linewid);
hold on ; 
if abs( mat.Threshpoint(2) - mat.Threshpoint(4) ) > 1
    scatter (mat.bestsvec(mat.Threshpoint(2):mat.Threshpoint(4)) + basic.t0, mat.muvec(mat.Threshpoint(2):mat.Threshpoint(4)),'k');
    hold on ; 
end
% plot upper Stage 1 boundary if there exists a Stage 1
if abs( mat.Threshpoint(3) - mat.Threshpoint(1)) > 1
    scatter (mat.bestsvec(mat.Threshpoint(3):mat.Threshpoint(1)) + basic.t0, mat.muvec(mat.Threshpoint(3):mat.Threshpoint(1)),'k');
    hold on ; 
end 
ylabel('Prior / Posterior Mean', 'Fontsize', advanced.bigfontsize);
xlabel('n_0 + t', 'Fontsize', advanced.bigfontsize);
ishoriz = false; % set to true if muvec is on horizontal axis, false if on vertical axis
plot_upper = mat.muvec(mat.Threshpoint(1)) + min(mat.muvec(mat.Threshpoint(1))/4,- mat.muvec(mat.Threshpoint(2))/4) ; % set upper limit of plot region ( + /  
plot_lower =  mat.muvec(mat.Threshpoint(2)) - min(-mat.muvec(mat.Threshpoint(2))/4,mat.muvec(mat.Threshpoint(1))/4)  ; % set lower limit of plot region
%advanced.smallfontsize = 11;
plot_lower_tmp=plot_lower-500;
plot_upper_tmp=plot_upper+500;
ylim( [ plot_lower_tmp, plot_upper_tmp] ) ;
xlim( [ 0, 2500] ) ;
line( [ basic.t0 + basic.tau, basic.t0 + basic.tau ], [ plot_lower_tmp, plot_upper_tmp ], 'Linestyle', ':', 'Color', 'k' ) ;
line( [ basic.t0 + basic.TMax, basic.t0 + basic.TMax ], [ plot_lower_tmp, plot_upper_tmp ], 'Linestyle', ':', 'Color', 'k' ) ;
line( [ basic.t0, basic.t0 ], [ plot_lower_tmp, plot_upper_tmp ], 'Linestyle', ':', 'Color', 'k' ) ;
if advanced.DOPDE
    line( [basic.t0 + basic.tau, basic.t0 + basic.tau ], [ mat.muvec(mat.Threshpoint(4)),mat.bndlower(1) ], 'Linestyle', '--', 'Color', 'k', 'LineWidth',linewid*1.1) ;
end
%line( [ basic.t0 + basic.tau + basic.TMax, basic.t0 + basic.tau + basic.TMax ], [ plot_lower_tmp, plot_upper_tmp ], 'Linestyle', ':', 'Color', 'k' ) ;
tmpfontsize = round(advanced.smallfontsize * 1.1);
text(  basic.t0 + 2.5,  plot_lower_tmp + 650, 'n_0', 'Fontsize', tmpfontsize) ;
text(  basic.t0 + basic.TMax - 180,  plot_lower_tmp + 650  , 'n_0 + T_{max}', 'Fontsize', tmpfontsize) ;
text(  basic.t0 + basic.tau + 2.5,  plot_lower_tmp + 650 , 'n_0 + \tau', 'Fontsize', tmpfontsize) ;
%text(  basic.t0 + basic.tau + basic.TMax - 550,  plot_lower_tmp + 600 , 'n_0 + \tau + T_{max}', 'Fontsize', tmpfontsize) ;
if length(graphicsuffix) == 0
    text(  ( basic.t0 + basic.tau ) / 5.0 , (mat.muvec(mat.Threshpoint(1)) + plot_upper_tmp)/2 , 'No trial', 'Fontsize', tmpfontsize); 
    text(  ( basic.t0 + basic.tau ) / 5.0,  (mat.muvec(mat.Threshpoint(2)) + plot_lower_tmp)/2 , 'No trial', 'Fontsize', tmpfontsize) ;
    %annotation( 'textarrow', [ 0.19, 0.16 ], [ 0.92, 0.97 ], 'String', 'No trial', 'Fontsize', tmpfontsize ) ;
    %annotation( 'textarrow', [ 0.19, 0.16 ], [ 0.17, 0.135 ], 'String', 'No trial', 'Fontsize', tmpfontsize ) ;
    text(  ( basic.t0 + basic.tau ) / 7.0 , -5000 , 'Sequential trial', 'Fontsize', tmpfontsize  ) ;
    text(  ( basic.t0 + basic.tau ) / 7.0 , -6500 , 'recruitment', 'Fontsize', tmpfontsize  ) ;
    %text(  ( basic.t0 + basic.tau ) / 6.2 , 5700 , 'Maximum extent of stage I', 'Fontsize', tmpfontsize) ;
    %text(  ( basic.t0 + basic.tau ) * 1.05,  5700 , 'Maximum extent of stage II', 'Fontsize', tmpfontsize) ;
    %text(  ( basic.t0 + basic.tau ) * 1.95,  5700 , 'Maximum extent of stage III', 'Fontsize', tmpfontsize) ;
end
ist0 = true ; % set to true if ABCD to be plotted at t0, false if to be plotted on horizontal axis
UtilPlotABCD( basic, advanced, mat, ishoriz, ist0 );  
%plot( paths( 1:outtime1, 1 )', paths( 1 : outtime1, 2 )', 'k', paths( outtime1 + 1:decidetime1, 1 )', paths( outtime1 + 1:decidetime1, 2 )', ':k' ) ; 
%plot( paths( 1:outtime2, 1 )', paths( 1 : outtime2, 3 )', 'k', paths( outtime2 + 1:decidetime2, 1 )', paths( outtime2 + 1:decidetime2, 3 )', ':k' ) ; 
%plot( paths( 1:outtime3, 1 )', paths( 1 : outtime3, 4 )', 'k', paths( outtime3 + 1:decidetime3, 1 )', paths( outtime3 + 1:decidetime3, 4 )', ':k' ) ; 
%plot( paths( 1:outtime7, 1 )', paths( 1 : outtime7, 8 )', 'k', paths( outtime7 + 1:decidetime7, 1 )', paths( outtime7 + 1:decidetime7, 8 )', ':k' ) ; 
%plot( paths( 1:outtime10, 1 )', paths( 1 : outtime10, 11 )', 'k', paths( outtime10 + 1:decidetime10, 1 )', paths( outtime10 + 1:decidetime10, 11 )', ':k' ) ; 
%text(  paths( decidetime3 ), outtimes( 3, 5) , '*1', 'Fontsize', advanced.smallfontsize) ;
%text(  paths( decidetime6 ), outtimes( 6, 5 ) , '*2', 'Fontsize', advanced.smallfontsize) ;
%text(  paths( decidetime7 -2 ), outtimes( 7, 5 ) , '*4', 'Fontsize', advanced.smallfontsize) ;
%text(  paths( decidetime7 -4 ), paths( decidetime7 , 8 ) - 10 , '4', 'Fontsize', advanced.smallfontsize) ;
%text(  paths( decidetime10 ), outtimes( 10, 5) , '*3', 'Fontsize', advanced.smallfontsize) ;
%xlim( [ 0 basic.TMax + basic.t0 + basic.tau ] ); % Baseline
%xlim( [ 0 3100 ] ); % Comparator
UtilStdizeFigure(fignum,advanced,true);
dirname = 'Figure' ;
UtilSaveFigEpsPdf(fignum,dirname,strcat('pathsStents',graphicsuffix),'-r600');
%savefig paths ; 
%print -r600 -deps paths ; 
%print -r600 -dpdf paths ; 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% For Other figures, use estimates based on a lot of sample
%%%%%% replications. so, the simulations need to be rerun.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

advanced.simNumReps = ProductionReps;  % set to specific number of replications required for this specific experiment.

%advanced.simFreqDeltaVec = (plot_lower:(plot_upper-plot_lower)/advanced.numinsimFreqDeltaVec:plot_upper);% P : EXPERIMENT

advanced.numinsimFreqDeltaVec=300;  % can change the number of points in freqdeltavec. the 11/3 in the next line is to get \pm 11/3 standard errors in width for vector's range
%advanced.simFreqDeltaVec = (basic.sigma/sqrt(basic.t0)) * (11/3) * (-ceil(advanced.numinsimFreqDeltaVec/2):ceil(advanced.numinsimFreqDeltaVec/2)) / ceil(advanced.numinsimFreqDeltaVec/2);
advanced.simFreqDeltaVec = (plot_lower:(plot_upper-plot_lower)/advanced.numinsimFreqDeltaVec:plot_upper);% P : EXPERIMENT

[ fignum, mat ] = DelaySimOverview( fignum, basic, advanced, mat ) ;
fignum = DelaySimOutput( fignum, basic, advanced, mat );  

% calculating the gain in terms of "realized reward" from using our method.
% MF 27/03: diff_trial is used in Section 5 because we have a reference
% experiment. diff_all dropped because is dealt with in Section 4.
% SEs are also included to see whether they make a difference
% when running 15,000 replications. These SE estimates are larger than they
% really are if CRN is used.
%diff_all = mat.outBayes.RealizedRewardAve - mat.outBayes.RealizedRewardAllAve; % difference between the "realized reward" of our methods and that of taking always a sample size equal to T_max
diff_trial = mat.outBayes.RealizedRewardAve - mat.outBayes.RealizedRewardTrialAve; % difference between the "realized reward" of our methods and that of taking always a sample size equal to the sample size in the trial we took the data from
diff_fixed = mat.outBayes.RealizedRewardAve - mat.outBayes.RealizedRewardFixedAve; % difference between the "realized reward" of our methods and that of taking always an optimal fixed sample size (i.e. ignoring Stage 2) 
diff_zero = mat.outBayes.RealizedRewardAve  - mat.outBayes.RealizedRewardAve  ; % for plotting a line at 0
%diff_all_approxSEifindep = sqrt(mat.outBayes.RealizedRewardSE.^2 + mat.outBayes.RealizedRewardAllSE.^2);
diff_trial_approxSEifindep = sqrt(mat.outBayes.RealizedRewardSE.^2 + mat.outBayes.RealizedRewardTrialSE.^2);
diff_fixed_approxSEifindep = sqrt(mat.outBayes.RealizedRewardSE.^2 + mat.outBayes.RealizedRewardFixedSE.^2);

%% MF 08/03: plot the differences in the expected rewards
fignum = fignum + 1; figure(fignum); % always increment fignum and create a new figure
hold off ; 
% MF 04/03: Plotting the gain in terms of realised reward using our method 
ishoriz = true ; % set to true because muvec is on horizontal axis
ist0 = false ; % set to true if ABCD to be plotted at t0, false if to be plotted on horizontal axis
% Plot the differences for fixed and all
plot(  advanced.simFreqDeltaVec, diff_fixed, 'k', advanced.simFreqDeltaVec, diff_trial, '-.k','LineWidth',linewid ) ;  
%plot_upper = 15000 ; % set upper limit of plot region
%plot_lower = -15000  ; % set lower limit of plot region
xlim( [ plot_lower, plot_upper] ) ;
ylim( [ -75000000/10 , 75000000 ] ) ;  
xlabel( 'Prior mean', 'Fontsize', advanced.bigfontsize ) ; 
ylabel( 'Difference between expected rewards', 'Fontsize', advanced.bigfontsize) ; 
legend( 'Optimal Bayes One Stage', 'Fixed (= trial sample size)', 'Location', 'northwest' ) ; 
hold on ; 
plot( advanced.simFreqDeltaVec, diff_zero, ':k','LineWidth',linewid ) ; 
% MF 08/03: Steve's code for SEs: these correctly compute SE estimate when
% output is run with CRN across experiments.
if advanced.keepAllOutput   % if all the simulation output was kept, try to produce some CIs for differences based on the use of CRN (if it was used)
    difftrialvecAVE =zeros(length(advanced.simFreqDeltaVec));
    difftrialvecSTD =difftrialvecAVE;
    %diffallvecAVE = difftrialvecAVE;
    %diffallvecSTD = difftrialvecAVE;
    difffixedvecAVE =  difftrialvecAVE;
    difffixedvecSTD = difftrialvecAVE;
    for i=1:length(advanced.simFreqDeltaVec)
        %diffall = mat.simBayesOut(i).RealizedReward - mat.simBayesOut(i).RealizedRewardAll;
        %diffallvecAVE(i) = mean(diffall);
        %diffallvecSTD(i) = std(diffall)/sqrt(length(diffall));
        difffixed = mat.simBayesOut(i).RealizedReward - mat.simBayesOut(i).RealizedRewardFixed;
        difffixedvecAVE(i) = mean(difffixed);
        difffixedvecSTD(i) = std(difffixed)/sqrt(length(difffixed));
        difftrial = mat.simBayesOut(i).RealizedReward - mat.simBayesOut(i).RealizedRewardTrial;
        difftrialvecAVE(i) = mean(difftrial);
        difftrialvecSTD(i) = std(difftrial)/sqrt(length(difftrial));
    end
    %fignum = fignum + 1; figure(fignum);
    %plot( advanced.simFreqDeltaVec, difffixedvecAVE + zval * difffixedvecSTD, ':k',advanced.simFreqDeltaVec, difffixedvecAVE - zval * difffixedvecSTD, ':k' ) ;
    %hold on ; 
    %title(sprintf('Difference in Outcomes between Sequential and All Pairs Trials\n%s',advanced.titlestring),'FontSize',advanced.bigfontsize,'FontName',advanced.fontname)
    %xlabel('posterior mean per patient','FontSize',advanced.smallfontsize,'FontName',advanced.fontname);
    %ylabel('population benefit from sequential trial','FontSize',advanced.smallfontsize,'FontName',advanced.fontname);
    %legend('difference',sprintf('difference +/- %3.2f SE',zval));
    %UtilStdizeFigure( fignum, advanced ); 
    %fignum = fignum + 1; figure(fignum);
    %plot ( advanced.simFreqDeltaVec, difftrialvecAVE + zval * difftrialvecSTD, ':k', advanced.simFreqDeltaVec, difftrialvecAVE - zval * difftrialvecSTD, ':k') ;
    %hold on ; 
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
lohivec = [ 0, 75000000 ] ;
UtilPlotABCD( basic, advanced, mat, ishoriz, ist0, lohivec );  
UtilStdizeFigure(fignum,advanced);
dirname = 'Figure' ;
UtilSaveFigEpsPdf(fignum,dirname, strcat('diffsStents',graphicsuffix),'-r600');

% SC's commands for plotting the differences in realized rewards (would expect these
% differences to be always positive over the range of prior mean
% considered)
%fignum = fignum + 1; figure(fignum);
%plot (advanced.simFreqDeltaVec, diff_trial, '-', advanced.simFreqDeltaVec, diff_fixed,'-')
%xlim([mat.muvec(mat.Threshpoint(2)) mat.muvec(mat.Threshpoint(1))]);
%legend('trial','fixed')
%UtilStdizeFigure( fignum, advanced )
 
%% Plotting the probability of making the correct decision 
fignum = fignum + 1; figure(fignum);
hold off ; 
ishoriz = true; % set to true if muvec is on horizontal axis, false if on vertical axis
ist0 = false ;
%fignum = fignum + 1; figure(fignum);
plot(advanced.simFreqDeltaVec, mat.outBayes.PrTrueBestSelectedAve, '--k', advanced.simFreqDeltaVec, mat.outBayes.PrTrueBestSelectedTrialAve, '-.k', advanced.simFreqDeltaVec, mat.outBayes.PrTrueBestSelectedFixedAve, '-k','LineWidth',linewid ); 
hold on ; 
%plot_upper = 15000 ; % set upper limit of plot region ( + / 
%plot_lower = -15000  ; % set lower limit of plot region
xlim( [ plot_lower, plot_upper] ) ;
ylim( [ 0.9, 1 ] ) ; 
hold on ; 
xlabel( 'Prior mean', 'Fontsize', advanced.bigfontsize ) ; 
ylabel( 'Proportion of correct decisions', 'Fontsize', advanced.bigfontsize ) ; 
legend( 'Optimal Bayes Sequential',  'Fixed (= trial sample size)', 'Optimal Bayes One Stage', 'Position',[220,65,10,10] ) ;
hold on ; 
lohivec = [ 0.9, 1 ] ;
UtilPlotABCD( basic, advanced, mat, ishoriz, ist0, lohivec ); 
UtilStdizeFigure(fignum,advanced);
UtilSaveFigEpsPdf(fignum,dirname, strcat('pr_select_bestStents',graphicsuffix),'-r600');

% Plotting the expected number seen
fignum = fignum + 1; figure(fignum);
ishoriz = true; % set to true if muvec is on horizontal axis, false if on vertical axis
ist0 = false ;
hold off ; 
ENumSeenAve = mat.outBayes.ENumSeenAve ; 
Optfix_4 =  [ mat.simBayesOut.NoStatsOpt_fix4 ] ;
p_basic.numpairs = basic.numpairs .* ones( 1, length( advanced.simFreqDeltaVec ) ) ; 
plot(  advanced.simFreqDeltaVec, ENumSeenAve', '--k', advanced.simFreqDeltaVec, p_basic.numpairs, '-.k', advanced.simFreqDeltaVec, Optfix_4, 'k','LineWidth',linewid ) ; 
hold on ; 
legend( 'Optimal Bayes Sequential',  'Fixed (= trial sample size)', 'Optimal Bayes One Stage', 'Location', 'northeast' ) ;
hold on ; 
ylim( [ 0, max([max(ENumSeenAve) max(p_basic.numpairs)  max(Optfix_4)])*1.1 ] ) ; 
%plot_upper = 15000 ; % set upper limit of plot region ( + / 
%plot_lower = -15000  ; % set lower limit of plot region
xlim( [ plot_lower, plot_upper] ) ;
hold on ; 
xlabel( 'Prior mean', 'Fontsize', advanced.bigfontsize ) ; 
ylabel( 'Expected sample size' , 'Fontsize', advanced.bigfontsize) ; 
hold on ; 
lohivec = [ 0, max([max(ENumSeenAve) max(p_basic.numpairs)  max(Optfix_4)])*1.1 ] ;
UtilPlotABCD( basic, advanced, mat, ishoriz, ist0, lohivec ); 
hold on ; 
UtilStdizeFigure(fignum,advanced);
UtilSaveFigEpsPdf(fignum,dirname, strcat('E_num_seenStents',graphicsuffix),'-r600');

fignum = fignum + 1; figure(fignum);
hold off ; 
ishoriz = true; % set to true if muvec is on horizontal axis, false if on vertical axis
ist0 = false ;
%fignum = fignum + 1; figure(fignum);
plot(advanced.simFreqDeltaVec, mat.outBayes.PrReversalAve, '--k','LineWidth',linewid);
hold on ; 
%plot_upper = 15000 ; % set upper limit of plot region ( + / 
%plot_lower = -15000  ; % set lower limit of plot region
xlim( [ plot_lower, plot_upper] ) ;
ylim( [ 0, 0.05] ) ;
hold on ; 
xlabel( 'Prior mean', 'Fontsize', advanced.bigfontsize ) ; 
ylabel( 'Probability of reversing decision with pipeline data', 'Fontsize', advanced.bigfontsize ) ; 
%legend( 'Optimal Bayes Sequential',  'Fixed', 'Optimal Bayes One Stage', 'Location', 'north' ) ;
hold on ; 
lohivec = [ 0.0, min(1.0, 1.2*max(mat.outBayes.PrReversalAve)) ] ;
lohivec = [ 0.0, min(1.0, 0.07) ] ;
UtilPlotABCD( basic, advanced, mat, ishoriz, ist0, lohivec ); 
UtilStdizeFigure(fignum,advanced);
UtilSaveFigEpsPdf(fignum,dirname, strcat('pr_reversalStents',graphicsuffix),'-r600');


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

    figout = fignum;


