function [figout] = DelaySimOutput( fignum, basic, advanced, mat )
%DelaySimOutput: Takes output mat structure generated by DelaySimComputer
%and generates several plots and statistical outputs, for the delay
%sequential sampling project of Martin, Paolo and Steve. 
%
% Typically this function will not be caled directly, but will be called
% via DelaySimOverview().
%
% The output is based on vectors of the following simout structure, one
% vector for Bayesian statistics, one vector for frequentist statistics.
% Each of these vectors is of lenth advanced.simFreqDeltaVec. The Bayesian
% vector is called mat.simBayesOut and represent realizations of trials
% which have means whose realizations are drawn from the prior distribution
% of the unknown mean with mean advanced.simFreqDeltaVec(i) and variance
% basic.sigma/sqrt(basic.t0) (that is, the prior variance together with a
% variety of different means). The frequentisti results assume the mean is
% actually equal to advanced.simFreqDeltaVec(i). Thus, this routine 
% provides both Bayesian output and frequentist outputs (such as power 
% curves).
%
% 2014 Apr 27: Created by Steve
%
if ~isfield(mat,'outBayes')
    warning('DelaySimOutput: outBayes field not found in mat, need to call DelaySimOverview first');
else
	codePDE = '--xk';codePDEFr = '--om';
    
    if advanced.saveplot % If files to be saved, write information about the version of code used.
        UtilSaveVersionFile( advanced.dirstring );
    end

    fields=fieldnames(mat.outBayes);
    for i=1:numel(fields)
        a = strfind( fields{i}, 'Ave' );
        if ~isempty(a)          % if 'Ave' is in name, 
            basename = fields{i}(1:(a-1));
            b = strfind( basename, 'All' );  % naming: [Pr]name{All}[Ave/SE]. All is optional.
            if isempty(b)       % ... and All is not... then do the plots. Plot the sequential data (and one-stage if it is available)
                fignum=fignum+1; figure(fignum); hold off;
                
                avevec = mat.outBayes.(fields{i});
                sevec = mat.outBayes.(strcat(basename,'SE'));
                avevecfreq = mat.outFreq.(fields{i});
                sevecfreq = mat.outFreq.(strcat(basename,'SE'));
                doOneStagePlot = isfield(mat.outBayes,[basename 'AllAve']);  % check if there is an analogous field for the other data
                if doOneStagePlot     
                    avevec2 = mat.outBayes.(strcat(basename,'AllAve'));
                    sevec2 = mat.outBayes.(strcat(basename,'AllSE'));
                    avevecfreq2 = mat.outFreq.(strcat(basename,'AllAve'));
                    sevecfreq2 = mat.outFreq.(strcat(basename,'AllSE'));
                end

                plot(advanced.simFreqDeltaVec,avevec,'-xr',advanced.simFreqDeltaVec,avevecfreq,'-vg','LineWidth',2);
                hold on
                legendvec = {'Sequential, MC Bayes prior mean','Sequential, MC freq. true mean'};
                if doOneStagePlot
                    plot(advanced.simFreqDeltaVec,avevec2,'-+b',advanced.simFreqDeltaVec,avevecfreq2,'-^c','LineWidth',1.5);
                    legendvec = [ legendvec 'One stage, MC Bayes prior mean' 'One stage, MC freq. true mean'];
                end
                %legendvec = [legendvec {legendval}];
                if advanced.PDEWithSimPlots
                    switch basename
                      case 'ENumSeen'   %% mean number of samples as function of prior mean (works).
                          % Note: The simulation runs even for means outside of
                          % the normal stopping time - this serves as a test of
                          % the PDE approximation - indeed one should be
                          % willing to simulate a bit more than the PDE would
                          % suggest at extreme values of the prior mean.
                        tauplusstageII = basic.tau + mat.ENumSamps;
                        tauplusstageII(mat.ENumSamps == 0) = 0;
                        plot(mat.muvec(:),mat.optENumSamps(:),codePDE);
                        legendvec = [legendvec 'Optimal (Stages I, II, III)'];
                        plot(mat.muvec(:),tauplusstageII,'x');
                        legendvec = [legendvec 'Seq., PDE Bayes (Stage II)'];
    %                  case 'DecisionMean'      % has slope of less than 45 degrees. Regret?
    %                  case 'DecisionMeanAll' % has slope of 45 degrees
    %                  case 'RealizedMean'
                      case 'RealizedReward' % case 'RealizedRewardAll'
                        plot(mat.muvec(:),mat.B0hat(:),codePDE);
                        legendvec = [legendvec 'Sequential, PDE Bayes'];
%                        if doOneStagePlot  % BayesOneshotExpectedReward is
%                        really for 'ExpectedReward', not 'RealizedReward'
%                            plot(mat.muvec(:),mat.BayesOneshotExpectedReward(:),codePDEFr);
%                            legendvec = [legendvec 'One stage, Bayes (analytical)'];
%                        end
                      case 'ExpectedReward'        % FIX: PDE solution is 'higher' than the simulaiton solution. Discounting effect and fact taht sims run even if stopping is initially called for??
                        plot(mat.muvec(:),mat.B0hat(:),codePDE);
                        legendvec = [legendvec 'Sequential, PDE Bayes'];
                        if doOneStagePlot
                            plot(mat.muvec(:),mat.BayesOneshotExpectedReward(:),codePDEFr);
                            legendvec = [legendvec 'One stage, Bayes (analytical)'];
                        end
                      case 'PrHitTop'   % Prob (top stopping boundary is hit before bottom boundary)(works).
                        plot(mat.muvec(:),mat.Puppermat(:,1),codePDE);
                        legendvec = [legendvec 'Sequential, PDE Bayes'];
                        tmp = axis; tmp(3)=0.0; tmp(4)=1.01; axis(tmp);
                      case 'PrNewSelected'          % Prob (new alternative is selected) (works).
                        plot(mat.muvec(:),mat.Pnewmat(:,1),codePDE);
                        legendvec = [legendvec 'Sequential, PDE Bayes'];
                        tmp = axis; tmp(3)=0.0; tmp(4)=1.01; axis(tmp);
%                  case 'PrNewSelectedAll'
                      case 'PrTrueBestSelected'     % PCS, given sequential procedure 
                        plot(mat.muvec(:),mat.Bayespcsvec(:),codePDE);
                        legendvec = [legendvec 'Sequential, PDE Bayes'];
                        if doOneStagePlot
                            plot(mat.muvec(:),mat.OneshotPCS(:),codePDEFr);
                            legendvec = [legendvec 'One stage, Bayes (analytical)'];
                        end
                        tmp = axis; tmp(3)=0.4; tmp(4)=1.01; axis(tmp);
                    end % switch
                end %PDEwithSimPlots

                % plot the legend
%                legend(legendvec,'Location','NorthWest');  % use this command for versions prior to Matlab R2017a
                legend(legendvec,'Location','NorthWest','AutoUpdate','off'); % FIX: Use AutoUpdate off for Matlab versions R2017a and later
                % then the error bounds
                plot(advanced.simFreqDeltaVec,avevec+sevec,'-r',advanced.simFreqDeltaVec,avevec-sevec,'-r','LineWidth',0.5);
                plot(advanced.simFreqDeltaVec,avevecfreq+sevecfreq,'-.g',advanced.simFreqDeltaVec,avevecfreq-sevecfreq,'-.g','LineWidth',0.5);
                if doOneStagePlot     % of the two-stage data is available, also plot it
                    plot(advanced.simFreqDeltaVec,avevec2+sevec2,'-b',advanced.simFreqDeltaVec,avevec2-sevec2,'-b','LineWidth',0.5);
                    plot(advanced.simFreqDeltaVec,avevecfreq2+sevecfreq2,'-.c',advanced.simFreqDeltaVec,avevecfreq2-sevecfreq2,'-.c','LineWidth',0.5);
                end

                % do some formatting and save if necessary
                UtilStdizeFigure( fignum, advanced );
                title(sprintf('%s\n%s (sample mean +/- std err)',advanced.titlestring,basename),'FontSize',advanced.bigfontsize,'FontName',advanced.fontname);
                xlabel('Mean reward (prior or true)','FontSize',advanced.smallfontsize,'FontName',advanced.fontname); 
                ylabel(basename,'FontSize',advanced.smallfontsize,'FontName',advanced.fontname);
                if strncmp('Pr',fields{i},2)   % if field is for a probability, set y-axis to be 0-1
                    plot([min(mat.bndlower) max(mat.bndupper)],[1 1],'xk');
                else
                    plot([min(mat.bndlower) max(mat.bndupper)],[0 0],'xk');
                end
                %axis('square');
                if advanced.saveplot
                    texttodifferentiate = sprintf('Sim%s',basename);
                    UtilSaveFigFile(fignum, advanced.dirstring, advanced.filestring, texttodifferentiate, advanced.graphicextension);
                end
            end %isempty(b)
        end %~isempty(a)
    end % for i = 1:numfields
end

%fignum=fignum+1; figure(fignum);
% FIX: PUT SOME STUFF ABOUT REGRET. VALIDATE THAT OUTPUT OF PLOTS MAKES
% SENSE.

figout = fignum;

end

