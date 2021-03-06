function [ figout, matout ] = DelaySimOverview( fignum, basic, advanced, mat )
%DelaySimOverview: Call this function to drive the generation of Monte 
% Carlo simulations for delay sequential sampling project of Martin, Paolo
% and Steve. Assumes that the basic and advanced structures have been
% created and validated, and that the mat structure has been filled in by
% running the dynamic programming analysis implemented by DelayCurvesRecur
% and DelayStageOne (or as generated by the TestDelayIterate function).
%
% IMPORTANT: Assumes sampling distribution is gaussian with known sampling
% variance and unknown mean. See DelaySimUnkOverview.m for the case of
% unknown sampling variances.
%
% If advanced.simNumReps is at least 2 (so that std dev estimates can be
% obtained), then the routine will run the simulations. The simulations
% generate advanced.simNumReps sample paths for each of the following
% scenarios:
%   - Bayesian: sample unknown mean from prior distribution, generate sample
%   paths accordingly.
%   - Frequentist: for each 'true mean' in vector advanced.simFreqDeltaVec,
%   generate sample paths which can be used for power calculations
%
% On output, creates a new mat structure with all of the attributes of the
% mat structure input, plus several structures which contain all of the
% simulation output statistics for the Bayesian and each frequentist run.
%
% Called by: 
%   - DelayDriver, 
%   - DelayDriverUnknown, 
%   - DelayPaperExperiments,
%   - TestDelayIterate
%
% 2014 Apr 27: Created by Steve
%
MYEPS = 10e-7;
matout = mat;
NREPS = advanced.simNumReps;                % number of replications to run

if NREPS < 2
    warning('DelaySimOverview should be called with advanced.simNumReps >= 2 to get Monte Carlo Results');
else
%	if advanced.UnkVariance
%            [ fignum, mat ] = DelaySimUnkOverview( fignum, basic, advanced, mat );
%    else
        origPLOTSIMS = advanced.PLOTSIMS;
        
        if advanced.CRN   % if CRN desired across multiple calls to DelaySimOverview, such as for doing across mat structure in a matvec from TestDelayIterate, reset the stream structure
            stream = RandStream('twister','Seed',int32(advanced.CRNAcrossExperiment)); % need int32() in case advanced.CRNAcrossExperiment is thought to be a logical.
            RandStream.setGlobalStream(stream);
        end

        normalnoise = normrnd(0,basic.sigma,basic.TMax,NREPS);  % compute the sampling noise here, to implement CRN across power tests and Bayesian means
        if exist('lhsnorm','file')              % if lhsnorm is available (part of statistics toolbox; not usual installation) use latin hypercube for bayesian means
            % Note: lhsnorm expects a variance, whereas normrn expects a standard deviation!
            bayesmeandelta = lhsnorm(0,(basic.sigma)^2/basic.t0,NREPS);    % compute deltas for realized Bayes means
        else                                    % otherwise usu the usual random samples from posterior distribution
            bayesmeandelta = normrnd(0,basic.sigma/sqrt(basic.t0),NREPS,1);        % compute deltas for realized Bayes means
        end

        % Do Bayes stuff first
        for i=1:length(advanced.simFreqDeltaVec)
            advanced.simTmpFlag = true;        % this is a temporary flag for use in sim computer, true for bayesian stopping bound, false for frequentist
            delta = advanced.simFreqDeltaVec(i);
            advanced.PLOTSIMS = ( (abs(abs(delta) - basic.sigma/sqrt(basic.t0)) < MYEPS) | (abs(delta) < MYEPS) ) & origPLOTSIMS;
            if (~advanced.CRN) && (i ~= 1)    % if CRN across bayes and frequent values of delta desired, keep same noise, otherwise get independently sampled noise
                normalnoise = normrnd(0,basic.sigma,basic.TMax,NREPS);  % compute the sampling noise here, to implement CRN across power tests and Bayesian means
            end
            if (~advanced.CRNAcrossBayesMu) && (i ~= 1)    % if not CRN for means get new set of means
                if exist('lhsnorm','file')              % if lhsnorm is available (part of statistics toolbox; not usual installation) use latin hypercube for bayesian means
                    % Note: lhsnorm expects a variance, whereas normrn expects a standard deviation!
                    bayesmeandelta = lhsnorm(0,(basic.sigma)^2/basic.t0,NREPS);    % compute deltas for realized Bayes means
                else                                    % otherwise usu the usual random samples from posterior distribution
                    bayesmeandelta = normrnd(0,basic.sigma/sqrt(basic.t0),NREPS,1);        % compute deltas for realized Bayes means
                end
            end
            [fignum, simOut] = DelaySimComputer(fignum, basic,advanced,mat,delta,delta+bayesmeandelta,normalnoise);  % do simulation analysis for bayesian predictive distribution
            if i==1
                tmp.simOut = repmat(simOut,length(advanced.simFreqDeltaVec),1);  % preallocate vector of simulation output structures for frequentist analysis
            else
                tmp.simOut(i) = simOut;  % preallocate vector of simulation output structures for frequentist analysis
            end
        end
        if advanced.keepAllOutput
            matout.simBayesOut = tmp.simOut;
        end
        [matout.outBayes] = DelaySimAnalysis( basic, advanced, tmp.simOut );

        % Do frequentist stuff next
        for i=1:length(advanced.simFreqDeltaVec)
            advanced.simTmpFlag = false;        % this is a temporary flag for use in sim computer, true for bayesian stopping bound, false for frequentist
            delta = advanced.simFreqDeltaVec(i);
            advanced.PLOTSIMS = ( (abs(abs(delta) - basic.sigma/sqrt(basic.t0)) < MYEPS) | (abs(delta) < MYEPS) ) & origPLOTSIMS;

            if ~advanced.CRN    % if CRN across bayes and frequent values of delta desired, keep same noise, otherwise get independently sampled noise
                normalnoise = normrnd(0,basic.sigma,basic.TMax,NREPS);  % compute the sampling noise here, to implement CRN across power tests and Bayesian means
            end
            meanvec = delta*ones(NREPS,1);        % for frequentist stats if the mean is known to be delta
            [fignum, tmp.simOut(i)] = DelaySimComputer(fignum, basic,advanced,mat,delta,meanvec,normalnoise);
        end
        if advanced.keepAllOutput
            matout.simFreqOut = tmp.simOut;
        end
        [matout.outFreq] = DelaySimAnalysis( basic, advanced, tmp.simOut );
%    end
end

figout = fignum;

end

