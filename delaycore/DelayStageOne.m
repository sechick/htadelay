function [ matout ] = DelayStageOne ( basic, advanced, mat )
%DelayStageOne implements the stage 1 of the backward recursion for the
%delay clinical trial paper with Chick, Forster, Pertile
% INPUTS
%   basic, advanced: created by DelayInputConstructor method, as in one of
%       the Set*.m files, to describe basic parameters of the sequential
%       trial, and for computational and output purposes. Once created by
%       the constructor, values can be tweaked, then should be run through
%       DelayInputValidator to insure all data is ok, and to compute
%       intermediate parameters which are required by this function.
%   mat: a structure of output values of multiple sorts, with following fields:
%           mat.B0mat;          % value function
%           mat.Puppermat;  % Probability of stopping due to hitting/exceeding upper boundary
%           mat.Pnewmat;      % Probability that if one stopped at a given (t, mu), that new technology is adopted
%           mat.tmat;            % a row vector of time values for columns of the above matrices
%           mat.muvec;       % a column vector of mu values for rows of above matrices
%           mat.tvec;        % tvec is time values for the boundaries: may be different from tmat, as tvec tries to compress
%           mat.bndupper;    % info requirements to get at curves of boundary, whereas tmat is equally spaced, in order
%           mat.bndlower;    % to help contour plot, bndupper and bndlower are upper and lower stopping boundaries
%           mat.B0vec;       % first column of B0mat (might trim this as redundant), expected reward to go, at muvec, at start of stage II assuming tau samples pending
%
%           mat.RewardPITmaxIII; % Expected reward with perfect information at time t=Tmax
%           mat.RewardTmaxIII; % Expected reward, assuming tau samples pending at time t=Tmax
%           mat.RewardPIMatII; % Expected reward with perfect information at time t in stage II
%           mat.RewardMatII; % Expected reward, assuming tau samples pending at time t in stage II
%           mat.RewardPITmaxI; % Expected reward with perfect information at time t=0 in stage I
%           note that B0hat is computed below to give the expected reward
%           of optimally moving forward.
%               Note: Regret will be RewardPI - Reward for each of these three stages
%
% OUTPUTS:
%   mat: a copy of the mat input structure, together with some new fields.
%       mat.B0hat: vector of optimal rewards, same length as muvec
%       mat.bestsvec: vector of number of samples to take in stage 1. Will be 0 if
%           it is optimal to select an alterative directly without sampling, will
%           be tau+1 if it is optimal to take tau samples and move to stage 2, and
%           will be between 0 and tau if it is optimal to take a few samples in
%           stage 1, and make a decision when the info from those samples comes in,
%           (as opposed to not sampling at all, or not sampling enough to enter
%           stage 1). If it is optimal to sample for at least one sample in stage
%           II, then bestsvec will have an entry of tau+1, do differentiate it from
%           the alternative of sampling exactly tau but no more.
%       mat.Threshpoint; % indices into muvec for points A, B, C, D
%           delineate the areas where one does not sample, where one
%           samples less than tau, and where one continues to stage II
%           (when these points are defined, otherwise indices are set to -1)
%
% Designed for Delay Sequential Trials paper of Chick, Forster, Pertile (alpha order).
%
% (c) 2014, S Chick
% Created: 14 April 2014
% Updated: 11 Feb 2015 to include threshpoints. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 1: Initialize variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%MYEPS = 10e-7;
muvec = mat.muvec'; %
sigma = basic.sigma;
t0 = basic.t0;
theta = basic.theta;
tau = basic.tau;
TMAX = basic.TMax;
%PPatients = basic.PPatients;
ICost = basic.ICost;
c = basic.c;
B0vec = mat.B0vec;
advanced2 = advanced; advanced2.DoRegretIntegral = true;    % force integral for regret calculations for one stage procedure

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 2: compute the value of having perfect information at time 0 and making the
% best decision
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numpatients = basic.PPatients + (1-advanced.fixedP)*TMAX;
predvarsample = numpatients^2 * sigma^2 / t0;

% NOTE FOR SPECIAL CASE OF UKNOWN VARIANCE:
% Using UnkVariance to modify pu and pd (probability in trinomial
% tree) helps with getting right variance in unknown mean below, but the
% code below does not account for the learning about the unknown variance
% with samples. When running code with pu and pd influencing mean and not
% variance, the upper and lower boundaries have odd behaviors (they curve
% back in toward I/P, and even cross each other in some cases), and the VOI
% is underestimated - implying that Stage I sampling, when correctly
% computed, has more value relative to the estimate from Stage II sampling.
% This incoherency means that it is probably better to plug in an estimator
% of the variance for the upper and lower boundaries when handling the case
% of Unknown Variance. This is best avoided in BOTH Phase II (DelayCurvesRecur.m) and
% Phase I (DelayStageOne.m) routines.
% Thus, we issue a warning message here to indicate that even though code
% is here for pu and pd to handle the mean, we do not use it in
% computations. It is left here in order to enable further development on
% the calculations for unknown variance at a later date.

%if advanced.UnkVariance
%    warning('DelayCurvesRecur: Unknown variance not implemented in Phase I analysis, fudge factor applied later in monte carlo simulations');
%    advanced.UnkVariance = false;
%end

if advanced.UnkVariance  % don't define enddof if variance of sampling distribution is known
    if advanced.UnkVarianceShape == -1
        dof = t0;         % use if there is a noninformative prior for mu, sigma^2 of sampling distribution
    else
        dof = 2*advanced.UnkVarianceShape;    % use if the prior is otherwise specified
    end
    RewardPIT0I = TerminalRewardFunctionUnk(muvec*numpatients-ICost,1.0,predvarsample,false,dof);
else
    RewardPIT0I = TerminalRewardFunction(muvec*numpatients-ICost,1.0,predvarsample,false);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 3: For cases of delay or no delay in samples, check for best
% one-stage allocation.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if tau > 0              % IF tau > 0, then we have the opprtunity to have a delay
    if advanced.StageOneChecks <= 0 % SEC: Updated 19 nov 2014 to allow for user to specify number of checks for optimal stage 1 sampling size
        deltas = 1;                 % Give 0 or a negative number to check for s=1, 2, 3, ..., \tau. 
                                    % Give positive integer to specify the number of equally spaced values of time to check
    else
        deltas = tau/advanced.StageOneChecks;   % check at times t0 + [ tau, 2 tau, ... tau*tau] / tau
    end
    svec=(t0+deltas):deltas:t0+tau; % if sampling to occur in stage 1, must take at least 1 sample
    svec(length(svec))=t0+tau;      % insure that the last entry is exactly t0+tau

    PPatientvec = basic.PPatients*ones(size(svec)) + (1-advanced.fixedP)*(TMAX - svec);  % number of patients to be treated following adoption
    discountvector = theta.^(svec - t0 + tau);  % decision will be made tau units of time after last of s observations made
    predvarvec = (PPatientvec.^2 * sigma^2) .* (svec - t0) ./ (t0*svec);   % predictive variance of posterior mean to be seen after svec samples
    PPMax = basic.PPatients + (1-advanced.fixedP)*TMAX;       % total number of patients in contract plus max in trial who could be switched

    G0base = zeros(length(muvec),length(svec));     % preallocate size of matrix for G0
    for i=1:length(svec)            % for each time value s for sampling, compute terminal expected reward for each
        % This 'RegretPenalty' code, if changed, needs to be changed in several
        % places: one place in DelayStageOne, two places in DelayCurvesRecur
        postvar = (PPatientvec(i)^2 * sigma^2) /  (svec(i));      % posterior variance after all samples arrive which can be seen by time of decision
        if advanced.UnkVariance
            G0base(:,i) = TerminalRewardFunctionUnk(muvec*PPatientvec(i)-ICost,discountvector(i),predvarvec(i),false,dof);
            RewardI = TerminalRegretUnk(muvec*PPatientvec(i) - ICost, discountvector(i)^(1- advanced.nochangeallowed), predvarvec(i), postvar, advanced2, dof );  % compute expected regret of a decision, assuming one may have option to continue sampling. This would be the expected value of perfect information, given the information state at the time of stopping, assuming one can continue
        else
            G0base(:,i) = TerminalRewardFunction(muvec*PPatientvec(i)-ICost,discountvector(i),predvarvec(i),false);
            RewardI = TerminalRegret(muvec*PPatientvec(i) - ICost, discountvector(i)^(1- advanced.nochangeallowed), predvarvec(i), postvar, advanced2 );  % compute expected regret of a decision, assuming one may have option to continue sampling. This would be the expected value of perfect information, given the information state at the time of stopping, assuming one can continue
        end
        Regretcontin = RewardPIT0I - RewardI;
        G0base(:,i) = G0base(:,i) - abs(advanced.RegretPenalty) * Regretcontin(:);        % set initial estimate of reward to go at terminal time tHoriz
    end

    if theta<1
        B0hat= -c *(1-theta^tau)/(1-theta) + theta^tau * B0vec;
        G0initcost = ones(length(muvec),1) * (-c*(1-theta.^(svec-t0))/(1-theta));
        if basic.online
            B0hat = B0hat + muvec*(1-theta^tau)/(1-theta);
            G0initcost = G0initcost + muvec(:) * (1-theta.^(svec-t0))/(1-theta);
        end
    else
        B0hat= - c * tau + B0vec;
        G0initcost = - c * ones(length(muvec),1) * (svec-t0);
        if basic.online
            B0hat = B0hat + muvec * tau;
            G0initcost = G0initcost + muvec(:) * (svec-t0);
        end
    end
    G0 = G0initcost + G0base;
    [beststage1, bests] = max(G0,[],2);
    bestsvec = min(svec(bests)-t0,tau);      % get optimal number of samples

    if advanced.DOPLOT
        figure(9); hold off;
        plot(muvec,B0hat,'-',muvec,beststage1','-.',muvec,max(muvec*PPatientvec(1)-ICost,0),'--');
        legend('value of continuing to phase II','optimal stage I value','stop now');
    end

    % check if taking 0 samples is optimal among policies which do not
    % advance to Stage II.
    if advanced.nochangeallowed % In stage I, no change allowed implies that stopping now before new data arrives
                % requires that one pick the best. One can only change if one
                % continues to stage II.
        beststage1 = max(muvec'*PPMax-ICost,0);      % check if s=0 is even better than positive s in (0, tau]
        bestsvec(beststage1 == max(muvec'*PPMax-ICost,0)) = 0;
    else  % The default is advanced.nochangeallowed=false, the following code, which permits stopping
                % before stage II after a few samples have been taken, and
                % allows one to decide which alternative is best.
        beststage1 = max(beststage1,max(muvec'*PPMax-ICost,0));      % check if s=0 is even better than positive s in (0, tau]
        bestsvec(beststage1 == max(muvec'*PPMax-ICost,0)) = 0;
    end

    % if optimal to continue into stage II sampling, then set optimal number of
    % samples to tau + 1. 
    % Note: formally this is correct when the variance is known. When the
    % variance is unknown we use a kludge: we test to see if the VOI from
    % one stage is maximizes on [0, tau] at tau and consider that to be a
    % decision to go to stage II (the reason is that in the variance
    % unknown case, the VOI in stage I is calculated with student
    % distributions, where as for stage II diffusion we assume known
    % variance (variance is known for diffusion sampled on a continuous
    % interval). Thus, the B0hat is not directly commensurate with the
    % stage I VOI due to the known/unknown variance mismatch.
    if advanced.UnkVariance
        bestsvec(bestsvec(:) >= tau) = tau + 1;    % this is the kludge: go to stage II if it appears VOI is increasing at tau samples
    else
        bestsvec(beststage1(:) < B0hat(:)) = tau + 1;    
    end
    %figure(plot(muvec(:),beststage1(:)-B0hat(:))
    %tmpvec=beststage1(:)-B0hat(:)

    B0hat = max(B0hat(:),beststage1(:));  %again this is right for known variance, an approx (lower bound) for unknown variance.

    if advanced.verbose
        bestsrealzero = (bestsvec==0);
        bestszero = (bestsvec<=tau) & (bestsvec>0);
        beststau = (bestsvec > tau);
        sum(bestsrealzero);
        sum(bestszero);
        sum(beststau);
    end

    % set output variables, including adding fields to the output structure.
    matout = mat;
    matout.B0hat = B0hat;
    matout.bestsvec = bestsvec;
    matout.optENumSamps = bestsvec;
    matout.optENumSamps(bestsvec == (tau + 1)) = tau + mat.ENumSamps(bestsvec == (tau + 1));      % compute mean optimal number of samples over both stages
else        % tau = 0: this means the data from stage II is what matters most - there is no stage I
    matout = mat;
    matout.B0hat = B0vec;           % no cost in stage 1
    bestsvec =  zeros(size(mat.ENumSamps));
    matout.bestsvec = bestsvec;
    matout.bestsvec(mat.ENumSamps>0) = tau + 1;   
    matout.optENumSamps = mat.ENumSamps;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 4: Now compute statistics for one stage trial of all TMax patients
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PPMax = basic.PPatients;
discountfactor = theta^(TMAX+(1-advanced.nochangeallowed)*tau);
predvarsample = (PPMax^2 * sigma^2 * (TMAX - advanced.nochangeallowed*tau)) / ((t0 + TMAX - advanced.nochangeallowed * tau) * t0);
postvar = (PPMax^2 * sigma^2) /  (t0 + TMAX - advanced.nochangeallowed*tau);      % posterior variance after all samples arrive which can be seen by time of decision
if advanced.UnkVariance
    if advanced.UnkVarianceShape == -1
        dof = t0 + TMAX - advanced.nochangeallowed*tau;         % use if there is a noninformative prior for mu, sigma^2 of sampling distribution
    else
        dof = 2*advanced.UnkVarianceShape + TMAX - advanced.nochangeallowed*tau;    % use if the prior is otherwise specified
    end
    [OneShotReward, OneshotPCS] = TerminalRegretUnk(muvec*PPMax - ICost, discountfactor^(1- advanced.nochangeallowed), predvarsample, postvar, advanced2, dof );  % compute expected regret of a decision, assuming one may have option to continue sampling. This would be the expected value of perfect information, given the information state at the time of stopping, assuming one can continue
else
    [OneShotReward, OneshotPCS] = TerminalRegret(muvec*PPMax - ICost, discountfactor^(1- advanced.nochangeallowed), predvarsample, postvar, advanced2 );  % compute expected regret of a decision, assuming one may have option to continue sampling. This would be the expected value of perfect information, given the information state at the time of stopping, assuming one can continue
end
OneshotRegret = RewardPIT0I - OneShotReward;
BayesOneshotExpectedReward = OneShotReward(:) - abs(advanced.RegretPenalty) * OneshotRegret(:);        % set initial estimate of reward to go at terminal time tHoriz

% need to take out sampling costs!
if basic.theta < 1                    % subtract off the sampling costs, depending on whether there is discounting or not
    rewarddelta2 = (basic.online * muvec - basic.c) * (1-basic.theta^basic.TMax)/(1-basic.theta);
else
    rewarddelta2 = (basic.online * muvec - basic.c) * basic.TMax;
end
BayesOneshotExpectedReward = BayesOneshotExpectedReward' + rewarddelta2;

matout.OneshotRegret = OneshotRegret;
matout.OneshotPCS = OneshotPCS;
matout.BayesOneshotExpectedReward = BayesOneshotExpectedReward;
zfreq = abs(muvec*PPMax - ICost)/(PPMax*sigma/sqrt(basic.TMax));       % SEC: fixed one shot pcs and eoc, 30 dec 2014.
if advanced.UnkVariance
    if advanced.UnkVarianceShape == -1
        dof = t0;
    else
        dof = 2*advanced.UnkVarianceShape;
    end
    matout.FreqOneShotPCS = tcdf(zfreq,dof);
    matout.FreqOneShotEOC = PPMax*(sigma/sqrt(basic.TMax))* PsiNormUV(zfreq,dof);
else
    matout.FreqOneShotPCS = normcdf(zfreq);
    matout.FreqOneShotEOC = PPMax*(sigma/sqrt(basic.TMax)) * PsiNorm(zfreq);
end

matout.Threshpoint = UtilGetABCDVector(bestsvec, tau);
if advanced.verbose
	display( 'bestsvec' ) ; 
	display( bestsvec ) ; 
	display('matout.Threshpoint') ; 
	display( matout.Threshpoint ) ; 
end

matout.RewardPIT0I = RewardPIT0I;
tvectotest = 0:1:basic.TMax;
[ matout ] = DelayOptimalBayesOneStage ( basic, advanced, tvectotest, matout );

tstval = (matout.bestsvec(1) == 0) && (matout.bestsvec(end) == 0);
if ~tstval
    warning('DelayStageOne: range of nonzero optimal allocations may exceed [mumin, mumax], consider making basic.mumax larger and/or basic.mumin smaller');
%    muminmax = [basic.mumin basic.mumax];
%    boundrange = [min(mat.bndlower) max(mat.bndupper)];
%    warning('%s, ',muminmax);
%    warning('%s, ',boundrange);
end

end

