function [ matout ] = DelayOptimalBayesOneStage ( basic, advanced, tvectotest, mat )
%DelayOptimalBayesOneStage tries to find the optimal one-stage experiment
%for a Bayesian one-shot trial, among times in the vector tvectotest, one
%for each value of mu in muvec
% INPUTS
%   basic, advanced: created by DelayInputConstructor method, as in one of
%       the Set*.m files, to describe basic parameters of the sequential
%       trial, and for computational and output purposes.
%   tvectotest: row vector or column vector of values to test for number of
%      samples in one-shot trial, such that tvectotest=0:10:(basic.Tmax*10), for example,
%      tests with 0 samples, 10 samples, etc. up to 10*(the max number of
%      patients to be sampled in the sequential trial)
%   mat: assumes that DelayCurvesRecur (stage II) has already been called
%      with basic, advanced, and that the output is mat. It is also ok of
%       other routines have been called after that (such as DelayStageOne or
%       DelaySim...)
% OUTPUTS:
%   matout: a copy of the mat input structure, together with some new fields.
%       matout.OptOneShotLength: the optimal number of samples to take in a
%          one-shot sampling test, for each value of mu, among potential test durations in the vector tvectotest 
%       matout.OptOneShotReward: expected reward for that optimal duraction
%           one-shot sampling test (vector same length as muvec)
%       matout.OneShotTimesTested: keeps a column vector version of
%          tvectotest
% Designed for Delay Sequential Trials paper of Chick, Forster, Pertile (alpha order).
%
% (c) 2014, S Chick
% Created: 8 Oct 2014
% Last touched: 8 Oct 2014
% 
muvec = mat.muvec'; %
sigma = basic.sigma;
t0 = basic.t0;
theta = basic.theta;
tau = basic.tau;
TMAX = basic.TMax;
%PPatients = basic.PPatients;
ICost = basic.ICost;
c = basic.c;
advanced2 = advanced; advanced2.DoRegretIntegral = true;    % force integral for regret calculations for one stage procedure

% compute the value of having perfect information at time 0 and making the
% best decision
numpatients = basic.PPatients + (1-advanced.fixedP)*TMAX;
predvarsample = numpatients^2 * sigma^2 / t0;

if advanced.UnkVariance  % don't define enddof if variance of sampling distribution is known
    if advanced.UnkVarianceShape == -1
        dof = t0;         % use if there is a noninformative prior for mu, sigma^2 of sampling distribution
    else
        dof = 2*advanced.UnkVarianceShape;    % use if the prior is otherwise specified
    end
end

if advanced.UnkVariance  % don't define enddof if variance of sampling distribution is known
    RewardPIT0I = TerminalRewardFunctionUnk(muvec*numpatients-ICost,1.0,predvarsample,false,dof);
else
    RewardPIT0I = TerminalRewardFunction(muvec*numpatients-ICost,1.0,predvarsample,false);
end

NUMCHECKS = length(tvectotest);% number of checks for values of s in interval (0, tau]
svec= t0 + tvectotest(:)';      % create a row vector which will have posterior mean after all samples seen


PPatientvec = basic.PPatients*ones(size(svec)) + (1-advanced.fixedP)*(TMAX - svec);  % number of patients to be treated following adoption
if min(PPatientvec) <= 0
    'warning: fixedP is positive and caused DelayOptimalOneStage to consider negative number of patients: check PPatients and maximum one-stage samples size'
end
discountvector = theta.^(svec - t0 + tau);  % decision will be made tau units of time after last of s observations made
discountvector(svec==t0) = 1.0;  % except if decision is made immediately, in which case no discounting occures
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
    G0initcost = ones(length(muvec),1) * (-c*(1-theta.^(svec-t0))/(1-theta));
    if basic.online
        G0initcost = G0initcost + muvec(:) * (1-theta.^(svec-t0))/(1-theta);
    end
else
    G0initcost = - c * ones(length(muvec),1) * (svec-t0);
    if basic.online
        G0initcost = G0initcost + muvec(:) * (svec-t0);
    end
end
G0 = G0initcost + G0base;
[beststage1, bests] = max(G0,[],2);
bestsvec = svec(bests)-t0;      % get optimal number of samples

if advanced.DOPLOT
    figure(9); hold off;
    plot(muvec,beststage1','-.',muvec,max(muvec*PPatientvec(1)-ICost,0),'--',muvec,RewardPIT0I,'-');
    legend('optimal E[one-shot reward]','E[reward] of no samples','EVPI');
    figure(9999); hold off;
    plot(muvec(:),bestsvec(:),'-');
    legend('num samples of optimal one-shot trial (bayes)');
end

% check to see if taking ZERO samples is best
beststage1 = max(beststage1,max(muvec'*PPMax-ICost,0));      % check if s=0 is even better than positive s in (0, tau]
bestsvec(beststage1 == max(muvec'*PPMax-ICost,0)) = 0;

% set output variables, including adding fields to the output structure.
matout = mat;
matout.OptOneShotReward = beststage1(:);
matout.OptOneShotLength = bestsvec(:);
matout.OneShotTimesTested = tvectotest(:);

end

