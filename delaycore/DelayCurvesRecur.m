function [retval, mat] = DelayCurvesRecur( basic, advanced )
% For delay health tech assessment optimal stopping project with Martin, Paolo, Steve.  
% This implements the stage II recursion of the delay sequential trials
% paper.
%
% Code is 'as is' and no guarantees for correctness.  Search for 'FIX'
%
% INPUTS
%   basic, advanced: created by DelayInputConstructor method, as in one of
%       the Set*.m files, to describe basic parameters of the sequential
%       trial, and for computational and output purposes. Once created by
%       the constructor, values can be tweaked, then should be run through
%       DelayInputValidator to insure all data is ok, and to compute
%       intermediate parameters which are required by this function.
% OUTPUTS
%   retval: true if upper and lower boundaries are within range of muvec,
%           and false if not (this signals that mumax and mumin should be wider) 
%   mat: a structure of output values of multiple sorts, with following fields:
%           mat.B0mat;          % value function
%           mat.Puppermat;  % Probability of stopping due to hitting/exceeding upper boundary
%           mat.Pnewmat;      % Probability that if one stopped at a given (t, mu), that new technology is adopted
%           mat.PCSmat;      % Probability of correctly selecting the best alternative.
%           mat.tmat;            % a row vector of time values for columns of the above matrices
%           mat.muvec;       % a column vector of mu values for rows of above matrices
%           mat.bndupper;    % info requirements to get at curves of boundary, whereas tmat is equally spaced, in order
%           mat.bndlower;    % to help contour plot, bndupper and bndlower are upper and lower stopping boundaries
%           mat.B0vec;       % first column of B0mat (might trim this as redundant), expected reward to go, at muvec, at start of stage II assuming tau samples pending
%
%           mat.Bayespcsvec = pcsin;    % PCS at time stage II starts, given tau samples to arrive and stopping boundary used going forward
%           mat.RewardPITmaxIII; % Expected reward with perfect information at time t=Tmax
%           mat.RewardTmaxIII; % Expected reward, assuming tau samples pending at time t=Tmax
%           mat.RewardPIMatII; % Expected reward with perfect information at time t in stage II
%           mat.RewardMatII; % Expected reward, assuming tau samples pending at time t in stage II
%           mat.RewardPITmaxI; % Expected reward with perfect information at time t=0 in stage I
%           mat.RewardTmaxI; % Expected reward, assuming tau samples pending at time t=0 in stage I
%               Note: Regret will be RewardPI - Reward for each of these three stages
%
% DEPENDENCES:
%  CALLS:
%   DelayInputValidator(): validates input parameters (basic and advanced),
%       and long the way resets or initializes some paramerters as needed. 
%   TerminalRewardFunction(): terminal reward function.
%  CALLED BY:
%   DelayDriver: code that illustrates calling this routime,
%          getting Gittins index, etc.
%  OUTPUT USED BY:
%   mat: is used by DelayStageOne(), the stage I part of the dynamic
%       program.
%
% Designed for Delay Sequential Trials paper of Chick, Forster, Pertile
% (alpha order).
%
% (c) 2014, S Chick, M Forster, P Pertile
% Created: 14 April 2014
% Last touched: 12 June 2014
% 

%======= BLOCK 1 ========
% Validate the input parameters, resetting some values if needed, and put
% the values into the local space as needed. Those adjusted values can be
% returned to the original calling routine

[ basic, advanced, rval, messages ] = DelayInputValidator( basic, advanced );
if ~rval || ~isempty(messages)
    warning(messages);
end

% take values from parameter structures and set local variables for use.
online = basic.online;
sigma = basic.sigma;
tau = basic.tau;
%mu0=basic.mu0;
t0 = basic.t0;
c = basic.c;
P = basic.PPatients;
I = basic.ICost;
theta = basic.theta;
%TMax = basic.TMax;
[tvec, muvec] = UtilMakeTvecMuvec(basic, advanced); % tvec should run from t0+tau to approximately t0+Tmax in small increments
tlen = length(tvec);
dt = advanced.dt;
dmu = advanced.dmu;
DOPDE = advanced.DOPDE;

if theta==1     % put in discounted costs for added samples
    thetadtfactor = dt;
else
    thetadtfactor = -(1-theta^dt)/log(theta);
end

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
if advanced.UnkVariance
    warning('DelayCurvesRecur: PDE approach with unknown variance may underestimate value of stage II sampling');
end

%======= BLOCK 2 ========
% allocate parameters, vectors et for local computation
musize = length(muvec);
flipmu = fliplr(muvec);

rawboundtop=zeros(size(tvec));     % use this to store the upper 'boundary' of the stopping region, from grid
rawboundbot=zeros(size(tvec));     % use this to store the lower 'boundary' of the stopping region, from grid

dtdmu2ratio = dt / dmu^2;          % ratio which is used frequently

% INITIALIZE TERMINAL CONDITIONS
discountfactor = (theta^tau);   % amount of discounting from a given time until the data for that time arrives
i=tlen;
curt=tvec(i);
numpeopletreated = P;
predvarsample = (numpeopletreated^2 * sigma^2 * tau) / ((curt-tau) * curt);    % note: at time curt, the effective number of samples is curt-tau, because tvec starts at time t0+tau, and we have tau samples not yet arrived

if advanced.UnkVariance  % don't define enddof if variance of sampling distribution is known
    if advanced.UnkVarianceShape == -1
        enddof = curt - tau;         % use if there is a noninformative prior for mu, sigma^2 of sampling distribution
    else
        enddof = curt + 2*advanced.UnkVarianceShape - tau - t0;    % use if the prior is otherwise specified
    end
end

% This 'RegretPenalty' code, if changed, needs to be changed in several
% places: one place in DelayStageOne, two places in DelayCurvesRecur
if advanced.UnkVariance  % terminal regret and expected reward with perfect info are different if sampling variance is unknown
    Cinit = TerminalRewardFunctionUnk(muvec*numpeopletreated-I,discountfactor,predvarsample,advanced.nochangeallowed, enddof);
    postvar = (numpeopletreated^2 * sigma^2) /  (enddof + tau * (1- advanced.nochangeallowed) ); % posterior variance after all samples arrive which can be seen by time of decision
    [RewardTmaxIII, pcsstopnow] = TerminalRegretUnk(muvec*numpeopletreated - I, discountfactor^(1- advanced.nochangeallowed), ...
            predvarsample, postvar, advanced, enddof );  % compute expected regret of a decision, assuming one may have option to continue sampling. This would be the expected value of perfect information, given the information state at the time of stopping, assuming one can continue
    % Get expected reward with perfect information at time t=Tmax, assuming no discounting
    RewardPITmaxIII = TerminalRewardFunctionUnk(muvec*numpeopletreated-I,1.0,(numpeopletreated^2 * sigma^2) / (curt-tau),false,enddof);
else
    Cinit = TerminalRewardFunction(muvec*numpeopletreated-I,discountfactor,predvarsample,advanced.nochangeallowed);
    postvar = (numpeopletreated^2 * sigma^2) /  (curt + tau * (1- advanced.nochangeallowed) ); % posterior variance after all samples arrive which can be seen by time of decision
    [RewardTmaxIII, pcsstopnow] = TerminalRegret(muvec*numpeopletreated - I, discountfactor^(1- advanced.nochangeallowed), ...
            predvarsample, postvar, advanced );  % compute expected regret of a decision, assuming one may have option to continue sampling. This would be the expected value of perfect information, given the information state at the time of stopping, assuming one can continue
    % Get expected reward with perfect information at time t=Tmax, assuming no discounting
    RewardPITmaxIII = TerminalRewardFunction(muvec*numpeopletreated-I,1.0,(numpeopletreated^2 * sigma^2) / (curt-tau),false);
end
pcsstopprev = pcsstopnow; 

%           mat.RewardPITmaxIII; % Expected reward with perfect information at time t=Tmax
%           mat.RewardTmaxIII; % Expected reward, assuming tau samples pending at time t=Tmax
RegretIII = RewardPITmaxIII - RewardTmaxIII;  % compute expected regret of a decision, assuming stopping. This would be the expected value of perfect information, given the information state at the time of stopping, discounted to tau time steps later if needed, for commensurability with Cinit.
Cinit = Cinit - abs(advanced.RegretPenalty) * RegretIII;        % set initial estimate of reward to go at terminal time tHoriz

%stoprewardvec = Cinit;
Cin = Cinit;        % set initial estimate of reward to go at terminal time tHoriz
Cintmp = Cinit;     % done just to size Cintmp correctly - used in intermediate calculations

if advanced.DOPLOT
    figure(1); hold off;
    plot(muvec,Cinit,'-');
    legend('Terminal reward')
    xlabel('posterior mean');ylabel('Cin');
end

%Issue: Difference between accept
%option to decision in tau time steps or not to allow waiting for tau time
%steps / versus the stopping boundary issue which says wait for tau time
%steps before selecting the alternative
%[tstval, tstindx] =  max( Cinit > 0);        % find first element with positive value for expected reward at end of horizon
%if tstval > 0 
%    rawboundbot(i) = muvec(tstindx);     % use this to store the lower 'boundary' of the stopping region, from grid
%    rawboundtop(i) = muvec(tstindx);     % use this to store the upper 'boundary' of the stopping region, from grid
%else
    rawboundbot(i) = I/numpeopletreated;     % use this to store the lower 'boundary' of the stopping region, from grid
    rawboundtop(i) = I/numpeopletreated;     % use this to store the upper 'boundary' of the stopping region, from grid
%end

% the following three are other interesting values of interest to a
% Bayesian. These computations are valid everywhere in the stopping region,
% but they will require tweaks below to account for their dynamics in the
% continuation region. Keep these all as column vectors. The + before a
% logical expression converts the 0s and 1s to double precision rather than
% logical values.
Pupperin = +(muvec >= rawboundtop(i));   % if on or above upper boundary, prob of stopping on/above top is 1, otherwise it is 0, at end of time horizon
Pupperin(muvec==rawboundtop(i)) = 0.5;  % call it 50/50 at equality
if advanced.UnkVariance % handle case of unknown variance and known variance separately
    Ppicknewin = TerminalProbPickNewUnk(muvec*numpeopletreated - I, predvarsample,advanced.nochangeallowed, enddof); % find probability of picking new alternative, given 
else
    Ppicknewin = TerminalProbPickNew(muvec*numpeopletreated - I, predvarsample,advanced.nochangeallowed); % find probability of picking new alternative, given
end
ENumSampsin = 0*muvec;                    % initialize the expected number of samples until stopping to 0: at end of horizon no more samples can be taken

% start to build up the matrices with return values. Time values in
% columns, different mu values in rows
B0mat = Cinit(:);
Puppermat = Pupperin(:);
Pnewmat = Ppicknewin(:);
PCSmat = pcsstopnow(:);
maxindx = 1; minindx = length(Cintmp); nonmonotoniclower = 0; nonmonotonicupper = 0;

RewardPIMatII = RewardPITmaxIII(:); % Expected reward with perfect information at time t=Tmax
RewardMatII = RewardTmaxIII(:); % Expected reward, assuming tau samples pending at time t=Tmax
tmat = tvec(i);

for i=(tlen-1):-1:1
    curt = tvec(i);     % set current time index (starting from t0+tau up to t0+TMAX
    numpeopletreated = P + (1-advanced.fixedP) * (tvec(tlen) - curt);
    predvarsample = (numpeopletreated^2 * sigma^2 * tau) / ((curt-tau) * curt);    % note: at time curt, the effective number of samples is curt-t0, because tvec starts at time t0, and we have tau samples not yet arrived
    if advanced.UnkVariance  % don't define enddof if variance of sampling distribution is known
        if advanced.UnkVarianceShape == -1
            enddof = curt - tau;         % use if there is a noninformative prior for mu, sigma^2 of sampling distribution
        else
            enddof = curt + 2*advanced.UnkVarianceShape - tau - t0;    % use if the prior is otherwise specified
        end
        fudgefactor = sqrt( enddof / (enddof - 2) );
        postvar = (numpeopletreated^2 * sigma^2) /  (enddof + tau * (1- advanced.nochangeallowed) );  % posterior variance after all samples arrive which can be seen by time of decision
    else
        fudgefactor = 1.0;
        postvar = (numpeopletreated^2 * sigma^2) /  (curt + tau * (1- advanced.nochangeallowed) );  % posterior variance after all samples arrive which can be seen by time of decision
    end
    
    % This 'RegretPenalty' code, if changed, needs to be changed in several
    % places: one place in DelayStageOne, two places in DelayCurvesRecur
    if advanced.UnkVariance  % don't define enddof if variance of sampling distribution is known
        stoprewardvec = TerminalRewardFunctionUnk(muvec*numpeopletreated-I,discountfactor,predvarsample,advanced.nochangeallowed,enddof);
        [Rewardmax, pcsstopnow] = TerminalRegretUnk(muvec*numpeopletreated - I, discountfactor^(1- advanced.nochangeallowed), ...
                predvarsample, postvar, advanced, enddof );  % compute expected regret of a decision, assuming one may have option to continue sampling. This would be the expected value of perfect information, given the information state at the time of stopping, assuming one can continue
        if abs(advanced.RegretPenalty) ~= 0
            RewardPI = TerminalRewardFunctionUnk(muvec*numpeopletreated-I,1.0,(numpeopletreated^2 * sigma^2) / (curt-tau),false,enddof);
            Regretcontout = RewardPI - Rewardmax; 
            stoprewardvec = stoprewardvec - abs(advanced.RegretPenalty) * Regretcontout;        % set initial estimate of reward to go at terminal time tHoriz
        end
    else
        stoprewardvec = TerminalRewardFunction(muvec*numpeopletreated-I,discountfactor,predvarsample,advanced.nochangeallowed);
        [Rewardmax, pcsstopnow] = TerminalRegret(muvec*numpeopletreated - I, discountfactor^(1- advanced.nochangeallowed), ...
                predvarsample, postvar, advanced );  % compute expected regret of a decision, assuming one may have option to continue sampling. This would be the expected value of perfect information, given the information state at the time of stopping, assuming one can continue
        if abs(advanced.RegretPenalty) ~= 0
            RewardPI = TerminalRewardFunction(muvec*numpeopletreated-I,1.0,(numpeopletreated^2 * sigma^2) / (curt-tau),false);
            Regretcontout = RewardPI - Rewardmax; 
            stoprewardvec = stoprewardvec - abs(advanced.RegretPenalty) * Regretcontout;        % set initial estimate of reward to go at terminal time tHoriz
        end
    end
    
    pu = fudgefactor * dtdmu2ratio * sigma^2 / ((curt - tau) *((curt-tau) + dt)) / 2;   % probability of going up depends on t now
    pd = pu;                    % probability of going down
    psame = 1 - pu - pd ;        % probability of going straight (less a factor from the discounting
    
    if DOPDE  % Use free boundary iteration to compute Cintmp, the vector containing value of continuing
        % at present, this code is debugged for the case of known variance.
        % for the case of unknown variance, the results are not computing
        % well. THus, we have set DOPDE false above so that the KG* naive
        % one-step look-ahead approximation is used.
        if psame < 0.001
            warning('DelayCurvesRecur: oops, psame smaller than expected');
        end
        Cintmp(2:(musize-1)) = pu*Cin(3:musize) + psame*Cin(2:(musize-1)) + pd*Cin(1:(musize-2));
        Cintmp = theta^dt * Cintmp + thetadtfactor * (muvec*online-c);       % and add in sampling cost and the online sampling reward if appropriate
    else % Use KG* type approach to compute Cintmp, the vector containing value of continuing
        numrepsleft = basic.t0 + basic.TMax - curt;
        % FIX: KGset hardcode below could probably be replaced with 2 or 3
        % paramters in the advanced structure, but ...
%        KGset = unique([0 2.^(0:min(2,max(1,floor(log2(numrepsleft))))/3) numrepsleft]); % check powers of 2 up to number of remaining samples, and all remaining samples
        KGset = unique([dt 2.^((0:.25:max(1.5,min(100,ceil(log2(numrepsleft))))))]); % check powers of 2 up to number of remaining samples, and all remaining samples
        for j=1:length(KGset)
            % see what happens if we take an extra KGset(j) samples before
            % stopping to observe the remaining samples
            predtst = (numpeopletreated^2 * sigma^2 * (tau + KGset(j))) / ((curt-tau) * (curt + KGset(j)));    % note: at time curt, the effective number of samples is curt-t0, because tvec starts at time t0, and we have tau samples not yet arrived
            if advanced.UnkVariance
            	posttst = (numpeopletreated^2 * sigma^2) /  (enddof + KGset(j) + tau * (1- advanced.nochangeallowed) );  % posterior variance after all samples arrive which can be seen by time of decision
                OneShotReward = TerminalRewardFunctionUnk(muvec*numpeopletreated - I, theta^KGset(j) * discountfactor^(1- advanced.nochangeallowed), predtst, advanced.nochangeallowed, enddof );  % compute expected regret of a decision, assuming one may have option to continue sampling. This would be the expected value of perfect information, given the information state at the time of stopping, assuming one can continue
            else
            	posttst = (numpeopletreated^2 * sigma^2) /  (KGset(j) + curt + tau * (1- advanced.nochangeallowed) );  % posterior variance after all samples arrive which can be seen by time of decision
                OneShotReward = TerminalRewardFunction(muvec*numpeopletreated - I, theta^KGset(j) * discountfactor^(1- advanced.nochangeallowed), predtst,advanced.nochangeallowed );  % compute expected regret of a decision, assuming one may have option to continue sampling. This would be the expected value of perfect information, given the information state at the time of stopping, assuming one can continue
            end
            if theta==1     % put in discounted costs for added samples
                OneShotReward = OneShotReward + KGset(j) * (muvec*online - c);
            else    % -(1-basic.theta^dt)/log(basic.theta) is integral of theta^t dt from 0 to KGset(j)
                OneShotReward = OneShotReward - ((1-theta^KGset(j))/log(theta)) * (muvec*online - c);
            end
            if j==1
                Cintmp = OneShotReward;
            else
                Cintmp = max(Cintmp,OneShotReward);
            end
        end
    end
    Cintmp(musize)=stoprewardvec(musize);   % now, the discounting and sampling cost might given incorrect values for the extremes of mu,
    Cintmp(1)=stoprewardvec(1);             % so we again must reset their values to indure that we remain in the stopping set

    %find stopping boundary
    Cout = max(Cintmp,stoprewardvec);          % implement maximizer for bellman equation
    stopped = ( Cintmp <= stoprewardvec );     % stopped(i) is true if it is optimal to stop when mean is muvec(i)
    oldmaxindx = maxindx; oldminindx = minindx;
    [minval, minindx] = min(stopped);         
    
    if (minval > 0)             % if there is not a nonempty stopping region ...
        [~, minindx] = min(muvec*numpeopletreated-I<0);
        maxindx = minindx; %length(muvec)+1-minindx; % Thanks Paolo! (SEC touched lines 203 and 205 on 7 oct 2014)
        rawboundbot(i) = muvec(minindx);%I/numpeopletreated;     % use this to store the lower 'boundary' of the stopping region, from grid
        rawboundtop(i) = muvec(maxindx);%flipmu(maxindx);%I/numpeopletreated;     % use this to store the upper 'boundary' of the stopping region, from grid
    else     % there is a nontrivial stopping region
        rawboundbot(i) = muvec(minindx);    % grid point in bottom boundary is in contin set
        flipstop = fliplr(stopped);
        [maxval, maxindx] = min(flipstop);   % this finds first point in the continuation set
        if maxval > 0           % if there not a value in the continuation set...
            rawboundtop(i) = muvec(length(muvec));  % ... assume stop bound is as high as possible (this condition should never happen, as there is at least one 1 in stopvec
        else
            rawboundtop(i) = flipmu(maxindx); % get first grid point in continuation set
        end
        maxindx = length(muvec)+1-maxindx;
    end
    % some test code to see if enforcing monotonicity of continuation set
    % results in an overcoming of some numerical stability issues
    if minindx > oldminindx
        nonmonotoniclower = nonmonotoniclower + 1;
        minindx = min(minindx,oldminindx);
        rawboundbot(i) = muvec(minindx);
    end
    if maxindx < oldmaxindx
        nonmonotonicupper = nonmonotonicupper + 1;
        maxindx = max(maxindx,oldmaxindx);
        rawboundtop(i) = muvec(maxindx);
    end
    
    % compute the prob(accept new technology) values, and prob(stop on or
    % above upper boundary)
    Pupperout = +(muvec >= rawboundtop(i));   % if above upper boundary, prob of stopping above top is 1, otherwise it is 0, at end of time horizon
    if advanced.UnkVariance % handle case of unknown variance and known variance separately
        if advanced.UnkVarianceShape == -1
            enddof = curt;         % use if there is a noninformative prior for mu, sigma^2 of sampling distribution
        else
            enddof = curt + 2*advanced.UnkVarianceShape - t0;    % use if the prior is otherwise specified
        end
        Ppicknewout = TerminalProbPickNewUnk(muvec*numpeopletreated - I, predvarsample, advanced.nochangeallowed, enddof ); % find probability of picking new alternative, given 
    else
        Ppicknewout = TerminalProbPickNew(muvec*numpeopletreated - I, predvarsample,advanced.nochangeallowed); % find probability of picking new alternative, given 
    end
    ENumSampsout = 0*muvec;                    % by default, expected num samps more to do is 0, reset below for case of contin region
    
    % for values in continuation region, compute expectation via
    % conditioning on what happens at time t+dt 
    Pupperout(minindx:maxindx)   = DelayContinExpectation(minindx,maxindx,pu,pd,psame,Pupperin);
    Ppicknewout(minindx:maxindx) = DelayContinExpectation(minindx,maxindx,pu,pd,psame,Ppicknewin);
    if advanced.UnkVariance
        ENumSampsout(minindx:maxindx)= dt + DelayContinExpectation(minindx,maxindx,pu,pd,psame,ENumSampsin);
    else
        ENumSampsout(minindx:maxindx)= dt + DelayContinExpectation(minindx,maxindx,pu,pd,psame,ENumSampsin);
    end
    Pupperin = Pupperout;               % now, update the variables for the recursion
    Ppicknewin = Ppicknewout;
    ENumSampsin = ENumSampsout;
    pcsout = pcsstopnow;
    pcsout(minindx:maxindx) = DelayContinExpectation(minindx,maxindx,pu,pd,psame,pcsstopprev);   

	if (mod(i,advanced.PlotCheck)==1) || (i == tlen-1)
        % This stores the matrices for the output values, and sets a new
        % time vector to manage them all. FIX: maybe some speed
        % optimization can be done here by preallocating these things?
        B0mat = [Cout(:) B0mat];
        Puppermat = [Pupperin(:) Puppermat];
        Pnewmat = [Ppicknewin(:) Pnewmat];
        PCSmat = [pcsstopnow(:) PCSmat];
        if abs(advanced.RegretPenalty) == 0     % this was computed earlier for the case of regretpenalty ~= 0
            if advanced.UnkVariance % handle case of unknown variance and known variance separately
                RewardPI = TerminalRewardFunctionUnk(muvec*numpeopletreated-I,1.0,(numpeopletreated^2 * sigma^2) / (curt-tau),false,enddof);
            else
                RewardPI = TerminalRewardFunction(muvec*numpeopletreated-I,1.0,(numpeopletreated^2 * sigma^2) / (curt-tau),false);
            end
        end
        RewardPIMatII = [RewardPI(:) RewardPIMatII]; % Expected reward with perfect information at time t=Tmax
        RewardMatII = [Rewardmax(:) RewardMatII]; % Expected reward, assuming tau samples pending at time t=Tmax
        tmat = [tvec(i) tmat];

        if advanced.DOPLOT
            figure(3); hold off;
            plot(muvec,ENumSampsout,'-')
            legend('E[Num samples until stop]','Location','NorthWest')
            xlabel('posterior mean');ylabel('E[Num Samples]');
            titlestr=strcat('Time = ',num2str(curt),'; Bound \approx ',num2str(rawboundtop(i)),'; pu = ',num2str(pu));
            title(titlestr);

            figure(4); hold off;
            plot(muvec,Cout,'-',muvec,stoprewardvec,'-.',muvec,(discountfactor^(1-advanced.nochangeallowed))*max(muvec*numpeopletreated - I,0),'--')
            legend('Cout','stopreward','discounted (if needed) max(0,upside)','Location','NorthWest')
            xlabel('posterior mean');ylabel('expected reward');
            titlestr=strcat('Time = ',num2str(curt),'; Bound \approx ',num2str(rawboundtop(i)),'; pu = ',num2str(pu));
            title(titlestr);
            
            figure(5); hold off;
            plot(muvec,Pupperout','-',muvec,Ppicknewout','-.',muvec,pcsstopnow',muvec(minindx)-dmu,1,'x',muvec(maxindx)+dmu,1,'x')
            legend('Pr(stop on upper)','Pr(pick new alternative)','P(CS)','contin set boundary','Location','NorthWest')
            xlabel('posterior mean');ylabel('probability');
            titlestr=strcat('Time = ',num2str(curt),'; Bound \approx ',num2str(rawboundtop(i)),'; pu = ',num2str(pu));
            title(titlestr);

            figure(6); hold off;
            plot(muvec,RewardPI,'-',muvec,Rewardmax,'-.')
            legend('Reward with Perfect Info Now','Reward if contin optimally','Location','NorthWest')
            xlabel('posterior mean');ylabel('expected regret');
            titlestr=strcat('Time = ',num2str(curt),'; Bound \approx ',num2str(rawboundtop(i)),'; pu = ',num2str(pu));
            title(titlestr);

            figure(7); hold off;
            plot(muvec,pcsstopnow,'-',muvec,pcsout,'-.')
            legend('pcsin','pcsout','Location','SouthWest')
            xlabel('posterior mean');ylabel('probability');
            titlestr=strcat('Time = ',num2str(curt),'; Bound \approx ',num2str(rawboundtop(i)),'; pu = ',num2str(pu));
            title(titlestr);

            figure(8)
            plot(muvec,Cintmp-stoprewardvec,'-.',rawboundbot(i) ,0,'x',rawboundtop(i),0,'o')
            xlabel('posterior mean');ylabel('value of continuing');
            titlestr=strcat('Time = ',num2str(curt),'; cont set \approx ',num2str(rawboundbot(i)),' to ',num2str(rawboundtop(i)));
            title(titlestr);
        end
	end
    
    pcsstopprev = pcsout;
	Cin = Cout;
end

%%% COMPUTE BOUNDARY (with smoothing if needed)

if ~advanced.smoothed        % no smoothing desired
    mintvec = tvec;
    bndupper = rawboundtop;
    bndlower = rawboundbot;
    if advanced.DOPLOT
        figure(9)
        plot(tvec,rawboundtop,'--',tvec,rawboundbot,'-.');
        legend('raw stopping boundary','Location','SouthEast');
    end
else                % do smoothing, using the sparser boundary as a basis
    tighten = true;        % false : keep all time values, true: try to reduce number of time values
    if advanced.UnkVariance 
        [mintvec1, bndupper] = LinearSmoothifier( tvec, -rawboundtop, tighten );
        bndupper = - bndupper;
    else
        [mintvec1, bndupper] = LinearSmootherLower( tvec, rawboundtop, tighten );
    end
    [mintvec2, bndlower] = LinearSmootherLower( tvec, rawboundbot, tighten );
    mintvec = sort(unique([mintvec1 mintvec2]));
    bndupper = interp1(mintvec1,bndupper,mintvec,'linear');
    bndlower = interp1(mintvec2,bndlower,mintvec,'linear');
    gridest = max(1,1+round( (min(bndlower)-basic.mumin)/advanced.dmu));
    ENumSampsin(gridest) = max(ENumSampsin(gridest),0.01);
    gridest = min(length(ENumSampsin),1+round( (max(bndupper)-basic.mumin)/advanced.dmu));
    ENumSampsin(gridest) = max(ENumSampsin(gridest),0.01);
    if advanced.DOPLOT
        figure(2)
        plot(tvec,rawboundtop,'-',tvec,rawboundbot,'-');
        hold on
        plot(mintvec,bndupper,'-.',mintvec,bndlower,'-.');        
        legend('raw stopping boundary','smoothed boundary','Location','SouthEast');
    end
end

% Finish setting output values
%B0vec = Cin;
% Build a structure which has all the extra outputs, if desired
mat.B0mat = B0mat;          % value function
mat.Puppermat = Puppermat;  % Probability of stopping due to hitting/exceeding upper boundary
mat.Pnewmat = Pnewmat;      % Probability that if one stopped at a given (t, mu), that new technology is adopted
mat.PCSmat = PCSmat;      % Probability that if one stopped at a given (t, mu), that new technology is adopted
mat.ENumSamps = ENumSampsin;
mat.Bayespcsvec = pcsout; 
mat.tmat = tmat;            % should be a row vector, gives time indices for the above matrices for contour plots

mat.muvec = muvec(:);       % should be a column vector of the mu values

mat.tvec = mintvec;            % tvec is for the boundaries: may be different from tmat, as tvec tries to compress
mat.bndupper = bndupper;    % info requirements to get at curves of boundary, whereas tmat is equally spaced, in order
mat.bndlower = bndlower;    % to help contour plot
mat.B0vec = Cin;                  % first column of B0mat (might trim this as it is extra unneeded info

% Get expected reward with perfect information at time t=Tmax, assuming no discounting
mat.RewardPITmaxIII = RewardPITmaxIII(:);
mat.RewardTmaxIII = min(RewardTmaxIII(:),RewardPITmaxIII(:)); % the 'min' corrects numerical stability issues with RewardTmaxIII
mat.RewardPIMatII = RewardPIMatII; % Expected reward with perfect information at time t in stage II
mat.RewardMatII = RewardMatII; % Expected reward, assuming tau samples pending at time t in stage II

retval = (max(bndupper) < max(muvec)) || (min(bndlower) > min(muvec));

if nonmonotoniclower + nonmonotonicupper > 0
    warning('forced boundaries to retain monotonicity properties');
end

end
