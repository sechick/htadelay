% For use on project: Chick and Frazier, WSC 09 and ensuing paper.
% This file allows one to run numerical tests for:
%   A) with k=1: compute...
%       a) optimal expected reward as estimated from PDE
%       b) upper bound on expected reward, assuming perfect info at no cost
%       c) lower bound on expected reward, using \underline(OER)
%       proposition in paper
%   B) with k > 1: compute...
%       a) upper bound on expected reward, assuming perfect info at no cost
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load in information from PDE solution first. This assumes that testrecurv3.m has been 
% called and that the .mat files with the solution to the PDE in them has been already
% created.
% Other dependencies: 
%   CFBranchGetB1Val: used to load correct .mat as
%       appropriate and then return an estimate of Btilde(w,s).
%   CFBranchSetParams: sets c, alpha, beta, gamma from problem statement
%       parameters.
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% May 25, 2009: created. SEC.
%
% The following two inputs need to be set, based upon the files that are
% output by the program testrecurv3.m

% FIRST, KNOWN VARIANCE CASE

global B1filestring;
global accumsvec;
global accumbound;
global fileindx;
global fileloaded;
fileloaded=-1;
global sbvec;
global wsubvec;
global subBmatrix;
global megawvec;
global megasvec;
global megaBmatrix;

mysmallfontsize = 18;

B1filestring = 'BwsMat';  % test file - now uses 'fix' in iteration with bounds.
StartFileVal = 3;       % The smallest file index of name 'NewBnd*.mat' to be loaded
EndFileVal = 9;        % The largest file index of name 'NewBnd*.mat' to be loaded

[accumsvec accumbound fileindx] = CFBranchLoadFiles(B1filestring, StartFileVal, EndFileVal);

PLOTFIG = 101
figure(PLOTFIG)
semilogy(1./accumsvec,accumbound);
set(gca,'FontSize',mysmallfontsize);
tmp=axis;
tmp(2) = 1/min(accumsvec);
tmp(4) = max(accumbound);
xlabel('tau');ylabel('bound(tau)');
axis(tmp);
%hold on
figure(PLOTFIG+1)
loglog(accumsvec,accumbound);
xlabel('s');ylabel('bound(s)');

approxbnd = ApproxBoundW(accumsvec);

PLOTFIG = 105
figure(PLOTFIG)
semilogy(1./accumsvec,accumbound,1./accumsvec,approxbnd);
tmp=axis;
tmp(2) = 1/min(accumsvec);
tmp(4) = max(accumbound);
xlabel('tau');ylabel('bound(tau)');
axis(tmp);
%hold on
figure(PLOTFIG+1)
loglog(accumsvec,accumbound,accumsvec,approxbnd);
xlabel('s');ylabel('bound(s)');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% KNOWN VARIANCE CASE

% Compute the PDE approximation to Bwiggle(w,s) and translate that back
% into (y/t, t) coordinates.

PLOTBASE = 201
PLOTON = 0 %1
MYEPS=10^-6             % set a lower bound for m, for calculation purposes, to avoid numerical stability issues at s=0

hourcost = 1; %10;      
sigma = 10^5; %10^6;     % std dev during one time unit
mu0 = 0;
t0 = 1
sig0 = sigma / sqrt(t0) ; % for prior distr of unknown means
cpumin = 60;       % cpu time to get output of that unit
m = 0

[c gamma beta alpha] = CFBranchSetParams(hourcost, cpumin, sigma);

figure(PLOTBASE+1)
semilogy((1/gamma)./accumsvec, m + accumbound/beta,'--');
set(gca,'FontSize',mysmallfontsize);

%plot((1/gamma)./accumsvec, m + accumbound/beta,'--',(1/gamma)./accumsvec, m + accumbound/beta,'--')
xlabel('t');ylabel('bound(t)');


'tests with k = 1'
mu0vec = [0 0 0]
t0vec = [1 10 100]
Bvec = 0*mu0vec; % initialize array length
nreps=1:.01:5000;

for i=1:length(mu0vec)
    mu0 = mu0vec(i);
    t0 = t0vec(i);
    Bvec(i) = CFBranchGetB1Val(beta*mu0,1/(gamma*t0))/beta;     % Get solution in terms of mean of X
    OERupper(i) = (sigma/sqrt(t0)) * PsiNorm( abs(mu0-m) / (sigma/sqrt(t0))); % perfect info at no cost
    predstd = sigma*sqrt(nreps ./ ( t0 * (t0 + nreps)));
    predtest = predstd .* PsiNorm( abs(mu0-m) ./ predstd ) - c * nreps;
    OERlower(i) = max(predtest);                                % Use the lemma to get a lower bound
end
Bvec
OERupper
OERlower

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% KNOWN VARIANCE CASE

% Monte Carlo estimates of upper bound for optimal expected reward (which
% comes from assuming perfect information without any sampling costs)

'tests with k > 1'
sigma=10^5
m = 0
mu0 = 0
t0 = 1
%kvec = [2 3 5 10 20 50]
kvec = [2 3 5 10 20]
OERupperk = 0 * kvec; % for initialization
OERupperkSE = 0 * kvec; % for initialization
betatst = 10.^(0:.1:4);

NMAXVALS = 10^7;
NMAXVALS = 10^6;
for i = 1:length(kvec)
    k = kvec(i);
    numruns = floor(NMAXVALS / (k+1))
    bigmatrix = normrnd(mu0,sigma/sqrt(t0),k+1,numruns);
    bigmatrix(1,:) = m*ones(1,numruns);  % include the m option as an alternative
    macroreps=max(bigmatrix);
    macroreps2= m + (m - min(bigmatrix)); % antithetic variates
    testvals = (macroreps+macroreps2) / 2;
    OERupperk(i) = mean(testvals);
    OERupperkSE(i) = std(testvals)/sqrt(length(testvals));
    
    betavec=0*betatst;  % initialize
    for j=1:length(betatst)
        betavec(j) = mean( max(mu0 + (macroreps-mu0)*sqrt(betatst(j)/(t0+betatst(j))), m) - c*k*betatst(j));
    end
    [ta tb] = max(betavec);
    OERlowerk(i) = ta;
    OERnruns(i) = tb*k;
end
OERupperk
OERupperkSE
OERlowerk
OERnruns
OERupperk - OERlowerk
(OERupperk - OERlowerk)./OERupperk

%%%%
%sigma=10^5
%m = 0
%mu0 = 0
t0 = 100
%kvec = [2 3 5 10 20 50]
OERupperk = 0 * kvec; % for initialization
OERlowerk = 0 * kvec; % for initialization
OERupperkSE = 0 * kvec; % for initialization
betatst = 10.^(0:.1:4);

%NMAXVALS = 10^7;
NMAXVALS = 10^6;
for i = 1:length(kvec)
    k = kvec(i);
    numruns = floor(NMAXVALS / (k+1))
    bigmatrix = normrnd(mu0,sigma/sqrt(t0),k+1,numruns);
    bigmatrix(1,:) = m*ones(1,numruns);  % include the m option as an alternative
    macroreps=max(bigmatrix);
    macroreps2= m + (m - min(bigmatrix)); % antithetic variates
    testvals = (macroreps+macroreps2) / 2;
    OERupperk(i) = mean(testvals);
    OERupperkSE(i) = std(testvals)/sqrt(length(testvals));
    
    betavec=0*betatst;  % initialize
    for j=1:length(betatst)
        betavec(j) = mean( max(mu0 + (macroreps-mu0)*sqrt(betatst(j)/(t0+betatst(j))), m) - c*k*betatst(j));
    end
    [ta tb] = max(betavec);
    OERlowerk(i) = ta;
    OERnruns(i) = tb*k;
end
OERupperk
OERupperkSE
OERlowerk
OERnruns
OERupperk - OERlowerk
(OERupperk - OERlowerk)./OERupperk

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% FOR SECTION 6 OF CHICK AND FRAZIER PAPER, 2009

myfontsize=22
mysmallfontsize=18

% UNKNOWN VARIANCE CASE: 1 DO PLOTS

PLOTBASE = 301
PLOTON = 0 %1
MYEPS=10^-6             % set a lower bound for m, for calculation purposes, to avoid numerical stability issues at s=0


Nreps = 3;
retirement = 0;

desiredscaleforvariance = 10
sigmaformeanofvar = 10^5
meanofvar = (sigmaformeanofvar)^2
varofvar = meanofvar^2 / (desiredscaleforvariance-2); 
effectivenumreps = 5
eta = effectivenumreps

mu0 = 0

shape = 2 + meanofvar^2 / varofvar  % back fit these inverted gamma params using Bernardo & Smith, p. 431
scale = (meanofvar * (shape - 1))   % noting that for matlab, the mean of a gamma is shape*scale...
vecvar = 1./gamrnd(shape*ones(1,Nreps),1/scale);
mean(vecvar)
var(vecvar)
vecmean = normrnd(mu0,sqrt(vecvar / eta));
mean(vecmean)
var(vecmean)
meanofvar/eta       % this and the value of var(vecmean) should be close
maxOER = mean(max(vecmean,retirement))

% The next bit is an alternative way to get student t for unknown mean - as
% a check for prior code. this is not needed, as we also need the variance
altvect = mu0 + sqrt(scale/eta/(shape)) * trnd(2*shape,1,Nreps);    % again, shape*scale instead of scale/shape
%tst1 = mu0 + trnd(2*shape,1,Nreps);    % again, shape*scale instead of scale/shape
mean(altvect)
var(altvect)
altmaxOER = mean(max(altvect,retirement))

'tests with k = 1'
vecvar = 1./gamrnd(shape*ones(1,Nreps),1/scale);
vecmean = normrnd(mu0,sqrt(vecvar / eta));

MAXREPS = 5000 % maximum number of samples to take
tmprepvec = ones(MAXREPS,1);
nsamples = normrnd(tmprepvec*vecmean,tmprepvec*sqrt(vecvar));   %matrix: is MAXREPS x Nreps

tvec = eta + (1:MAXREPS);
numsamps = (1:MAXREPS)';    %row vector
sumx = cumsum(nsamples);        % put
posteriormean = (mu0*eta + sumx) ./ ((eta + numsamps)*ones(1,Nreps));
postmean = zeros(size(nsamples));
postvar = postmean;

hourcost = 1; %10;      
sigma = sqrt(meanofvar); %10^6;     % std dev during one time unit
%sig0 = sigma / sqrt(t0) ; % for prior distr of unknown means
cpumin = 60;       % cpu time to get output of that unit
%m = 0
[c gamma beta alpha] = CFBranchSetParams(hourcost, cpumin, sigma);

cfactor = c^(1/3)
for i=1:Nreps
       % Initialize parameters for inference
    myscale = scale;
    myshape = shape;
    myeta = eta;
    mymu0 = mu0;
    samples(i) = 0;
    trueboundapprox = ApproxBoundW(1/gamma./tvec)/beta;
    boundapprox = trueboundapprox;
        % do the inference until maximum number of samples observed, or
        % until outside stopping set
    inContinuation = true;
    for j=1:MAXREPS
        % observe next sample for run
        myeta = myeta + 1;
        myshape = myshape + 1/2;
        myscale = myscale + (nsamples(j,i) - mymu0)^2 * myeta / (myeta+1) / 2;
        mymu0 = mymu0 + (nsamples(j,i) - mymu0)/(myeta);
        postmean(j,i) = mymu0;
        postvar(j,i) = myscale / (myshape - 1);
        sval = postvar(j,i)^(1/3) / cfactor^2 / myeta;
        bndapprox = cfactor * postvar(j,i)^(1/3) * ApproxBoundW(sval);
        boundapprox(j) = bndapprox;
    end
    EstimatedBoundDuringRun = cfactor * postvar(:,i).^(1/3) .* ApproxBoundW(postvar(:,i).^(1/3) / cfactor^2 ./ tvec');
    truevar=vecvar(i)
    boundwithvarknown = cfactor * truevar^(1/3) * ApproxBoundW(truevar^(1/3) / cfactor^2 ./ tvec');

%    svec = postvar(:,i).^(1/3) ./ cfactor^2 ./ tvec';
%    boundapprox = cfactor * postvar(:,i).^(1/3) .* ApproxBoundW(svec);
    figure(PLOTBASE+i)
    plot(tvec,boundapprox,'--',tvec,postmean(:,i),'-',tvec,boundwithvarknown,'-.',tvec,-boundapprox,'--')
    tmp=axis;tmp(1)=min(tvec);tmp(2)=max(tvec)/5;axis(tmp);
    set(gca,'FontSize',mysmallfontsize);
    legend('Estimated Boundary','Posterior Mean','Boundary if Realized Variance Were Known');
    xlabel('Effective number of samples','FontSize',myfontsize,'FontName','Times'); 
    ylabel('Posterior mean of unknown mean','FontSize',myfontsize,'FontName','Times')
%    title('Approximation of continuation set with a sample path of posterior mean','FontSize',myfontsize,'FontName','Times')
    pause
    mytitle = strcat(strcat('UnkVarSampBound',int2str(i)),'.eps');
	print('-deps',mytitle);	

    figure(PLOTBASE+Nreps+i)
    meanofvar
    boundwithmeanvar = cfactor * meanofvar^(1/3) * ApproxBoundW(meanofvar^(1/3) / cfactor^2 ./ tvec');
    loglog(tvec,EstimatedBoundDuringRun,'--',tvec,boundwithvarknown,'-',tvec,boundwithmeanvar,'-.')
    tmp=axis;tmp(1)=min(tvec);tmp(2)=max(tvec);axis(tmp);
    set(gca,'FontSize',mysmallfontsize);
    legend('Estimated Boundary','Boundary if Realized Variance Were Known','Boundary With Mean Variance')
    xlabel('Effective number of samples','FontSize',myfontsize,'FontName','Times'); 
    ylabel('Posterior mean of unknown mean','FontSize',myfontsize,'FontName','Times');
    pause
    mytitle = strcat(strcat('CompareContinSets',int2str(i)),'.eps');
    print('-deps',mytitle);	
end

figure(PLOTBASE+1+2*Nreps)
plot(tvec,postvar(:,i));
set(gca,'FontSize',mysmallfontsize);
tmp=axis;tmp(1)=min(tvec);tmp(2)=max(tvec);axis(tmp);
xlabel('Effective number of samples, \eta','FontSize',myfontsize,'FontName','Times'); 
ylabel('Posterior mean of unknown variance','FontSize',myfontsize,'FontName','Times');
pause
mytitle = strcat(strcat('UnkVarSampVar',int2str(i)),'.eps');
print('-deps',mytitle);	


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% FOR SECTION 6 OF CHICK AND FRAZIER PAPER, 2009

% UNKNOWN VARIANCE CASE: 2 NO PLOTS - GO FOR MORE NUMBERS

PLOTBASE = 321
PLOTON = 0 %1
MYEPS=10^-6             % set a lower bound for m, for calculation purposes, to avoid numerical stability issues at s=0

Nreps = 40000;
Nreps = 100000;

retirement = 0;
desiredscaleforvariance = 10
sigmaformeanofvar = 10^5
meanofvar = (sigmaformeanofvar)^2
varofvar = meanofvar^2 / (desiredscaleforvariance-2); 
effectivenumreps = 5
eta = effectivenumreps

mu0 = 0

shape = 2 + meanofvar^2 / varofvar  % back fit these inverted gamma params using Bernardo & Smith, p. 431
scale = (meanofvar * (shape - 1))   % noting that for matlab, the mean of a gamma is shape*scale...
vecvar = 1./gamrnd(shape*ones(1,Nreps),1/scale);
mean(vecvar)
var(vecvar)
vecmean = normrnd(mu0,sqrt(vecvar / eta));
mean(vecmean)
var(vecmean)
meanofvar/eta       % this and the value of var(vecmean) should be close


'tests with k = 1'
samples = zeros(1,Nreps); simreward = samples; oppcost = samples; simnetreward = samples; 
    % allocate memory for total number of samples, fo
    % for reward upon stopping not counting sampling costs, and for total
    % reward including cost of sampling

hourcost = 1; %10;      
sigma = sqrt(meanofvar); %10^6;     % std dev during one time unit
%sig0 = sigma / sqrt(t0) ; % for prior distr of unknown means
cpumin = 60;       % cpu time to get output of that unit

[c gamma beta alpha] = CFBranchSetParams(hourcost, cpumin, sigma);

cfactor = c^(1/3);
for i=1:Nreps
    if (mod(i,10000)==0) 
        i
    end
       % Initialize parameters for inference
    myscale = scale;
    myshape = shape;
    myeta = eta;
    mymu0 = mu0;
    samples(i) = 0;
        % do the inference until maximum number of samples observed, or
        % until outside stopping set
    inContinuation = true;
    while inContinuation
        samples(i) = samples(i) + 1;
        Xij = normrnd(vecmean(i),sqrt(vecvar(i)));
        % observe next sample for run
        myeta = myeta + 1;
        myshape = myshape + 1/2;
        myscale = myscale + (Xij - mymu0)^2 * myeta / (myeta+1) / 2;
        mymu0 = mymu0 + (Xij - mymu0)/(myeta);
        postmean = mymu0;
        postvar = myscale / (myshape - 1);
        sval = postvar^(1/3) / cfactor^2 / myeta;
        bndapprox = cfactor * postvar^(1/3) * ApproxBoundW(sval); %* (myeta / myeta-2);
        if postmean > retirement + bndapprox    % went above continuation region
            simreward(i) = postmean;
            oppcost(i) = max(vecmean(i),retirement) - simreward(i);
            simnetreward(i) = postmean - c*samples(i);
            inContinuation = false;
        elseif postmean < retirement - bndapprox % wend below continuation region
            simreward(i) = retirement;
            oppcost(i) = max(vecmean(i),retirement) - retirement;
            simnetreward(i) = retirement - c*samples(i);
            inContinuation= false;
        end % else, still in continuation region, get another sample
    end
end
maxOERwithmean=mean(max(retirement,vecmean))

meanreward=mean(simreward)
stder=std(simreward)/sqrt(length(simreward))
relerrrew = stder/meanreward

meanoppcost=mean(oppcost)
stder=std(oppcost)/sqrt(length(oppcost))
relerropcost = stder/meanoppcost

meansamples=mean(samples)
stder=std(samples)/sqrt(length(samples))
relerrsamples = stder/meansamples

meanNetReward=mean(simnetreward)
stder=std(simnetreward)/sqrt(length(simnetreward))
relerrnet = stder/meanNetReward

prec = scale / shape / eta;
maxOERwithT=sqrt(prec)*PsiNormUV((mu0-retirement)/sqrt(prec),2*shape )

%best one-stage sampling
nreps=1:1:3000;
predstd = sqrt(scale * nreps ./ ( eta * shape * (eta + nreps)));
predtest = predstd .* PsiNormUV( abs(mu0-retirement) ./ predstd,2*shape ) - c * nreps;
[minOER indx ] = max(predtest)        ;                        % Use the lemma to get a lower bound
minOER
nreps(indx)
maxOERwithT - minOER - c*nreps(indx)

'tests with k = 1'
Bvec = 0*vecvar(1:min(length(vecvar),5000)); % initialize array length
maxOERvec = Bvec;
minOERvec = Bvec;

for i=1:length(Bvec)
    gamma = (cfactor^2 / vecvar(i) ) ^(1/3);
    beta = 1 / (cfactor * vecvar(i))^(1/3);
    Bvec(i) = CFBranchGetB1Val(beta*mu0,1/(gamma*eta))/beta;     % Get solution in terms of mean of X
%    OERupper(i) = (sigma/sqrt(t0)) * PsiNormUV( abs(mu0-m) / (sigma/sqrt(t0))); % perfect info at no cost
    maxOERvec(i) = sqrt(vecvar(i)/eta) * PsiNorm( abs(mu0-retirement) / sqrt(vecvar(i)/eta) );
    predstd = sqrt(vecvar(i)*nreps ./ ( eta * (eta + nreps)));
    predtest = predstd .* PsiNorm( abs(mu0-retirement) ./ predstd ) - c * nreps;
    minOERvec(i) = max(predtest);                                % Use the lemma to get a lower bound
end
MeanWithAverageOverUnknownGivenOptimalForKnown=mean(Bvec)
maxOERval = mean(maxOERvec)
minOERval = mean(minOERvec)

%values with stopping rule based upon estimated optimal stopping boundary.
maxOERwithT - meansamples - meanoppcost
std(oppcost+samples)/sqrt(length(samples))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% FOR SECTION 6 OF CHICK AND FRAZIER PAPER, 2009

% UNKNOWN VARIANCE: k > 1 MULTIPLE SYSTEMS

'tests with k > 1'
retirement = 0;
desiredscaleforvariance = 10
sigmaformeanofvar = 10^5
meanofvar = (sigmaformeanofvar)^2
varofvar = meanofvar^2 / (desiredscaleforvariance-2); 
effectivenumreps = 5
eta = effectivenumreps

mu0 = 0

shape = 2 + meanofvar^2 / varofvar  % back fit these inverted gamma params using Bernardo & Smith, p. 431
scale = (meanofvar * (shape - 1))   % noting that for matlab, the mean of a gamma is shape*scale...
c=1

%kvec = [2 3 5 10 20 50]
kvec = [2 3 5 10 20 50 100 ]
OERupperk = 0 * kvec; % for initialization
OERupperkSE = 0 * kvec; % for initialization
%betatst = 1:.1:2000;
betatst = 10.^(0:.025:4);


NMAXVALS = 3*10^6;        % maximum size of vector of means
for i = 1:length(kvec)
    k = kvec(i);
    numruns = floor(NMAXVALS / (k+1));  % number of macroreplications to do
    bigmatrix = mu0 + sqrt(scale/eta/shape) * trnd(2*shape,k+1,numruns);  %was normrnd(mu0,sigma/sqrt(t0),k+1,numruns);   % generate the unknown means. 
    bigmatrix(1,:) = retirement*ones(1,numruns);  % include the m option as an alternative
    macroreps=max(bigmatrix);
    macroreps2= retirement + (retirement - min(bigmatrix)); % antithetic variates
    testvals = (macroreps+macroreps2) / 2;
    OERupperk(i) = mean(testvals);
    OERupperkSE(i) = std(testvals)/sqrt(length(testvals));
    
    betavec=0*betatst;  % initialize
    for j=1:length(betatst)
        betavec(j) = mean( max(mu0 + (macroreps-mu0)*sqrt(betatst(j)/(eta+betatst(j))), retirement) - c*k*betatst(j));
    end
    [ta tb] = max(betavec);
    OERlowerk(i) = ta;
    OERnruns(i) = betatst(tb)*k;
end
OERupperk
OERupperkSE
OERlowerk
toOERnruns
OERupperk - OERlowerk
(OERupperk - OERlowerk)./OERupperk
OERupperkSE./OERupperk
clear bigmatrix;
OERupperk-OERlowerk-OERnruns

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

myfontsize=16
mysmallfontsize=14

% UNKNOWN VARIANCE CASE: k >= 1 NUMBERS

%EOCckStop = true    % SET TO TRUE TO HAVE EOC_{ck} STOP RULE, AND FALSE TO HAVE EOC_{c1} STOP RULE
EOCckStop = 1   ;
EOC1kStop = 2   ;
KGstarstop = 3  ;
KGonestop = 4   ;
GitStop = 5     ;    % This is the PDE stopping rule
StopRule = GitStop    % Pick a stopping rule.

LLAlloc = 1     ;
KGAlloc = 2     ;
KGOneAlloc = 3  ;
GitAlloc = 4    ;    % This is the PDE allocation rule
AllocRule = KGOneAlloc    % Pick an allocation rule.

MAXBETASTAR = 10^4

PLOTBASE = 700
PLOTON = 0 %1
MYEPS=10^-9             % set a lower bound for m, for calculation purposes, to avoid numerical stability issues at s=0

Nreps = 10;
retirement = 0;

nreps=round([1; (2:2:10)' ; 10.^(1.05:.1:3.5)']);      % these determine the betas to test for the KG*

mu0=0
desiredscaleforvariance = 10
sigmaformeanofvar = 10^5
meanofvar = (sigmaformeanofvar)^2
varofvar = meanofvar^2 / (desiredscaleforvariance-2); 
effectivenumreps = 5
eta = effectivenumreps

shape = 2 + meanofvar^2 / varofvar  % back fit these inverted gamma params using Bernardo & Smith, p. 431
scale = (meanofvar * (shape - 1))   % noting that for matlab, the mean of a gamma is shape*scale...
c=1
cpumin=60;

kvec = [1 2 3 5 10]
kvec = [2 3 5 10 20 50 100 ]
%kvec = [100]
truemeanReward = 0 * kvec;
SEtruemeanReward = 0 * kvec;
meanReward = 0 * kvec;
SEmeanReward = 0 * kvec;
meanOppCost = 0 * kvec;
SEmeanOppCost = 0 * kvec;
meanBudgetUsed = 0*kvec;
SEmeanBudgetUsed = 0*kvec;
cputim = 0*kvec;
varmeanReward = 0*kvec;
vartruemean = 0*kvec;
covmeanandtrue = 0*kvec;
OverallNet = 0*kvec;

NMAXVALS = 1000;        % maximum size of vector of means
NMAXVALS = 200;        % maximum size of vector of means
NMAXVALS = 15000;        % maximum size of vector of means
for i = 1:length(kvec) %1:length(kvec)
    k = kvec(i);
    numruns = floor(NMAXVALS / k);  % number of macroreplications to do
    selectedvalue = zeros(numruns,1);
    truemeanselected = zeros(numruns,1);
    oppcostvec = zeros(numruns,1);
    budgetspent = zeros(numruns,1);
    totbudgetvec = [3 5 10.^(1.05:.1:3.0) 1 7 9];    % total budget can depend upon k
%    totbudgetvec = round(10.^(0:.2:3));
    tic;
    RedoBeta = -1;  % set to negative value of all the betastars need recalculation,
                    % set to a positive index if only that one betastar
                    % needs recalculation
    
    for ell=1:numruns %1:numruns 
        %initialize prior distributoin for alternatives , with special attention to 'known' alternative
        mu0vec = mu0 * ones(k+1,1); mu0vec(k+1)=retirement;
        scalevec = scale * ones(k+1,1);
        shapevec = shape * ones(k+1,1);
        scalevec(k+1) = scalevec(k+1) / 10^20;  % make the variance for the 'known' system very small
        etavec = effectivenumreps*ones(k+1,1); etavec(k+1) = 10^5; % 
        notdone = true;
        cvec = c*ones(k+1,1); %cvec(k+1) = 10^3;

        truevarvec = 1./gamrnd(shapevec,1./scalevec);
        truemeanvec = normrnd(mu0vec,sqrt(truevarvec./etavec));
        truemeanvec(k+1) = retirement;
        
        predstartest = zeros(k+1,1);
        betastarvec = 10*ones(k+1,1);

        fractst = [.01 .02 .05 .1 .2 .4 .8];
        Curbudget = totbudgetvec(1);
        while notdone       % iterate until stopping rule satisfied
            % pick 'best' mean so far and second best mean
            [bstmean bstindx] = max(mu0vec);
            tmpvec = mu0vec; tmpvec(bstindx) = min(mu0vec) - 1; 
            [scndmean scndindx] = max(tmpvec);
            testvec = bstmean*ones(k+1,1); % create vector of 'best of the others'
            testvec(bstindx) = scndmean; % used for test statistics
            diffbestvec = abs(mu0vec - testvec);     % use formula from prop in appendix b of paper, using plug in estimator for variance

            if (StopRule == KGstarstop) || (AllocRule == KGAlloc)
                if RedoBeta < 0            % if all betastars need recalculating...
                    uvec = (abs(testvec-mu0vec) .* etavec ./ (scalevec ./ (shapevec-1))).^2;
                    betastartstvec = (etavec / 4) .* (uvec -1 + sqrt((uvec - 1).^2 + 6*uvec + 1 ) );
%                    fractst=0.9;
%                    betastarvec2 = etavec .* ( -1 + (etavec + uvec/fractst) .* normpdf(sqrt(uvec/fractst),etavec) ./ (2*normcdf(-sqrt(uvec/fractst),etavec)) ./ (etavec - 1) );
%                    betastarvec2 = etavec .* ( -1 + normpdf(sqrt(uvec/fractst),etavec) ./ (2*normcdf(-sqrt(uvec/fractst),etavec)) );
%                    betastarvec2 = etavec .* ( -1 + tpdf(sqrt(uvec/fractst),etavec) ./ (2*tcdf(-sqrt(uvec/fractst),etavec)) );
                    for ijk=1:k
%[ta tb tc]=fminsearch('PsifuncPredOpt',10,optimset('Display','off','TolFun',1,'MaxFunEvals',20),abs(testvec(ijk)-mu0vec(ijk)),scalevec(ijk)/shapevec(ijk), etavec(ijk), 2*shapevec(ijk),cvec(ijk));
%                        tmpval = max(1,abs(testvec(ijk)-mu0vec(ijk)) / sqrt(scalevec(ijk)/shapevec(ijk)));
%                        nreps = [1 4 9 (1:.5:2.5)*etavec(ijk)*(tmpval^2 + tmpval) max(betastarvec(ijk)-1,1)]';
                        betastarvec2 = etavec(ijk) .* ( -1 + normpdf(sqrt(uvec(ijk)./fractst)) ./ (2*normcdf(-sqrt(uvec(ijk)./fractst))) );
                        nreps = max(1,[betastarvec(ijk) betastarvec2]);
%                        nreps = max(1,betastarvec2(ijk)*(.8:.05:1.2));
                        if (max(nreps) == min(nreps))
                            nreps = nreps(1)*ones(1,1);
                        end
                        lamvecijk=sqrt((scalevec(ijk)/shapevec(ijk)/etavec(ijk)) .* nreps ./ (etavec(ijk) + nreps));
                        tmpvec = lamvecijk .* PsiNormUV(abs(testvec(ijk)-mu0vec(ijk))./lamvecijk,etavec(ijk)) - cvec(ijk)*nreps;
                        tmpvec = tmpvec ./ nreps ./ cvec(ijk);
                        [ta tb] = max(tmpvec);
                        if ta < 0   % if we don't have a result with positive value, then...
                            [ta tb tc]=fminbnd('PsifuncPredOpt',1,1500,optimset('Display','off','TolX',1,'MaxFunEvals',10),abs(testvec(ijk)-mu0vec(ijk)),scalevec(ijk)/shapevec(ijk), etavec(ijk), 2*shapevec(ijk),cvec(ijk));
                            predstartest(ijk) = -tb;
                            betastarvec(ijk) = ta;
                        else
                            predstartest(ijk) = ta;
                            betastarvec(ijk) = nreps(tb);
                        end

                        lamvecijk=sqrt((scalevec./shapevec./etavec) .* betastartstvec ./ (etavec + betastartstvec));
                        tmpvec = lamvecijk .* PsiNormUV(abs(testvec-mu0vec)./lamvecijk,etavec) - cvec.*betastartstvec;
                        tmpvec = tmpvec ./ betastartstvec ./ cvec;

                    end
                    betastarvec(k+1) = -1;
                    predstartest(k+1) = min(predstartest)-1;  % prevent the known alternative from being selected
                elseif RedoBeta <= k        % only one betastar needs recalculating
                    ijk = RedoBeta;
%[ta tb tc]=fminsearch('PsifuncPredOpt',10,optimset('Display','off','TolFun',1,'MaxFunEvals',20),abs(testvec(ijk)-mu0vec(ijk)),scalevec(ijk)/shapevec(ijk), etavec(ijk), 2*shapevec(ijk),cvec(ijk));
%                    tmpval = max(1,abs(testvec(ijk)-mu0vec(ijk)) / sqrt(scalevec(ijk)/shapevec(ijk)));
                    uvec(ijk) = (abs(testvec(ijk)-mu0vec(ijk)) .* etavec(ijk) ./ (scalevec(ijk) ./ (shapevec(ijk)-1))).^2;
                    betastartstvec(ijk) = (etavec(ijk) / 4) .* (uvec(ijk) -1 + sqrt((uvec(ijk) - 1).^2 + 6*uvec(ijk) + 1 ) );
%                    nreps = [1 4 9 (1:.5:2.5)*etavec(ijk)*(tmpval^2 + tmpval) max(betastarvec(ijk)-1,1)]';
%                    nreps = max(1,betastartstvec(ijk)*[1 1.2 1.44 1.9]);
%                    nreps = max(1,betastarvec2(ijk)*(.8:.05:1.2));
                    betastarvec2 = etavec(ijk) .* ( -1 + normpdf(sqrt(uvec(ijk)./fractst)) ./ (2*normcdf(-sqrt(uvec(ijk)./fractst))) );
                    nreps = max(1,[betastarvec(ijk) betastarvec2]);
                    if (max(nreps) == min(nreps))
                        nreps = nreps(1)*ones(1,1);
                    end
                    lamvecijk=sqrt((scalevec(ijk)/shapevec(ijk)/etavec(ijk)) .* nreps ./ (etavec(ijk) + nreps));
                    tmpvec = lamvecijk .* PsiNormUV(abs(testvec(ijk)-mu0vec(ijk))./lamvecijk,etavec(ijk)) - cvec(ijk)*nreps;
                    tmpvec = tmpvec ./ nreps ./ cvec(ijk);
                    [ta tb] = max(tmpvec);
                    if ta < 0   % if we don't have a result with positive value, then...
                      [ta tb tc]=fminbnd('PsifuncPredOpt',1,1500,optimset('Display','off','TolX',1,'MaxFunEvals',10),abs(testvec(ijk)-mu0vec(ijk)),scalevec(ijk)/shapevec(ijk), etavec(ijk), 2*shapevec(ijk),cvec(ijk));
                        predstartest(ijk) = -tb;
                        betastarvec(ijk) = ta;
                    else
                        predstartest(ijk) = ta;
                        betastarvec(ijk) = nreps(tb);
                    end
%                    lamvecijk=sqrt((scalevec(ijk)/shapevec(ijk)/etavec(ijk)) .* nreps ./ (etavec(ijk) + nreps));
%                    tmpvec = lamvecijk .* PsiNormUV(abs(testvec(ijk)-mu0vec(ijk))./lamvecijk,etavec(ijk)) - cvec(ijk)*nreps;
%                    [ta tb] = max(tmpvec);
%                    predstartest(ijk) = ta;
%                    betastarvec(ijk) = nreps(tb);
                else
                    'invalid value of redobeta'
                end
            end
            if (StopRule == KGonestop) || (AllocRule == KGOneAlloc)
                betaonevec = ones(size(diffbestvec)); 
                predonestd = sqrt(scalevec .* betaonevec ./ ( etavec .* shapevec .* (etavec + betaonevec)));
                predonetest = predonestd .* PsiNormUV( abs(testvec-mu0vec) ./ max(10^-200,predonestd), 2*shapevec ) ./ (cvec .* betaonevec) ;
                predonetest(k+1) = min(predonetest)-1;
            end
            if (StopRule == GitStop) || (AllocRule == GitAlloc)
                [ctmp gavec bevec alvec] = CFBranchSetParams(cvec, cpumin, sqrt(scalevec./(shapevec-1)));
                bndapprox = ApproxBoundW(1./gavec./etavec)./bevec;
                bndapprox(k+1) = 0;
                okformore = bndapprox - abs(mu0vec-testvec);   % max(bndapprox - abs(mu0vec-testvec),0);
                okformore(k+1) = min(okformore) - 1;        %make sure that system k+1, the retirement option, is not selected for samples
                AllocTestValue = okformore .* bevec;
    %%% START TEST CODE FOR PETER'S ALLOC RULE: THE ONE THAT WOULD CONTINUE FOR THE LARGEST c
                cmin = c; numeligmin = 0;
                while numeligmin < 1
                    cmin = cmin / 2;
                    [ctmp gavec bevec alvec] = CFBranchSetParams(cmin, cpumin, sqrt(scalevec./(shapevec-1)));
                    bndapprox = ApproxBoundW(1./gavec./etavec)./bevec; bndapprox(k+1) = 0;
                    mintestvec = bndapprox - abs(mu0vec - testvec);
                    numeligmin = sum( mintestvec > 0 );
                end
                cmax = c; numeligmax = 2;       % find a "big" c such that no more than one alternative would be chosen
                while numeligmax > 1
                    cmax = cmax * 2;
                    [ctmp gavec bevec alvec] = CFBranchSetParams(cmax, cpumin, sqrt(scalevec./(shapevec-1)));
                    bndapprox = ApproxBoundW(1./gavec./etavec)./bevec; bndapprox(k+1) = 0;
                    numeligmax = sum( bndapprox > abs(mu0vec - testvec));
                end
                cmid = (cmin + cmax) / 2;
                [ctmp gavec bevec alvec] = CFBranchSetParams(cmid, cpumin, sqrt(scalevec./(shapevec-1)));
                bndapprox = ApproxBoundW(1./gavec./etavec)./bevec; bndapprox(k+1) = 0;
                iseligiblevec = (bndapprox - abs(mu0vec - testvec));
                numeligmid = sum( iseligiblevec > 0);
                ccheckindx=0;
                while (numeligmid ~= 1)
                    ccheckindx=ccheckindx+1;
                    if numeligmid > 1
                        cmin = cmid;
                        mintestvec = iseligiblevec;
                    elseif numeligmid < 1
                        cmax = cmid;
                    end
                    cmid = (cmin + cmax) / 2; 
                    [ctmp gavec bevec alvec] = CFBranchSetParams(cmid, cpumin, sqrt(scalevec./(shapevec-1)));
                    bndapprox = ApproxBoundW(1./gavec./etavec)./bevec; bndapprox(k+1) = 0;
                    iseligiblevec = bndapprox - abs(mu0vec - testvec);
                    numeligmid = sum( iseligiblevec > 0);
                    if ccheckindx > 20
                        break
                    end
                end
                AllocTestValue = iseligiblevec';
                if max(AllocTestValue) <= 0
                    AllocTestValue = mintestvec;
                end
    %            ccheckindx
    %%% END TEST CODE FOR PETER'S ALLOC RULE: THE ONE THAT WOULD CONTINUE FOR THE LARGEST c
            end

            if AllocRule == LLAlloc
                [alloc, mydof] = calcrepslinminunkvar(Curbudget,-mu0vec',(scalevec./shapevec)',etavec',1); 
            elseif (AllocRule == KGAlloc)
                [ta, tb] = max( predstartest );
                alloc = zeros(1,k+1);
                alloc(tb) = 1;    % allocate one replication that has the highest knowledge gradient (over budgets tested)
                mydof = 2*shapevec';
            elseif AllocRule == KGOneAlloc
                [ta, tb] = max( predonetest );
                alloc = zeros(1,k+1);
                alloc(tb) = 1;    % allocate one replication that has the highest knowledge gradient (over budgets tested)
                mydof = 2*shapevec';
            elseif AllocRule == GitAlloc
                [ta, tb] = max( AllocTestValue );
                alloc = zeros(1,k+1);
                if ta == 0
                    tb = min(k,max(1,ceil(k*rand)));        % pick system randomly if Git stopping rule would say to stop
                end
                alloc(tb) = 1;    % allocate one replication that has the highest knowledge gradient (over budgets tested)
                mydof = 2*shapevec';
            end
            
            if StopRule == EOCckStop   % sets stopping criterion
                [evi, errcode, EVIvec] = evi_ll(bstmean - mu0vec', (scalevec./shapevec)', bstindx, etavec', mydof', alloc);
                stopcritval = evi - alloc*cvec; % EOC_ck stopping rule
            elseif StopRule == EOC1kStop
                [evi, errcode, EVIvec] = evi_ll(bstmean - mu0vec', (scalevec./shapevec)', bstindx, etavec', mydof', alloc);
                stopcritval = max(EVIvec) - (alloc*cvec);    % EOC_1k stopping rule
            elseif StopRule == KGstarstop
                [stopcritval asdf] = max( predstartest );
                Curbudget = max(1, betastarvec(asdf));
            elseif StopRule == KGonestop
                stopcritval = max( predonetest ) - 1;
            elseif StopRule == GitStop
                if max(okformore) > 0
                    stopcritval=1;
                else
                    stopcritval=-1;
                end
            end
            if stopcritval <= 0 
                % try to find a different budget that works - starting from big
                % budgets and working down
%                        'seeking new alloc'
                notdone = false;
                if (StopRule == EOC1kStop) || (StopRule == EOCckStop)       % if not KG* stopping rule, check other budget values...
                    for j=1:length(totbudgetvec);%:-1:1
                        %check if this budget works
                        Curbudget=totbudgetvec(j);
                        [alloc,mydof] = calcrepslinminunkvar(Curbudget,-mu0vec',(scalevec./shapevec)',etavec',1); % calculate allocation for a big number of samples
                        [evi, errcode, EVIvec] = evi_ll(bstmean - mu0vec', (scalevec./shapevec)', bstindx, etavec', mydof', alloc);
                        if StopRule == EOCckStop   % sets stopping criterion
                            stopcritval = evi - alloc*cvec; % EOC_ck stopping rule
                        elseif StopRule == EOC1kStop
                            stopcritval = max(EVIvec) - (alloc*cvec);    % EOC_1k stopping rule
                        end
                        if stopcritval > 0 % check if we found an allocation that works
                            % found a budget!
                            Curbudget = sum(alloc);
                            notdone = true;
                            break;
                        end
                    end %for
                end%if
            end%if

            % check if stopping rule satisfied
            if notdone
                %take a sample
                [maxrepsval sysindx] = max( alloc );    % pick the alternative that would get the most replications
                if sysindx == k+1
                    alloc
                    'error with allocating to system 1'
                    beep
%                    break;
                end
                Xij = normrnd( truemeanvec(sysindx), sqrt( truevarvec(sysindx)) );

                % update posterior
                scalevec(sysindx) = scalevec(sysindx)+ (Xij - mu0vec(sysindx))^2 * etavec(sysindx) / (etavec(sysindx)+1) / 2;
                mu0vec(sysindx) = mu0vec(sysindx) + (Xij - mu0vec(sysindx))/(etavec(sysindx));
                etavec(sysindx) = etavec(sysindx)+1;
                shapevec(sysindx) = shapevec(sysindx)+1/2;
                
                %(sysindx == bstindx) || (max(mu0vec) == bstmean) || 
                if (sysindx == bstindx) || (max(mu0vec) ~= bstmean)  % if best was simulated or there is a new best...
                    RedoBeta = -1;      % if we simulated the best system, or something tied for best, or have a new best...
                else                % otherwise we only need to redo the beta for the alternative that was simulated
                    RedoBeta = sysindx;
                end
%                RedoBeta = -1;
                % update budget info
                Curbudget = max(1, Curbudget - 1);
            end % check if stopping rule is satisfied.
        end % while notdone
        % record the statistics at the stopping time
        [selexpval selindx] = max(mu0vec);
        selectedvalue(ell) = selexpval;
        oppcostvec(ell) = max(truemeanvec) - truemeanvec(selindx);
        truemeanselected(ell) = max(truemeanvec);
        budgetspent(ell) = (etavec(1:k) - effectivenumreps)' * cvec(1:k);
    end % for each macrorep
    cputim(i) = toc
    meanReward(i) = mean(selectedvalue)
    SEmeanreward(i) = sqrt(var(selectedvalue)/numruns)
    meanOppCost(i) = mean(oppcostvec)
    SEmeanOppCost(i) = sqrt(var(oppcostvec)/numruns)
    truemeanReward(i) = mean(truemeanselected)
    SEtruemeanreward(i) = sqrt(var(truemeanselected)/numruns)
    SEsubopt(i) = sqrt(var((oppcostvec + budgetspent))/numruns);
    varmeanReward(i) = var(selectedvalue - budgetspent)
    vartruemean(i) = var(truemeanReward)
    tmptmp = cov(selectedvalue - budgetspent,truemeanselected);
    covmeanandtrue(i) = tmptmp(1,2)
    meanBudgetUsed(i) = mean(budgetspent)
    SEmeanBudgetUsed(i) = sqrt(var(budgetspent)/numruns)
    OverallNet = OERupperk - meanBudgetUsed - meanOppCost
end % for each value for the number of alternatives
mytitle = sprintf('BigKSimAlloc%sStop%s.mat',int2str(AllocRule),int2str(StopRule))
cputim
kvec
meanReward
SEmeanreward 
meanReward - meanBudgetUsed
truemeanReward
SEmeanBudgetUsed
SEmeanOppCost
sqrt(SEmeanOppCost.^2 + SEmeanBudgetUsed.^2)./OverallNet
OERupperk - OERlowerk
subopt = meanOppCost + meanBudgetUsed
SEsubopt
OverallNet = OERupperk - meanBudgetUsed - meanOppCost
meanBudgetUsed
meanOppCost
sum(cputim)*15000/60/60/NMAXVALS
save(mytitle);


del = 2
lam = sqrt(scale / eta * shape)
sval = del / lam
dof = 2*eta
num= (scale / shape / lam) * (dof + sval^2) * tpdf(sval,dof) / (dof - 1) / 2
denom = lam * Psifunc(sval,dof)
reptest = 0.1:1:2000;
lamtest = sqrt((scale * reptest) ./ (shape * eta * (eta + reptest)));
testtest= lamtest .* Psifunc(del ./ lamtest,dof) - reptest;
plot(reptest,testtest)

OERupperk =
  1.0e+005 *
    0.1758    0.3009    0.3941    0.5212    0.6915    0.8497    1.0432
EOC-ckOverall =
  1.0e+004 *
    1.7253    2.9412    3.8046    4.9176    6.3206    7.5134    9.0484
EOC-1kOverall =
  1.0e+004 *
    1.7255    2.9447    3.7872    4.9394    6.2483    7.5005    8.2652
OERlowerk =
  1.0e+004 *
    1.7136    2.9313    3.8321    5.0571    6.6555    8.0992    9.7392
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

TO HERE???

testvals = macroreps;
    OERupperk(i) = mean(testvals);
    OERupperkSE(i) = std(testvals)/sqrt(length(testvals));
    
    betavec=0*betatst;  % initialize
    for j=1:length(betatst)
        betavec(j) = mean( max(mu0 + (macroreps-mu0)*sqrt(betatst(j)/(eta+betatst(j))), retirement) - c*k*betatst(j));
    end
    [ta tb] = max(betavec);
    OERlowerk(i) = ta;
    OERnruns(i) = betatst(tb)*k;
end


mu0 = 0

shape = 2 + meanofvar^2 / varofvar  % back fit these inverted gamma params using Bernardo & Smith, p. 431
scale = (meanofvar * (shape - 1))   % noting that for matlab, the mean of a gamma is shape*scale...
vecvar = 1./gamrnd(shape*ones(1,Nreps),1/scale);
mean(vecvar)
var(vecvar)
vecmean = normrnd(mu0,sqrt(vecvar / eta));
mean(vecmean)
var(vecmean)
meanofvar/eta       % this and the value of var(vecmean) should be close
maxOER = mean(max(vecmean,retirement))

% The next bit is an alternative way to get student t for unknown mean - as
% a check for prior code. this is not needed, as we also need the variance
altvect = mu0 + sqrt(scale/eta/(shape)) * trnd(2*shape,1,Nreps);    % again, shape*scale instead of scale/shape
%tst1 = mu0 + trnd(2*shape,1,Nreps);    % again, shape*scale instead of scale/shape
mean(altvect)
var(altvect)
altmaxOER = mean(max(altvect,retirement))

'tests with k = 1'
Nreps = 5;
vecvar = 1./gamrnd(shape*ones(1,Nreps),1/scale);
vecmean = normrnd(mu0,sqrt(vecvar / eta));

MAXREPS = 5000 % maximum number of samples to take
tmprepvec = ones(MAXREPS,1);
nsamples = normrnd(tmprepvec*vecmean,tmprepvec*sqrt(vecvar));   %matrix: is MAXREPS x Nreps

tvec = eta + (1:MAXREPS);
numsamps = (1:MAXREPS)';    %row vector
sumx = cumsum(nsamples);        % put
posteriormean = (mu0*eta + sumx) ./ ((eta + numsamps)*ones(1,Nreps));
postmean = zeros(size(nsamples));
postvar = postmean;

hourcost = 1; %10;      
sigma = sqrt(meanofvar); %10^6;     % std dev during one time unit
%sig0 = sigma / sqrt(t0) ; % for prior distr of unknown means
cpumin = 60;       % cpu time to get output of that unit
%m = 0
[c gamma beta alpha] = CFBranchSetParams(hourcost, cpumin, sigma);

cfactor = c^(1/3)
for i=1:Nreps
       % Initialize parameters for inference
    myscale = scale;
    myshape = shape;
    myeta = eta;
    mymu0 = mu0;
    samples(i) = 0;
    trueboundapprox = ApproxBoundW(1/gamma./tvec)/beta;
    boundapprox = trueboundapprox;
        % do the inference until maximum number of samples observed, or
        % until outside stopping set
    inContinuation = true;
    for j=1:MAXREPS
        % observe next sample for run
        myeta = myeta + 1;
        myshape = myshape + 1/2;
        myscale = myscale + (nsamples(j,i) - mymu0)^2 * myeta / (myeta+1) / 2;
        mymu0 = mymu0 + (nsamples(j,i) - mymu0)/(myeta);
        postmean(j,i) = mymu0;
        postvar(j,i) = myscale / (myshape - 1);
        sval = postvar(j,i)^(1/3) / cfactor^2 / myeta;
        bndapprox = cfactor * postvar(j,i)^(1/3) * ApproxBoundW(sval);
        boundapprox(j) = bndapprox;
    end
%    svec = postvar(:,i).^(1/3) ./ cfactor^2 ./ tvec';
%    boundapprox = cfactor * postvar(:,i).^(1/3) .* ApproxBoundW(svec);
    figure(PLOTBASE+i)
    plot(tvec,boundapprox,'--',tvec,-boundapprox,'--',tvec,postmean(:,i),'-',tvec,trueboundapprox,'-.')
    xlabel('Effective number of samples, \eta','FontSize',myfontsize,'FontName','Times'); 
    ylabel('Posterior mean of unknown mean','FontSize',myfontsize,'FontName','Times')
    title('Approximation of continuation set with a sample path of posterior mean','FontSize',myfontsize,'FontName','Times')
%    pause

    figure(PLOTBASE+Nreps+i)
    EstimatedBoundDuringRun = cfactor * postvar(:,i).^(1/3) .* ApproxBoundW(postvar(:,i).^(1/3) / cfactor^2 ./ tvec');
    truevar=vecvar(i)
    boundwithvarknown = cfactor * truevar^(1/3) * ApproxBoundW(truevar^(1/3) / cfactor^2 ./ tvec');
    meanofvar
    boundwithmeanvar = cfactor * meanofvar^(1/3) * ApproxBoundW(meanofvar^(1/3) / cfactor^2 ./ tvec');
    loglog(tvec,EstimatedBoundDuringRun,'-',tvec,boundwithvarknown,'--',tvec,boundwithmeanvar,'-.')
    legend('Estimated Boundary','Bound if Realized Variance Were Known','Bound With Mean Variance')
    xlabel('Effective number of samples, \eta','FontSize',myfontsize,'FontName','Times'); 
    ylabel('Posterior mean of unknown mean','FontSize',myfontsize,'FontName','Times');

    pause

    mytitle = strcat(strcat('CompareContinSets',int2str(i)),'.eps');
    print('-deps',mytitle);	
    figure(PLOTBASE+i)
    mytitle = strcat(strcat('UnkVarSampBound',int2str(i)),'.eps');
	print('-deps',mytitle);	
end

figure(PLOTBASE+1+2*Nreps)
plot(tvec,postvar(:,i))
xlabel('Effective number of samples, \eta','FontSize',myfontsize,'FontName','Times'); 
ylabel('Posterior mean of unknown variance','FontSize',myfontsize,'FontName','Times');
pause
mytitle = strcat(strcat('UnkVarSampVar',int2str(i)),'.eps');
print('-deps',mytitle);	

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

oldstuff

% PARAMETERS FOR STOPPING RULE                
STOP_NOBOUND = 0; % no stopping boundary, only a time-limiting boundary
STOP_PDE = 1;     % stopping boundary from PDE solution
STOP_BOUND = 2;   % stopping boundary from lower bound on index
STOP_RULE = STOP_PDE
MAXDAYS = 24.5   % for all stopping rules, need an overall total limit for number of reps to prevent infinite loop

% PARAMETERS FOR ALLOCATION RULE
ALLOC_EQUAL = 0;         % allocation replications equally to systems in round robin fashion
ALLOC_PDE = 1;           % Use PDE approximation to determine allocation index
ALLOC_BOUND = 2;         % Use bound for Gittins index to determine allocation index
ALLOC_EOC = 3;            % Use EOC rule
ALLOC_RULE = ALLOC_PDE

STOP_LOWCHECK = ((STOP_RULE == 1) + (STOP_RULE == 2)) * ((ALLOC_RULE == 1) + (ALLOC_RULE == 2));

WIPEOUTSTART = 0
if WIPEOUTSTART
    INITREPS = round(n0base);       % number of reps to take from each system, no matter what
else
    INITREPS = n0base;
end

basename = 'NoWipeGetRepsD01n04';
MAXREPS = round(MAXDAYS*24*60/cpumin) %;   %40000;
MAXREPS*gamma;

maxk = 5;
nmacroreps = 1600;
FRACSTOP = MYEPS;    % fraction of relative error for Gittins index versus value of stopping, ...
                    % if stopping is to occur, using the extra bounds
% Allocation For 
totalreps=zeros(maxk,nmacroreps);
lowgitreps=zeros(maxk,nmacroreps);
truemeanreward=zeros(maxk,nmacroreps);
gotbest=zeros(maxk,nmacroreps);
eocloss=zeros(maxk,nmacroreps);
discmeanreward=zeros(maxk,nmacroreps);
mingit=zeros(maxk,nmacroreps);
eval=zeros(1,maxk);
disceval=zeros(1,maxk);
truedisceval=zeros(1,maxk);
discse=zeros(1,maxk);
truese=zeros(1,maxk);
stopped=zeros(1,maxk);

tlist = n0base + (1:MAXREPS);
% THIS ASSUMES THAT KAPPAINVVEC and BETA ARE SAME FOR ALL SYSTEMS!
stopvec = (BranchBndFromS(1/gamma./tlist,accumsvec,accumbound))/beta; % find stopping boundary as function of number of replications
stop0 = (BranchBndFromS(1/gamma/n0base,accumsvec,accumbound)) / beta

fileloaded
max(0,BranchGetB1Val(0,1/gamma/n0base))/beta - c/gamma
[ta tb tc]=fminbnd('Bztauapproxc',0,1/gamma/n0base,optimset('TolFun',beta/sigma,'MaxIter',10^4),0,1/gamma/n0base);       % FIXED VERSION
-tb/beta - c/gamma
sig0*PsiNorm(-mu0/sig0)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%for k=1:maxk;                     % number of systems 
%for k=maxk:-1:1;                     % number of systems 
for k=maxk:-1:4;                     % number of systems 
    mu0vec = mu0 * ones(1,k);        % set up prior distribution for unknown means
    n0vec = n0base * ones(1,k);      
    sigvec = sigma * ones(1,k); % standard deviations of output for each system
%    cvec = c * ones(1,k);       %c;    % cost, by system
    if PLOTON
        figure(PLOTBASE+k); hold off
        plot(1./accumsvec/gamma,accumbound/beta);
        tmp=axis; tmp(1) = n0base / 2; tmp(2) = (10*24*60/cpumin) + tmp(1); tmp(3) = - sig0*4; tmp(4) = sig0*4;
        axis(tmp); hold on
    end
    
    for i=1:nmacroreps                              % SAMPLE A NEW PROBLEM STANTEMENT
        truemeanvec = mu0vec + randn(1,k) .* sig0;	% generate a random configuration from prior distribution
        sumoutput = INITREPS * truemeanvec + sqrt(INITREPS) * randn(1,k) .* sigvec;     % for statistics collection
        if WIPEOUTSTART
            nvec = zeros(1,k);          % number of observations from systems 1 ... k
        else
            nvec = INITREPS * ones(1,k);          % number of observations from systems 1 ... k
        end
        tvec = n0vec + nvec;        % initialize the statistics for the prior/posterier distributions
        bestybyt = -c/gamma;        % best of values for 'stopped' systems
        bestsys = 0;

        ybyt = (n0vec .* mu0vec + sumoutput) ./ tvec;
%        wwigglevec = (ybyt) * beta + kappainvvec;                     % convert initial statistics to standardized form
        wwigglevec = (ybyt - bestybyt) * beta;                     % convert initial statistics to standardized form
        swigglevec = 1 / gamma ./ tvec;
        stopped=zeros(1,k);          % toggles as to whether a given system has been 'selected'
        
        lastmean = ybyt;
        lastt = tvec;
        
% initialize the Gittins indices
%        Gitvec=max(0,BranchGetB1Val(wwigglevec,swigglevec,accumsvec,accumbound,fileindx))/beta - c/gamma;
%TEST OF USING LOWER BOUND ON GITTINS INDEX AS AN INDEX, RATHER THAN PDE
%SOLUTION
        Gitvec=zeros(size(wwigglevec));
        if ALLOC_RULE == ALLOC_PDE
            Gitvec=max(0,BranchGetB1Val(wwigglevec,swigglevec))/beta + bestybyt;
            if (i==1)*(k==maxk)
                fileloaded
            end
        else   %if ALLOC_RULE == ALLOC_BOUND
            Gitpos=zeros(size(wwigglevec));
%            fraccheck=0.01;
            for ell=1:k
%                [ta tb tc]=fminsearch('Bztauapproxc',fraccheck/gamma./n0vec(ell),optimset('TolFun',beta/sigvec(ell),'MaxIter',10^4),kappainvvec + beta*mu0vec(ell), 1/gamma./n0vec(ell));       % FIXED VERSION
%                [ta tb tc]=fminsearch('Bztauapproxc',fraccheck*swigglevec(ell),optimset('TolFun',beta/sigvec(ell),'MaxIter',10^4),wwigglevec(ell), swigglevec(ell));       % FIXED VERSION
                [ta tb tc]=fminbnd('Bztauapproxc',0,swigglevec(ell),optimset('TolFun',beta/sigvec(ell),'MaxIter',10^4),wwigglevec(ell), swigglevec(ell));       % FIXED VERSION
                Gitvec(ell) = -tb/beta + bestybyt;
                Gitpos(ell) = ta;
            end
        end
%        tmpgit = Gitvec - max(mu0vec+c/gamma,0);
%        Gitvec = 2 * tmpgit + max(mu0vec+c/gamma,0);           % This is a kludge, while waiting for better estimates
        Gitstart(k) = Gitvec(1);
        notdone = 1;
        mingit(k,i) = max(Gitvec);

        while notdone                               % apply the simulation selection algorithm
            % COMPUTE ALLOCATION
            if ALLOC_RULE == ALLOC_EQUAL
                [Gitval Indxval] = max(Gitvec);         % get index and Gittins index for system with biggest Gittins index
                [junk Indxval] = min(tvec);         % get index and number of reps for system with fewest observations
            elseif ALLOC_RULE == ALLOC_EOC
                [Gitval Indxval] = max(Gitvec);         % get index and Gittins index for system with biggest Gittins index
                [alloc] = calcrepslinminmeanvar(1,-ybyt,sigvec.^2,tvec);
                [junk Indxval] = max(alloc);
            else
                [Gitval Indxval] = max(Gitvec);         % get index and Gittins index for system with biggest Gittins index
            end
            mingit(k,i) = min(mingit(k,i),Gitval);  % record the minimum of (the best of the Gittins indices) so far

            if nvec(Indxval) > 0
                stopval = stopvec(round(nvec(Indxval))) + bestybyt;
            else
                stopval = stop0 + bestybyt;
            end
            totreps = sum(nvec);

            % CHECK STOPPING CONDITION
            if (STOP_RULE==STOP_PDE)*(ybyt(Indxval) >= stopval) + (STOP_RULE==STOP_BOUND)*(ybyt(Indxval) >= stopval)*( abs(Gitvec(Indxval) - ybyt(Indxval)) <= (FRACSTOP+MYEPS) * Gitvec(Indxval)) % return
%               'Hit stopping boundary', use first case if PDE is used,
%               use second case if Gittins index approximation boudnary is
%               used.
%                if nvec(Indxval) > 1           % add in a linear correction term
%                    correction(k,i) = 0;% - (ybyt(Indxval) - lastmean(Indxval)) * ( ybyt(Indxval) - stopval ) / ( (ybyt(Indxval) - stopval) - (lastmean(Indxval) - stopvec(nvec(Indxval)-1)) ) ;
%                else 
%                    correction(k,i) = 0;% - (ybyt(Indxval) - mu0vec(Indxval)) * ( ybyt(Indxval) - stopval ) / ( (ybyt(Indxval) - stopval) - (mu0vec(Indxval) - stop0) );
%                end
%                stopped=0*stopped;
                stopped(Indxval)=1;
                if (prod(stopped)==1)           % all systems are now stopped       FIX THE FOLLOWING FINAL REWARDS
                    'all stopped'
                    [obsval Indxval] = max(ybyt);
                    truemeanreward(k,i) = exp(-gamma * totreps) * (truemeanvec(Indxval)) - (1-exp(-gamma*totreps))*c/gamma;
                    discmeanreward(k,i) = exp(-gamma * totreps) * ybyt(Indxval) - (1-exp(-gamma*totreps))*c/gamma;
                    diffvec(k,i) = max(0,ybyt(Indxval) - stopval);
                    gotbest(k,i) = (truemeanvec(Indxval)==max(truemeanvec));
                    eocloss(k,i) = max(truemeanvec) - truemeanvec(Indxval);
                    totalreps(k,i) = totreps;
                    if (i==1)*(k==maxk)
                        fileloaded
                    end
                    notdone = 0;                      % break out of the loop
                else
%                    Gitvec(Indxval) = - c / gamma;      % make it so that this system will not be selected again
                    bestybyt = max(bestybyt,ybyt(Indxval));           % the 'best' is now improved to ybyt
%                    'removing ',Indxval
                    wwigglevec = (1-stopped).*(ybyt - bestybyt) * beta + stopped*(-2*sigma)*beta;                     % convert initial statistics to standardized form
                    if ALLOC_RULE == ALLOC_PDE
                        Gitvec=max(0,BranchGetB1Val(wwigglevec,swigglevec))/beta + bestybyt;
                        for ell=1:k
                            if stopped(ell)
                                Gitvec(ell)=Gitvec(ell)-1;
                            end
                        end
                    else   %if ALLOC_RULE == ALLOC_BOUND
                        Gitpos=zeros(size(wwigglevec));
                        for ell=1:k
                            if stopped(ell)
                                Gitvec(ell) = - c/gamma;
                            else
                                [ta tb tc]=fminbnd('Bztauapproxc',0,swigglevec(ell),optimset('TolFun',beta/sigvec(ell),'MaxIter',10^4),wwigglevec(ell), swigglevec(ell));       % FIXED VERSION
                                Gitvec(ell) = -tb/beta + bestybyt;
                                Gitpos(ell) = ta;
                            end
                        end
                    end
                end
            elseif totreps >= MAXREPS               % we have done the maximum permitted number of replications without hitting
%                'Maximum number of replications'
                [obsval Indxval] = max(ybyt);
                if (i==1)*(k==maxk)
                    fileloaded
                end
%                correction(k,i) = 0;                % the stopping boundary, so let's stop and take default action
                truemeanreward(k,i) = exp(-gamma*totreps) * max(-c/gamma,truemeanvec(Indxval)) - (1-exp(-gamma*totreps))*c/gamma;
                discmeanreward(k,i) = exp(-gamma*totreps) * max(-c/gamma,obsval) - (1-exp(-gamma*totreps))*c/gamma;
                diffvec(k,i) = max(0,obsval - stopval);
                [aw bw] = max(truemeanvec);
                gotbest(k,i) = (Indxval == bw);
                eocloss(k,i) = max(truemeanvec) - truemeanvec(Indxval);
                totalreps(k,i) = totreps;
                notdone = 0;                      % break out of the loop
            elseif STOP_LOWCHECK * ((Gitval <= (bestybyt+MYEPS)) + ( (ybyt(Indxval)-bestybyt)/(sigma/sqrt(tvec(Indxval))) < -3.5 ))                     % here, the maximum of the Gittin's indices evaluates to 0, and we haven't exhausted the computer budget
%                'Zero option is best, or very poor outlook'
%                correction(k,i) = 0;                % the stopping boundary, so let's stop and take default action
                [obsval Indxval] = max(ybyt);
                truemeanreward(k,i) = exp(-gamma*totreps) * max(-c/gamma,truemeanvec(Indxval)) - (1-exp(-gamma*totreps))*c/gamma;
                discmeanreward(k,i) = exp(-gamma*totreps) * max(-c/gamma,obsval) - (1-exp(-gamma*totreps))*c/gamma;
                diffvec(k,i) = max(0,obsval - stopval);
                [aw bw] = max(truemeanvec);
                gotbest(k,i) = (Indxval == bw) | ((obsval < 0) & (aw < 0));
                eocloss(k,i) = max(truemeanvec) - truemeanvec(Indxval);
                totalreps(k,i) = totreps;
                notdone = 0;                      % break out of the loop
            else        % run a replication
                lastt(Indxval) = tvec(Indxval);         % save last value of mean etc.
                lastmean(Indxval) = ybyt(Indxval);
                simrep = truemeanvec(Indxval) + randn * sigvec(Indxval);         % run a replication
                sumoutput(Indxval) = simrep + sumoutput(Indxval);
                nvec(Indxval) = nvec(Indxval) + 1;
                tvec(Indxval) = n0vec(Indxval) + nvec(Indxval);
                swigglevec(Indxval) = 1 / gamma / tvec(Indxval);
                ybyt(Indxval) = (n0vec(Indxval) * mu0vec(Indxval) + sumoutput(Indxval)) / tvec(Indxval);
%               wwigglevec(Indxval) = ybyt(Indxval) * beta + kappainvvec;                     % convert initial statistics to standardized form
                wwigglevec(Indxval) = (ybyt(Indxval)-bestybyt) * beta;                     % convert initial statistics to standardized form
                if ALLOC_RULE == ALLOC_PDE
                    Gitvec(Indxval) = max(0,BranchGetB1Val(wwigglevec(Indxval),swigglevec(Indxval)))/beta + bestybyt;
                elseif (ALLOC_RULE == ALLOC_BOUND)+(STOP_RULE == STOP_BOUND)  %TEST OF USING LOWER BOUND ON GITTINS INDEX AS AN INDEX, RATHER THAN PDE SOLUTION
%                    [ta tb tc]=fminsearch('Bztauapproxc',Gitpos(Indxval),optimset('TolFun',beta/sigvec(Indxval),'MaxIter',10^3),kappainvvec + beta*ybyt(Indxval), 1/gamma./tvec(Indxval));       % FIXED VERSION
                    [ta tb tc]=fminbnd('Bztauapproxc',0,swigglevec(ell),optimset('TolFun',beta/sigvec(ell),'MaxIter',10^4),wwigglevec(Indxval), swigglevec(Indxval));       % FIXED VERSION
                    Gitvec(Indxval) = -tb/beta + bestybyt;
                    if ta<= 0
                        Gitpos(Indxval) = 1/gamma/tvec(Indxval)/2;
                    else
                        Gitpos(Indxval) = ta;
                    end 
                end
 
                if PLOTON*(k==1)*(i<=20)
                    plot([lastt(Indxval) tvec(Indxval)], [lastmean(Indxval) ybyt(Indxval)],'LineWidth',Indxval)
                end
            end
        end
        if mod(i,100)==0 
            [ k i nvec ]
            ybyt
            truemeanvec
            stopped
            bestybyt
        end
    end

    truedisceval(k) = mean(truemeanreward(k,:));
    truese(k) = std(truemeanreward(k,:))/sqrt(nmacroreps);
    disceval(k) = mean(discmeanreward(k,:));
    discse(k) = std(discmeanreward(k,:))/sqrt(nmacroreps);

avediscmean=mean(discmeanreward') 
avedisctrue=mean(truemeanreward')
pcsvector=mean(gotbest')
eocvector=mean(eocloss')
fracvarious=[sum(totalreps'==MAXREPS) / nmacroreps  sum(discmeanreward'<=-c/gamma)/nmacroreps]
meandiffcorr=[mean(diffvec')]
avenumruns=mean(totalreps');
avedays = avenumruns * cpumin / 24 / 60

    figure(k)
    [f, xf] = ecdf(discmeanreward(k,:));
    stairs(xf,f);
    ylabel('Empirical CDF','FontSize',20,'FontName','Times')
    xlabel('E[NPV]','FontSize',20,'FontName','Times')
	mytitle = strcat(strcat('EmpDistrENPV',basename,int2str(k)),'.eps');
	print('-deps',mytitle);
	h=gcf;
	mytitle = strcat(strcat('EmpDistrENPV',basename,int2str(k)),'.fig'); 
    saveas(h,mytitle,'fig');

    if (ALLOC_RULE == ALLOC_PDE)+(ALLOC_RULE == ALLOC_BOUND)
        figure(maxk+k)
        [f, xf] = ecdf(mingit(k,:));
        stairs(xf,f);
        ylabel('Empirical CDF','FontSize',20,'FontName','Times')
        xlabel('Minimum Value of Gittins Index During Run','FontSize',20,'FontName','Times')
        mytitle = strcat(strcat('mingit',basename,int2str(k)),'.eps');
        print('-deps',mytitle);
        h=gcf;
        mytitle = strcat(strcat('mingit',basename,int2str(k)),'.fig'); 
        saveas(h,mytitle,'fig');
    end
    
    figure(maxk+1)
    hold on
    stairs(xf,f.^(1/k),'-.');
    hold off

end % for k = 1:maxk

Gitstart
avediscmean=mean(discmeanreward') 
avedisctrue=mean(truemeanreward')
pcsvector=mean(gotbest')
eocvector=mean(eocloss')
fracvarious=[sum(totalreps'==MAXREPS) / nmacroreps  sum(discmeanreward'<=-c/gamma)/nmacroreps]
meandiffcorr=[mean(diffvec')]
relerror=std(discmeanreward')/sqrt(nmacroreps/2) ./  mean(discmeanreward')
avenumruns=mean(totalreps')
avedays = avenumruns * cpumin / 24 / 60

sum(totalreps'>=MAXREPS)

figure(2*maxk+1)
plot(mean(totalreps')*cpumin/60/24)
    ylabel('Elapsed CPU Days','FontSize',20,'FontName','Times')
    xlabel('Number of Systems','FontSize',20,'FontName','Times')
	mytitle = strcat(strcat('EmpDistrCPU',basename),'.eps')
	print('-deps',mytitle);
	h=gcf;
	mytitle = strcat(strcat('EmpDistrCPUD5',basename),'.fig'); 
    saveas(h,mytitle,'fig');
figure(2*maxk+2)
plot(disceval)
    ylabel('E[E[NPV]]','FontSize',20,'FontName','Times')
    xlabel('Number of Systems','FontSize',20,'FontName','Times')
	mytitle = strcat(strcat('ENPVD5',basename),'.eps')
	print('-deps',mytitle);
	h=gcf;
	mytitle = strcat(strcat('ENPVD5',basename),'.fig'); 
    saveas(h,mytitle,'fig');

