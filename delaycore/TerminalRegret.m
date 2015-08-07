function [ereward, pcs] = TerminalRegret(valuevec,discountfactor,predvarsample,postvar,advanced)
% (c) 2014 Chick, Forster, Pertile
% This file is for input to matlab, and does calculations
% to support the paper of Steve, Paolo and Martin on delayed response in two-arm clinical trial.
% Edited 2014 to include improvement in computation value from Mills ratio
%
% Code is 'as is' and no guarantees for correctness.
%
% INPUTS:
% valuevec : vector of means of posterior distribution
% discountfactor: discountfactor to apply to reward after the patients are seen
% predvarsample: a scalar with the predictive variance of the posterior mean
% postvar: the posterior variance remaining at the time the samples are all in
% advanced: if included, then nochangeallowed and NumPointsQuadrature are used from
%   that variable in order to compute the various integrals, etc. If
%   advanced is not passed, then nochangeallowed is assumed to be false by
%   default.
%
% Notes: 
% a) if prevarsample is 0, then the regret is 0 (no value in information)
% b) if predvarsample is set as negative, it is interpreted as the negative
% of the actual predictive variance of the unknown mean, and the realized
% mean is taken (as an approximation) to be the posterior mean. This makes
% code run faster, but is an approximation.
%
% OUTPUTS:
% ereward : vector of outputs with expected reward
% pcs : vector of outputs with expected posterior probability of picking best

% validate inputs
if nargin < 5       % If last argument is missing, assume that no change is allowed.
    [~, advanced] = DelayInputConstructor();
end
nochangeallowed = advanced.nochangeallowed;
DOINTEGRAL = advanced.DoRegretIntegral;
NUMSTEPS = advanced.NumPointsQuadrature;
if nargin < 4
  'TerminalRegret: not enough arguments'
end

% process inputs
if nochangeallowed % if no change is allowed, then it implies that there is no variance in posterior decision, or pred variance is 0
    predvarsample = 0;
end
if predvarsample > 0
    poststd = sqrt(postvar);
    predstd = sqrt(predvarsample);

    testvec = -(valuevec) / predstd;
    ereward = discountfactor^(1-nochangeallowed) * predstd * PsiNorm(testvec);  % get expected reward of decision later: this is used to back out the regret later, by subtracting it from value of perfect info

    if DOINTEGRAL  % Do numerical integration in area where quadrature is potentially a poor approximation
        ZLIMHI = 4.75;
        ZLIMLO = 1.4; % SEC: in debug effor in Sept 2014, tried setting this to 0. 

        tmpvec=((abs(valuevec)/predstd)<=ZLIMHI)&((abs(valuevec)/predstd)>=ZLIMLO);
%        eoc=zeros(size(valuevec));

        pcs=zeros(size(valuevec));

        if sum(tmpvec) % in area where potentially poor approximation, do numerical integration
%            eoc(tmpvec) = integral(@(x) poststd*PsiNorm(abs(x+valuevec(tmpvec))/poststd)*normpdf(x,0,predstd), ...
%                -4*predstd, +4*predstd,'ArrayValued',true,'RelTol',1e-2);
            pcs(tmpvec) = integral(@(x) normcdf(abs(x+valuevec(tmpvec))/poststd)*normpdf(x,0,predstd), ...
                -4*predstd, +4*predstd,'ArrayValued',true,'RelTol',1e-2);
        end
        
        % do quadrature approx outside of that range, as it is much faster
        tmp=1/NUMSTEPS; % the tau samples arrive, and then over the realized mean given the posterior after those tau samples arrive.
        pvec=-tmp/2+(1:NUMSTEPS)*tmp;   % to do this, find a set of zvalues to do the integration, rather than calling a loop with 
        zvec=norminv(pvec,0,1);         % integrations over values of each element of valuevec.  hopefully matrix mult will be faster
        tmpvec = ~tmpvec;

        if sum(tmpvec)
            mulen = sum(tmpvec);
            tmpmu = valuevec(tmpvec);
            tmpmu = tmpmu(:)*ones(1,NUMSTEPS);
            tmpdel = ones(mulen,1)*zvec;
            bigmatrix=tmpmu'+ predstd*tmpdel';
            zmatrix = abs(bigmatrix) / poststd;

            Pmatrix = normcdf ( zmatrix );
            pcs(tmpvec) = mean(Pmatrix);
%            Lmatrix = poststd * PsiNorm(zmatrix,Pmatrix);
%            eoc(tmpvec) = mean(Lmatrix);    
        end
%        eoc = discountfactor^(1-nochangeallowed) * integral(@(x) poststd*PsiNorm(abs(x+valuevec)/poststd)*normpdf(x,0,predstd), ...
%            -4*predstd, +4*predstd,'ArrayValued',true,'AbsTol',2e-3);
%        pcs = integral(@(x) normcdf(abs(x+valuevec)/poststd)*normpdf(x,0,predstd), ...
%            -4*predstd, +4*predstd,'ArrayValued',true,'AbsTol',1e-3);
    else
        tmp=1/NUMSTEPS; % the tau samples arrive, and then over the realized mean given the posterior after those tau samples arrive.
        pvec=-tmp/2+(1:NUMSTEPS)*tmp;   % to do this, find a set of zvalues to do the integration, rather than calling a loop with 
        zvec=norminv(pvec,0,1);         % integrations over values of each element of valuevec.  hopefully matrix mult will be faster
        
        mulen = length(valuevec);
        tmpmu = valuevec(:)*ones(1,NUMSTEPS);
        tmpdel = ones(mulen,1)*zvec;
        bigmatrix=tmpmu'+ predstd*tmpdel';
        zmatrix = abs(bigmatrix) / poststd;
        Pmatrix = normcdf( zmatrix );
        pcs = mean(Pmatrix);
 %       Lmatrix = poststd * PsiNorm(zmatrix, Pmatrix);
 %       eoc = discountfactor^(1-nochangeallowed) * mean(Lmatrix);
    end
elseif predvarsample == 0       % this line indicates that if no change is allowed, the decision can be implemented without waiting for the delay of tau patients
    if postvar > 0
        testvec = -(valuevec) / sqrt(postvar);
        pcs = normcdf(abs(testvec));
        ereward = (discountfactor^(1-nochangeallowed)) * sqrt(postvar) * PsiNorm(testvec,pcs); 
    else
        ereward = 0*valuevec;  %max(valuevec,0);
        pcs = 0*valuevec;
    end
else    % if value is negative, treat it as positive and use an approximation. Here, we go for regret not expected loss, so we use absolute values here
    testvec = -(valuevec) / sqrt(-predvarsample);
    pcs = normcdf(abs(testvec));
    ereward = discountfactor^(1-nochangeallowed) * sqrt(-predvarsample) * PsiNorm(testvec,pcs);
	warning('TerminalRegret: called with negative predvarsample');
end
