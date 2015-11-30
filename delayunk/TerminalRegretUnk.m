function [ereward, pcs] = TerminalRegretUnk(valuevec,discountfactor,predvarsample,postvar,advanced,dof)
% (c) 2014 Chick, Forster, Pertile
% This file is for input to matlab, and does calculations
% to support the paper of Steve, Paolo and Martin on delayed response in two-arm clinical trial.
% 
% terminal regret with case of unknown variance, does not use mills ratio.
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
nochangeallowed = advanced.nochangeallowed;
DOINTEGRAL = advanced.DoRegretIntegral;
NUMSTEPS = advanced.NumPointsQuadrature;
if nargin < 6       % If last argument is missing, assume that no change is allowed.
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
    ereward = discountfactor^(1-nochangeallowed) * predstd * PsiNormUV(testvec, dof);  % get expected reward of decision later: this is used to back out the regret later, by subtracting it from value of perfect info

    if DOINTEGRAL  % Do numerical integration in area where quadrature is potentially a poor approximation
        ZLIMHI = 19;
        ZLIMLO = 20; % SEC: in debug effor in Sept 2014, tried setting this to 0. 

        tmpvec=((abs(valuevec)/predstd)<=ZLIMHI)&((abs(valuevec)/predstd)>=ZLIMLO);

        pcs=zeros(size(valuevec));

        if sum(tmpvec) % in area where potentially poor approximation, do numerical integration
            pcs(tmpvec) = integral(@(x) tcdf(abs(x+valuevec(tmpvec))/poststd, dof)*tpdf(x/predstd, dof), ... FIX
                -4*predstd, +4*predstd,'ArrayValued',true,'RelTol',1e-2);
        end
        
        % do quadrature approx outside of that range, as it is much faster
        tmp=1/NUMSTEPS; % the tau samples arrive, and then over the realized mean given the posterior after those tau samples arrive.
        pvec=-tmp/2+(1:NUMSTEPS)*tmp;   % to do this, find a set of zvalues to do the integration, rather than calling a loop with 
        zvec=tinv(pvec,dof);         % integrations over values of each element of valuevec.  hopefully matrix mult will be faster
        tmpvec = ~tmpvec;

        if sum(tmpvec)
            mulen = sum(tmpvec);
            tmpmu = valuevec(tmpvec);
            tmpmu = tmpmu(:)*ones(1,NUMSTEPS);
            tmpdel = ones(mulen,1)*zvec;
            bigmatrix=tmpmu'+ predstd*tmpdel';
            zmatrix = abs(bigmatrix) / poststd;

            Pmatrix = tcdf ( zmatrix, dof );
            pcs(tmpvec) = mean(Pmatrix);
        end
    else
        tmp=1/NUMSTEPS; % the tau samples arrive, and then over the realized mean given the posterior after those tau samples arrive.
        pvec=-tmp/2+(1:NUMSTEPS)*tmp;   % to do this, find a set of zvalues to do the integration, rather than calling a loop with 
        zvec=tinv(pvec,dof);         % integrations over values of each element of valuevec.  hopefully matrix mult will be faster
        
        mulen = length(valuevec);
        tmpmu = valuevec(:)*ones(1,NUMSTEPS);
        tmpdel = ones(mulen,1)*zvec;
        bigmatrix=tmpmu'+ predstd*tmpdel';
        zmatrix = abs(bigmatrix) / poststd;
        Pmatrix = tcdf( zmatrix, dof );
        pcs = mean(Pmatrix);
    end
elseif predvarsample == 0       % this line indicates that if no change is allowed, the decision can be implemented without waiting for the delay of tau patients
    if postvar > 0
        testvec = -(valuevec) / sqrt(postvar);
        pcs = tcdf(abs(testvec),dof);
        ereward = (discountfactor^(1-nochangeallowed)) * sqrt(postvar) * PsiNormUV(testvec, dof, pcs); 
    else
        ereward = 0*valuevec;  %max(valuevec,0);
        pcs = 0*valuevec;
    end
else    % if value is negative, treat it as positive and use an approximation. Here, we go for regret not expected loss, so we use absolute values here
    testvec = -(valuevec) / sqrt(-predvarsample);
    pcs = tcdf(abs(testvec), dof);
    ereward = discountfactor^(1-nochangeallowed) * sqrt(-predvarsample) * PsiNormUV(testvec, dof, pcs);
	warning('TerminalRegretUnk: called with negative predvarsample');
end
