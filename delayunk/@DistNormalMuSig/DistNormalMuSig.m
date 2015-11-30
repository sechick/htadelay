classdef DistNormalMuSig < handle
% DistNormalMuSig:
% Object to represent random variables which have normal distribution with
% unknown mean Mu and unknown standard deviation Sigma.
% 
% public properties:
%   Hyperparameters: the following are the hyperparameters which define the 
%       prior distribution for the unknown parameters and/or other fixed 
%       for the sampling distribution. 
%       For this object, these are:
%           mu0, t0, xi0, chi0
%       so the unknown mean W | sigma ~ Normal(mu0, sigma^2/t0), 
%       variance is Sigma ~ InvGamma(xi0, chi0),
%       and samples have distribution X | w, sigma ~ Normal(w, sigma^2)
%       (mean of Sigma is chi0/(xi0-1))
%
%   Thetavec: (read only) matrix of parameter values sampled from prior distribution
%      defined by the hyperparameters. One column per sampled parameter, and
%      one row per dimension in the parameter (if dim = 1 then matrix is a row vector).
%      Can be the empty matrix []. For normal with unknown mean and
%      variance, the first row is the set of unknown means, and the second
%      row is the unknown variances.
%   Thetacdf: (read only) cumulative distribution function of the sampled theta in
%       Thetavec, or -1 if not available. row vector.
%   Samplecdfmat: (read only) cumulative distribution function of samples,
%       matrix, one col for each col of Thetacdf, one row per sample. Can
%       be empty matrix if no samples have been taken.
%
% read-only properties
%   ESample: (realized) mean of sampling distribution
%   EVariance: (mean) variance of sampling distribution, marginal
%   over parameter uncertainty
%
% public methods: (beyond get/set)
%   Constructor: DistNormalMu(mu0, t0, xi0, chi0, defaultmu)
%       
%   SetMean: change the hyperprior parameters so that the mean value of the
%       sampled output is the given value (marginal distribution of output
%       given hyperparameter). Adjust parameters other than mean if needed.
%       
%   Sample: generate a sample from the given parameter. requires Thetavec
%      to be instantiated by one of the following two routines.
%
%	SampleThetaBayes: % sample a parameter from the prior distribution,
%   	given the hyperparameter (univariate or vector). samples in rows,
%   	different fields of unknown parameter in columns
%
%	SampleThetaFreq:  % default parameter value, given hyperparameter 
%       (univariate or vector). samples in rows, different fields of 
%       unknown parameter in columns
%
%   TestCase: % calls the routines above, might be useful to see
%       functionality and for use in validating proper operation of the
%       code.
%
    properties
        mu0
        t0
        xi0
        chi0
    end    
    properties (Dependent)
        ESample        % returns the mean of a sample, given the
        EVariance    % returns the mean of a sample, given the
    end
    properties (SetAccess=private,GetAccess=public) % read-only properties
        Thetavec        % matrix of sampled parameter values, one row per sampled parameters, one column per unknown parameter dimension.
        % for the case of normal distribution with unknown mean, one row
        % per sample, and one column with sampled (unknown) mean of
        % sampling distribution
        Thetacdfvec     % cdfs of the sampled Thetavec (if Bayes sampled and cdf is known, otherwise contains -1;
        Samplecdfmat     % cdfs of the sampled Thetavec (if Bayes sampled and cdf is known, otherwise contains -1;
    end
    properties (SetAccess=private) % properties only assessible by routines of the object
        EffectiveSumObs        % effective sum of observations, mu0 * t0. used for improving numerical stability
    end
    
    methods
        function obj = DistNormalMuSig(mu0, t0, xi0, chi0, defaultpar) 
           %Hyperparameter for Normal(W,Sig^2) with unknown mean and unknown
           %variance is taken to be the conjugate normal-invgamma prior
           % W | sigma ~ Normal(mu0,sigma^2/t0) 
           % Sigma^2 ~ InvGamma(xi0, chi0), with mean chi0/(xi0-1)
           % Default mu0 is 0. Default t0 is 4 (must be greater than 0).
           % Default chi0 is 16, default xi is 4.
           % Thus, (mu0, t0, xi0, chi0) is the hyperparameter.
           % defaultpar should be a column vector with a 'default' realized 
           % value of the unknown (mu, sigma^2), 
           % and can be used to generate samples. Defaults to empty value.
           % (can be modified with SampleThetaBayes or SampleThetaFreq).
            if nargin<1
                obj.mu0 = 0;        % default hypermean is 0
            else
                obj.mu0 = mu0;
            end
            if nargin<2
                obj.t0 = 4;        % default effective number of samples is 4
            else
                if t0 <= 0
                    warning('t0 must be greater than 0, reset to 4');
                    t0 = 4;
                end
                obj.t0 = t0;
            end
            if nargin<3
                obj.xi0 = 4;      % default xi0 is 4
            else
                if xi0 <= 1
                    warning('xi0 must be greater than 1, reset to 4');
                    xi0 = 4;
                end
                obj.xi0 = xi0;
            end
            if nargin<4
                obj.chi0 = 16;      % default xi0 is 4
            else
                if chi0 <= 0
                    warning('xi0 must be greater than 0, reset to 16');
                    chi0 = 16;
                end
                obj.chi0 = chi0;
            end
            if nargin<5
                obj.Thetavec = [];    % initialize the vector of sampled thetas to be the empty vector
            else
                obj.Thetavec = defaultpar(:); % unless it is specified (force it to be column vector)
            end
            obj.Samplecdfmat = [];
            % *TO ADAPT*: To set up a different distribution, make sure that
            % there are properties for all of the hyperparameters.
            % For sampling distributions with multiple parameters, each dimension of
            % the parameter for the sampling distribution should be in a
            % separate row of thetavec.
        end

         function obj = SetMean(obj, newmean) 
           %Reset hyperparameter so that the mean output is newmean, 
           if nargin == 2
               obj.mu0 = newmean;
           end
            % *TO ADAPT*: To set up a different distribution, it might be necessary
            % to reset several parameters in order to get the output to have
            % the specified mean: for example, for lognormal, the mean
            % output is a function of both the mean and variance of the
            % underlying normal distribution that generates the lognormal
            % random variable.
        end

        function ESample = get.ESample(obj)
            % If Thetavec is non-empty, then return a vector of means for
            % those parameters. if it is empty, then return the marginal
            % mean of samples, conditional on hyperparameters.
            tmp = obj.Thetavec;
            if length(tmp) == 0
                ESample = obj.mu0;
            else
                ESample = tmp(1,:);
            end
            % *TO ADAPT*: Implement a function which returns the mean value
            % of the sampling distribution.
        end

        function obj = set.ESample(obj,~)
            fprintf('%s%f\n','Mean is: ',obj.ESample)
            error('You cannot set ESample explicitly. Instead use SetMean method.');
        end

        function EVariance = get.EVariance(obj)
            EVariance = (obj.chi0 / obj.xi0) * ((obj.t0+1) / obj.t0);
            % *TO ADAPT*: Implement a function which returns the mean value
            % of the sampling distribution.
        end

        function obj = set.EVariance(obj,~)
            fprintf('%s%f\n','EVariance is: ',obj.EVariance)
            error('You cannot set EVariance explicitly. Instead use SetMean method or set hyperparameter properties');
        end

        function obj = set.mu0(obj,mu0)
            if nargin >= 2
                obj.mu0 = mu0;
                obj.EffectiveSumObs = obj.t0 * obj.mu0; % hidden, used to improve numerical stability of updating
            end
        end
        
        function obj = set.t0(obj,t0)
            if nargin >= 2
                obj.t0 = t0;
                obj.EffectiveSumObs = obj.t0 * obj.mu0; % hidden, used to improve numerical stability of updating
            end
        end

       %%%%%%%%%%%%%% ROUTINES FOR SAMPLING FROM PRIOR AND CONDITIONAL ON THETA %%%%%%%%%%%%%%

       function samps = Sample(obj,n,u)
           % Given the distribution object 'obj', generate samples from
           % that distribution accordig to its parameter values. That is,
           % sample X given obj.<parameters>.
           % by default, return a single sample for each theta in thetavec
           % If optional n are included, return an nxm matrix of samples,
           % where m is the length of thetavec.
           % If optional u is included, try to generate samples with invert
           % the cumulative method, using uniforms in u. assumes that u is
           % an nxm matrix of uniform (0,1) rv. Useful for CRN.
           % if u is not present, then generate samples independently. 
           m = length(obj.Thetavec);
           if nargin < 2   
               n=1;    
           end
           if nargin < 3
               u = -1;
           else
               [um, un] = size(u);
               if m > un   % if we have more uniforms than needed to cover realized parameters, it is ok, otherwise warn
                   warning('need length of u and n to be the same to get CRN to work, sampling independently')
                   u=-1;
               end
               if um > n   % if we are asking for more samples per parameter than we have uniforms to handle, then warn
                   warning('need length of u and n to be the same to get CRN to work, sampling independently')
                   u=-1;
               end
           end
           if m > 0
               samps = zeros(n,m);  % preallocate memory for the samples
               if u == -1
                    u = rand(n,m);
                    obj.Samplecdfmat = u;   % keep CDF of samples in case one wants to use CRN after
               end
               for i=1:m
                   samps(:,i) = SampleWorker(obj,i,n,u(:,i));
               end
           else
                warning('no realizations of parameters found: call SampleThetaBayes() or SampleThetaFreq() before Sample()')
               samps = [];      % if there are no thetas, then can't sample from them
               obj.Samplecdfmat = []; % reset so CRN for samples is goned
           end
            % *TO ADAPT*: This method does not need much adaptation. Instead, see SampleWorker method below. 
       end

       function obj = SampleThetaBayes(obj,m,u)
           % Given the distribution object 'obj', generate samples from the
           % prior distribution (given the hyperparamters). That is,
           % sample theta given obj.<hyperparameters>
           % by default, return a single theta
           % If optional m is included, return an dxm matrix of theta,
           % where d is the dimension of the unknown parameters to be
           % sampled. (e.g., if mean is unknown and variance is known then
           % d=1; if both are unknown then d=2).
           % If optional u is included, try to generate thetas with invert
           % the cumulative method, using uniforms in u. assumes that u is
           % an mxn matrix of uniform (0,1) rv. Useful for CRN. Is only
           % used when it is possible to use it - if CRN is not easy to
           % implement then return independent realizations.
           % if u is not present, then generate thetas independently. 
           if nargin < 2   
               m=1;
           end
           if m == 0  % allow ability to reset the set of sampled theta to the empty vector
               obj.Thetavec = [];
               obj.Thetacdfvec = -1;
           else
               if nargin < 3
                   u = -1;
               else
                   [~, un] = size(u);
                   if (un ~= m) && (un ~= 1)
                       u = -1;
                       warning('SampleThetaBayes: m, u not of compatible sizes');
                   end
               end
               obj.Thetavec = SampleThetaWorker(obj,m,u);
           end
            % *TO ADAPT*: This method does not need much adaptation. Instead, see SampleThetaWorker method below. 
       end

       function obj = SampleThetaFreq(obj,m,mu)
           % generate a dxm vector of values of theta, assuming the
           % default frequentist value for theta, given the hyperprior,
           % where d is the dimensionality of the unknown parameter(s) of 
           % the sampling distribution.
           % Here, theta should be a column vector with the MEAN and VARIANCE.
           if nargin < 2   
               m=1;
           end
           if nargin < 3   
               mu = obj.mu0;
           end
           if m == 0
               obj.Thetavec = [];
           else
               varest = obj.chi0 / obj.xi0;
               obj.Thetavec = [mu*ones(1,m); varest*ones(1,m)];
           end
           obj.Thetacdfvec = -1;
       end

       %%%%%%%%%%%%%% ROUTINES FOR INFERENCE %%%%%%%%%%%%%%

       function obj = BayesUpdate(obj,x)
           % Implement Bayes' rule, when observing 'x', the input (prior
           % distribution) is updated and the output returned is the
           % posterior distriution.
           %
           % for now, assumes obj is univariate, but x might be a vector of
           % samples
           %
           % at present, implementation uses hidden property
           % EffectiveSumObs in order to keep track of the sum of
           % observations, with goal of improving numerical stability for
           % update process
           FIX
           if length(x) > 0
               tmpsum = obj.EffectiveSumObs + sum(x); % do the updates on the local version
               obj.t0 = obj.t0 + length(x); % reassign to hyperparameter
               obj.mu0 = tmpsum / obj.t0;            % that reassignment should reset the key statistics
           end
       end
       
       %%%%%%%%%%%%%% ROUTINES FOR OUTPUT %%%%%%%%%%%%%%
        function displayObjectName(self)
            disp(inputname(1))
        end       
       
        function disp(obj)
            [m, n] = size(obj);
            if m*n == 1
                fprintf(1,'Distribution: %s\n',mfilename('class'));
                fprintf(1,'Hyperparameters:\n  mu0=%f, t0=%f, xi0=%f, chi0=%f\n',...
                    obj.mu0,obj.t0,obj.xi0,obj.chi0);
                fprintf(1,'Parameters:\n');
                obj.Thetavec
            else
                for i=1:m
                    for j=1:n
                        fprintf(1,'Distribution: %s (size %d x %d)\n',mfilename('class'),m,n);
                        fprintf(1,'Hyperparameters:\n  mu0=%f, t0=%f, xi0=%f, chi0=%f\n',...
                            obj(i,j).mu0,obj(i,j).t0,obj(i,j).xi0,obj(i,j).chi0);
                        obj(i,j).Thetavec
                    end
                end
            end
        end % disp
    
       %%%%%%%%%%%%%% ROUTINES FOR TESTING %%%%%%%%%%%%%%
        function obj = TestCase(obj)
            SampleThetaBayes(obj,5)        % Generate 5 independent realizations of parameters from the prior distribution defined by the hyperparaameters 
            meanval = obj.ESample
            varvar = obj.EVariance
            SampleThetaBayes(obj,5)        % Generate 5 more independent realizations of parameters from the prior distribution defined by the hyperparaameters 
            obj.Thetavec                   % display those realizations
            obj.Thetacdfvec                % display the CDF of the sampled thetas with respect to the prior distribution for those values.

            samples = Sample(obj,10)         % Generate 10 independent observations conditional on reaziations of parameters, for each of the realized parameter values

            samples = Sample(obj,10)         % Generate 10 more independent observations conditional on reaziations of parameters, for each of the realized parameter values
            umat = obj.Samplecdfmat;       % get the CDF of the most recently sampled values

            meanval = obj.ESample
            SampleThetaBayes(obj,5)        % Generate 5 more independent realizations of parameters from the prior distribution defined by the hyperparaameters 
            SetMean(obj, 5)               % reset the mean to 5
            SampleThetaBayes(obj,5,obj.Thetacdfvec) % Generate 5 realizations of parameters with the new mean, which are correlated to the theta that were previously generated
            meanval = obj.ESample
            Sample(obj,10,umat)           % generate samples with those new parameters, which are correlated with the samples previously generated for the other theta
        end       
 
    end

    methods (Access = 'private') % Access by class members only
       function sampledthetas = SampleThetaWorker(obj,m,u)
           % Worker routine for Sample method above. This is where the
           % 'guts' of generating samples from a given distribution
           % happens. If u is -1, then generate samples independently. If u
           % is a dxm matrix of uniforms, try to generate using invert the
           % cumulative if possible (if not possible, then try some other
           % method to take advantage of CRN, or, if all else fails,
           % generate independent randoms)
           if (nargin < 3)
               u = -1;
           end
           if u == -1   % generate independent random numbers
               u = rand(2,m);
           end % thetavec must have parameters in each column, number of rows is number of dimensions of parameter
           obj.Thetacdfvec = u;
           
%           sampledvars = 1./gamrnd(obj.xi0,1/obj.chi0,1,m);
%           sampledmeans = norminv(u,obj.mu0,sqrt(sampledvars/obj.t0));
           sampledvars = 1./gaminv(u(2,:),obj.xi0,1/obj.chi0);
           sampledmeans = norminv(u(1,:),obj.mu0,sqrt(sampledvars/obj.t0));
           sampledthetas = [ sampledmeans; sampledvars];
           % *ADAPT*: This is for generating one col per parameter from
           % prior distribution, one row for each dimension of unknown parameter
       end
       
        function samps = SampleWorker(obj,i,n,u)
           % Worker routine for Sample method above. This is where the
           % 'guts' of generating samples from a given distribution
           % happens. If u is -1, then generate samples independently. If u
           % is a matrix of uniforms, try to generate using invert the
           % cumulative if possible (if not possible, then try some other
           % method to take advantage of CRN, or, if all else fails,
           % generate independent randoms)
           % otherwise, generate n samples (in column vector form) conditional
           % on the i-th parameter in Thetavec.
           % assumes error checking has been done so that u is a nx1 matrix
           if u == -1   % generate independent random numbers
               u = rand(n,1); 
           end % samps must be a row vector, hence u(:)'
           samps = norminv(u(:),obj.Thetavec(1,i),sqrt(obj.Thetavec(2,i))); % *ADAPT*
           obj.Samplecdfmat(:,i) = u;
           % *ADAPT*: This is for generating a column vector of samples, given
           % the i-th parameter value
        end
       
    end
end
