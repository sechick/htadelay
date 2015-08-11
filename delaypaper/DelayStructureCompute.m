function[basic, advanced]=DelayStructureCompute(basic,advanced)


basic.tau = basic.delayinyears * basic.patientsperyear ;
basic.theta = exp( - log(1+basic.annualdiscountrate) / basic.patientsperyear ); %  per patient discount rate

% basic.sigma0 = basic.sigma/sqrt(basic.t0) ; no need to assign, will be
% computed
