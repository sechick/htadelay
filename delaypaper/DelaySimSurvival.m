function[ figout ] = DelaySimSurvival( fignum, basic, advanced, mat ) 
% Takes a vector of failure times from the simulation and performs nonparametric survival
% analysis, producing estimates of the cdf and the hazard function. Suppplements the 
% histograms for outtime that SC has already produced in DelaySimComputer. 
% Project with Chick, Forster, Pertile (alpha order)
% Written 31/03/2015 MF

%Some test commands which can be used to work directly with .mat output
%clear ; 
%load('c:\martin\chick\hta\trunk\delaypaper\Sinus.mat') ; 

TMax = basic.TMax ; 

xx = linspace(1,(TMax-1),TMax) ;

% preallocate structure
so = struct ;

for i=1:length( advanced.simFreqDeltaVec )
    simnumber = i ;
    so( i ).mu0 = advanced.simFreqDeltaVec( i ) ;
    % pick the vector of outtimes 
    obstime = mat.simBayesOut(i).ENumSeen ;
    obs_max= max( obstime ) ; 
    obs_min = min(obstime ) ; 
    %obstime_max = max( obstime ) ; 
    % calculate the number of `failures' (arrive at a decision) and number
    % of `censoreds' reach TMax witout a decision
    failed = obstime( obstime < TMax - 1 ); 
    nfailed = sum( failed ) ;
    censored = ( obstime >= TMax - 1 ); 
    ncensored = sum( censored );
    % empirical cdf (Kaplan-Meier estimate of CDF)
    [empF,x,empFlo,empFup] = ecdf( obstime, 'censoring' , censored );
    [empS,x,empSlo,empSup] = ecdf( obstime, 'censoring' , censored, 'Function', 'survivor' ) ;
    % Kernel density estimate of the hazard function 
    [npF,ignore,u] = ksdensity(obstime,xx,'cens',censored,'function','cdf');
    npF3 = ksdensity(obstime,xx,'cens',censored,'function','cdf','width',u/3);
    hazrate = ksdensity(obstime,xx,'cens',censored,'width',u/3) ./ (1-npF3);
    so(i).simnumber = simnumber ;
    so(i).obstime = obstime ;
    so(i).obsmax = obs_max ;
    so(i).obsmin = obs_min ; 
    so(i).nfailed = nfailed ;
    so(i).ncensored = ncensored ;
    so(i).empF = empF ;
    so(i).empS = empS ;
    so(i).x = x ; 
    so(i).xx = xx ;
    so(i).hazrate = hazrate ; 
end

fignum = fignum + 1; figure(fignum); % always increment fignum and create a new figure

% Compare two estimates of hazard and survivor functions, one on either
% side of mu_0 = 0. Choose which positions you would like to compare. See
% Matlab's HELP files for the code. 
upper_cell = 156; 
lower_cell = 126 ; 
hold off ; 
fignum( ) ; 
leg1 = num2str( so(upper_cell).mu0 ) ;
leg2 = num2str( so(lower_cell).mu0 ) ; 
plot( so(upper_cell).xx, so(upper_cell).hazrate, 'k', so(lower_cell).xx, so(lower_cell).hazrate, '--k' ) ; 
hold on; 
legend( ['mu0=' leg1], ['mu0=' leg2] );
xlim([0,TMax]) ; 
ylim([0,1/10]) ;
title( 'Kernel density estimate of hazard functions' ) ; 
%pause ; 
fignum = fignum + 1; figure(fignum);
hold off ;
fignum( ) ; 
plot( so(upper_cell).x, so(upper_cell).empS, 'k', so(lower_cell).x, so(lower_cell).empS, '--k' ) ; 
legend( ['mu0=' leg1], ['mu0=' leg2] );
xlim([0,TMax]) ;
title( 'Empirical Survival Function ' ) ; 
xlabel( 'Sample size' ) ; 
ylabel( 'Proportion not failed' ) ; 

figout = fignum;

end





