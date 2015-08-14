function figout = DoReversalPlot(fignum,basic,advanced,mat,basiccomp,advancedcomp,matcomp,graphicsuffix)

dirname = 'Figure'; % directory for putting figures

plot_upper = 6000 ; % set upper limit of plot region ( + / 
plot_lower = -6000  ; % set lower limit of plot region

maxval = max(max(mat.outBayes.PrReversalAve),max(mat.outBayes.PrReversalAve));
maxval = ceil(10*maxval)/10;

fignum = fignum + 1; figure(fignum);
hold off ; 
ishoriz = true; % set to true if muvec is on horizontal axis, false if on vertical axis
ist0 = false ;
%fignum = fignum + 1; figure(fignum);
plot(advanced.simFreqDeltaVec, mat.outBayes.PrReversalAve, '-k');
hold on ; 
plot(advancedcomp.simFreqDeltaVec, matcomp.outBayes.PrReversalAve, '--k');
xlim( [ plot_lower, plot_upper] ) ;
hold on ; 
xlabel( 'Prior mean', 'Fontsize', advanced.bigfontsize ) ; 
ylabel( 'Prob(reversing decision with pipeline data)', 'Fontsize', advanced.bigfontsize ) ; 
legend( 'Baseline',  'Comparator', 'Location', 'northeast' ) ;
hold on ; 
lohivec = [ 0.045, maxval ] ;
UtilPlotABCD( basic, advanced, mat, ishoriz, ist0, lohivec ); 
lohivec = [ 0.004, maxval ] ;
UtilPlotABCD( basiccomp, advancedcomp, matcomp, ishoriz, ist0, lohivec ); 
UtilStdizeFigure(fignum,advanced);
UtilSaveFigEpsPdf(fignum,dirname, strcat('pr_reversal_same',graphicsuffix),'-r600');

figout=fignum;

end
