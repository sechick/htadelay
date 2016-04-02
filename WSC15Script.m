LocalDelaySetPaths % Edit a local copy of this MACRO, don't check in on top of installed file.

NUMSIMREPS = 10;
NUMGRIDCONTOUR = 50;

NUMSIMREPS = 10000;
NUMGRIDCONTOUR = 120;


%%%%%%%% TEST KNOWN VERSUS UNKNOWN VARIANCE %%%%%%%%
if ~exist('fignum','var'), fignum = 20; end;
[basic, advanced] = SetStents_WSC();
advanced.filestring='WSCUKvar'
advanced.DOPDE = true;
advanced.NumGridsForContours = NUMGRIDCONTOUR;
advanced.simNumReps = NUMSIMREPS;
[basic, advanced, rval, msgs] = DelayInputValidator( basic, advanced );
if ~rval, msgs, end;
fieldname = 'UnkVariance';       % pick the name of the field whose values are to be changed
basicflag = false;           % set to true if the field is in the 'basic' structure, false if it is in the 'advanced' structure
fieldvec = [true false];    % create an array of values to rotate through for that field
subtitle = 'Known and unknown variance';    % give a short subtitle name for the figures
advanced.verbose = true;
fmodifier = 'UnkVec';      % used to diferentiate file name if it is used for saving plots
if advanced.verbose % the 'true' part of this keeps the various variables around. the 'else' part has it implemented in a subroutine, and thus does not keep the output for futher use
    % create vectors of basic and advanced structures
    [basicvec, advancedvec, legendvec, veclen] = UtilExperimentVectorCreate( basic, advanced, basicflag, fieldname, fieldvec);
    % run the analysis for each of the structures
    for i=1:veclen
%        advancedvec(i).DOPDE = ~advancedvec(i).UnkVariance; % only do pde of variance is known
        [~, mat] = DelayCurvesRecur(basicvec(i), advancedvec(i));        % Do stage II of DP
        [mat] = DelayStageOne(basicvec(i), advancedvec(i), mat );           % Do stage I of DP
        [ fignum, mat ] = DelaySimOverview( fignum, basicvec(i), advancedvec(i), mat );
%        fignum = DelaySimOutput( fignum, basicvec(i), advancedvec(i), mat );     % requires DelaySimOverview to have been called
        if i==1 matvec = mat; else matvec = [matvec mat]; end
    end
    % do the plots based on the structures
    [fignum] = UtilExperimentVectorPlot( fignum, basicvec, advancedvec, legendvec, matvec, subtitle, fmodifier);
    for i=1:veclen
        advancedvec(i).subtitle = legendvec{i};
        advancedvec(i).filestring = [advancedvec(i).filestring num2str(i)];
%        fignum = DelayDoContours(fignum, basicvec(i), advancedvec(i), matvec(i));        % Generate a bunch of contour plots
    end
else
    [fignum, basicvec, advancedvec, legendvec, matvec] = TestDelayIterate(fignum,basic,advanced,basicflag,fieldname,fieldvec,subtitle,fmodifier);
    [fignum] = UtilExperimentVectorPlot( fignum, basicvec, advancedvec, legendvec, matvec, subtitle,fmodifier);
end

%%%%%%%% TEST PDE VERSUS KG*, KNOWN VARIANCE %%%%%%%%
if ~exist('fignum','var'), fignum = 20; end;
[basic, advanced] = SetStents_WSC();
advanced.filestring='WSCKvarKGs'
advanced.DOPDE = false;
advanced.UnkVariance = false;
advanced.NumGridsForContours = NUMGRIDCONTOUR;
advanced.simNumReps = NUMSIMREPS;
%basic.theta = 1;
[basic, advanced, rval, msgs] = DelayInputValidator( basic, advanced );
if ~rval, msgs, end;
fieldname = 'DOPDE';       % pick the name of the field whose values are to be changed
basicflag = false;           % set to true if the field is in the 'basic' structure, false if it is in the 'advanced' structure
fieldvec = [false true];    % create an array of values to rotate through for that field
subtitle = 'KG* and PDE (known var)';    % give a short subtitle name for the figures
advanced.verbose = true;
fmodifier = 'KGs';      % used to diferentiate file name if it is used for saving plots
if advanced.verbose % the 'true' part of this keeps the various variables around. the 'else' part has it implemented in a subroutine, and thus does not keep the output for futher use
    % create vectors of basic and advanced structures
    [basicvec, advancedvec, legendvec, veclen] = UtilExperimentVectorCreate( basic, advanced, basicflag, fieldname, fieldvec);
    % run the analysis for each of the structures
    for i=1:veclen
        [~, mat] = DelayCurvesRecur(basicvec(i), advancedvec(i));        % Do stage II of DP
        [mat] = DelayStageOne(basicvec(i), advancedvec(i), mat );           % Do stage I of DP
        [ fignum, mat ] = DelaySimOverview( fignum, basicvec(i), advancedvec(i), mat );
%        fignum = DelaySimOutput( fignum, basicvec(i), advancedvec(i), mat );     % requires DelaySimOverview to have been called
        if i==1 matvec = mat; else matvec = [matvec mat]; end
    end
    % do the plots based on the structures
    [fignum] = UtilExperimentVectorPlot( fignum, basicvec, advancedvec, legendvec, matvec, subtitle, fmodifier);
    for i=1:veclen
        advancedvec(i).subtitle = legendvec{i};
        advancedvec(i).filestring = [advancedvec(i).filestring num2str(i)];
%        fignum = DelayDoContours(fignum, basicvec(i), advancedvec(i), matvec(i));        % Generate a bunch of contour plots
    end
else
    [fignum, basicvec, advancedvec, legendvec, matvec] = TestDelayIterate(fignum,basic,advanced,basicflag,fieldname,fieldvec,subtitle,fmodifier);
    [fignum] = UtilExperimentVectorPlot( fignum, basicvec, advancedvec, legendvec, matvec, subtitle,fmodifier);
end

%%%%%%%% TEST PDE VERSUS KG*, KNOWN VARIANCE %%%%%%%%
if ~exist('fignum','var'), fignum = 20; end;
[basic, advanced] = SetStents_WSC();
advanced.filestring='WSCUvarKGs'
advanced.DOPDE = false;
advanced.UnkVariance = true;
advanced.NumGridsForContours = NUMGRIDCONTOUR;
advanced.simNumReps = NUMSIMREPS;
%basic.theta = 1;
[basic, advanced, rval, msgs] = DelayInputValidator( basic, advanced );
if ~rval, msgs, end;
fieldname = 'DOPDE';       % pick the name of the field whose values are to be changed
basicflag = false;           % set to true if the field is in the 'basic' structure, false if it is in the 'advanced' structure
fieldvec = [false true];    % create an array of values to rotate through for that field
subtitle = 'KG* and PDE (unknown var)';    % give a short subtitle name for the figures
advanced.verbose = true;
fmodifier = 'KGsUnk';      % used to diferentiate file name if it is used for saving plots
if advanced.verbose % the 'true' part of this keeps the various variables around. the 'else' part has it implemented in a subroutine, and thus does not keep the output for futher use
    % create vectors of basic and advanced structures
    [basicvec, advancedvec, legendvec, veclen] = UtilExperimentVectorCreate( basic, advanced, basicflag, fieldname, fieldvec);
    % run the analysis for each of the structures
    for i=1:veclen
        [~, mat] = DelayCurvesRecur(basicvec(i), advancedvec(i));        % Do stage II of DP
        [mat] = DelayStageOne(basicvec(i), advancedvec(i), mat );           % Do stage I of DP
        [ fignum, mat ] = DelaySimOverview( fignum, basicvec(i), advancedvec(i), mat );
%        fignum = DelaySimOutput( fignum, basicvec(i), advancedvec(i), mat );     % requires DelaySimOverview to have been called
        if i==1 matvec = mat; else matvec = [matvec mat]; end
    end
    % do the plots based on the structures
    [fignum] = UtilExperimentVectorPlot( fignum, basicvec, advancedvec, legendvec, matvec, subtitle, fmodifier);
    for i=1:veclen
        advancedvec(i).subtitle = legendvec{i};
        advancedvec(i).filestring = [advancedvec(i).filestring num2str(i)];
%        fignum = DelayDoContours(fignum, basicvec(i), advancedvec(i), matvec(i));        % Generate a bunch of contour plots
    end
else
    [fignum, basicvec, advancedvec, legendvec, matvec] = TestDelayIterate(fignum,basic,advanced,basicflag,fieldname,fieldvec,subtitle,fmodifier);
    [fignum] = UtilExperimentVectorPlot( fignum, basicvec, advancedvec, legendvec, matvec, subtitle,fmodifier);
end

save WSC15script.mat
