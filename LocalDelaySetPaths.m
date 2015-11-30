function [ rval ] = LocalDelaySetPaths( basepath )
% LocalDelaySetPaths: 
%   PURPOSE: updates the matlab path variable so that the code for the
%   paper on Bayesian Sequential Optimization with Delayed Samples' (Chick,
%   Forster, Pertile). To be called before running code associated with
%   that paper.
%   INPUT: basepath is either empty (to default to the current directory in
%   the matlab environment) or a text string with the path name associated
%   with the routines. Do not terminate the directory name with a \ on
%   windows systems.
%   OUTPUT: always returns true
%   EXAMPLE USAGE:
%       LocalDelaySetPaths();
%       LocalDelaySetPaths('d:\users\papers\Forster\SeqDelay.git');
%   INTENDED WORK FLOW:
%       When downloading matlab code for the paper, this file should be
%       updated once to adapt to the local machine's requirements. This
%       update is not to be committed back to the code repo.
%
% Source provided 'as is' with no warrantees or claims provided or implied.
% 2015 S Chick

    if nargin < 1
        BASEDIR = pwd;  
    else
        BASEDIR = basepath; % 
    end

    %BASEDIR = 'd:\users\papers\Forster\hta\trunk\'

    delaycoredir = [ BASEDIR '\delaycore\'];
    delaypaperdir = [ BASEDIR '\delaypaper\'];
    delayunkdir = [ BASEDIR '\delayunk\'];
    basedelaydir = [BASEDIR];    

    addpath(genpath(delaycoredir));
    addpath(genpath(delaypaperdir));
    addpath(genpath(delayunkdir));
    addpath(basedelaydir);              % set up this way to preclude .git directory structure from being added to path
    
    rval = 1;
end