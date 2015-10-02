function [ rval ] =  UtilSaveVersionFile( dirname )
%UtilSaveVersionFile saves code version information, assuming that SVN was used to localize the code.
%
%   To locally adapt this, the path to where the svn command is found should be edited.
%
%   
% For use in Delay sequential sampling project.
if ~isempty(dirname)      
    tmp = dir(dirname);
    if ~length(tmp) %~isdir(dirname)
        mkdir(dirname);
    end
    fullpath = sprintf('%s/',dirname);
else
    fullpath = '';
end

versionfile = [fullpath 'version.txt'];
% This looks for svn first assuming svn is a known command, then if the.
% seems to be ok if svn is installed with cygwin. not sure yet how to get
% it to work with tortoisesvn
[errsvn, svncommand] = system('where svn');
if errsvn % if could not find svn, it might be entered manually below. E.g. svn might be installed with cygwin instead of in cmd environment
   svncommand = 'c:\cygwin64\bin\svn '; % for cygwin64's svn default installation
   [errcode, codeversion] = system([svncommand ' info']); % get information about version of code in repo
   if errcode
       svncommand = '/usr/bin/rapidsvn.sh ';    % for rapidsvn default installation
       [errcode, codeversion] = system([svncommand ' info']); % get information about version of code in repo
   end
else
    [errcode, codeversion] = system([svncommand ' info']); % get information about version of code in repo
end

[errstatus, codestatus] = system([svncommand ' status']); % check if there are any local uncommited changes to the repo
if ~(errcode + errstatus)
    fid = fopen(versionfile,'w');
    codeversioninfo = sprintf('%s\n%s\n',codeversion, codestatus);
    rval = fprintf(fid,'%s\n',codeversioninfo);
    fclose(fid);
else
    rval = false;
end

end
