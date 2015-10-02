function [rval, basic,advanced] = DelayInputModifier(basic, advanced, basicarray, advancedarray)
%DELAYINPUTMODIFIER 

% Designed for Delay Sequential Trials paper of Chick, Forster, Pertile
% (alpha order).
%
% (c) 2014, S Chick
% Created: 14 April 2014
% Last touched: 14 April 2014
% 
% The first two batches of parameters are most critical in terms of defining the
% problem structure of the clinical trial.
%
if nargin < 4
    advancedarray = {};
end
if nargin < 3
    basicarray = {};
end

basicarrlen = length(basicarray);
advarrlen = length(advancedarray);
rval = 1;
for i=1:(basicarrlen/2)
    if isfield(basic,basicarray{2*i-1}) || strcmp(basicarray{2*i-1},'mumax') ||  strcmp(basicarray{2*i-1},'mumin')
        basic.(basicarray{2*i-1}) = basicarray{2*i};
    else
        warning(sprintf('invalid basic field: %s',char(basicarray{2*i-1}))); 
%        basicarray{2*i-1}
        rval = 0;
    end
end
for i=1:(advarrlen/2)
    if isfield(advanced,advancedarray{2*i-1}) || ...
            strcmp(advancedarray{2*i-1},'dmu') || ...
            strcmp(advancedarray{2*i-1},'simFreqDeltaVec') 
        advanced.(advancedarray{2*i-1}) = advancedarray{2*i};
    else
        warning(sprintf('invalid advanced field: %s',char(basicarray{2*i-1}))); 
%        'invalid advanced field: ' 
%        advancedarray{2*i-1}
        rval = 0;
    end
end

end

