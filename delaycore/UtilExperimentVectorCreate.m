function [basicvec, advancedvec, legendvec, veclen] =  UtilExperimentVectorCreate( basic, advanced, basicflag, fieldname, fieldvec)
%UtilCreateExperimentVector creates a vector of basic and advanced structures, where a given field is varied
%in the elements of the vector by setting a given field to have the values in fieldvec.

% do a sanity check for the field name
exitflag = 0;
if basicflag
    if ~isfield(basic,fieldname)
        'Error: field not found in basic structure'
        exitflag=1;
    end
else
    if ~isfield(advanced,fieldname)
        'Error: field not found in advanced structure'
        exitflag=1;
    end
end

advanced.iterateextratext = '';   % for file name differentiation
if ~exitflag
    veclen = length(fieldvec);
    for i=1:veclen
       % make the return vector a bit larger
       legendval = sprintf('%s = %s',fieldname, num2str(fieldvec(i)));
       if i==1 
           basicvec = basic;
           advancedvec = advanced;
           legendvec = {legendval};
       else
           basicvec = [basicvec basic];
           advancedvec = [advancedvec advanced];
           legendvec = [legendvec {legendval}];
       end
       advancedvec(i).iterateextratext = sprintf('%s%d',fieldname, num2str(fieldvec(i)));   % for file name differentiation
       % set the values appropriately
       if basicflag
           basicvec(i).(fieldname) = fieldvec(i);
       else
           advancedvec(i).(fieldname) = fieldvec(i);
       end
       [basicvec(i), advancedvec(i), rval, msgs] = DelayInputValidator( basicvec(i), advancedvec(i) );
    end
end

end
