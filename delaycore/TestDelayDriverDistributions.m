% Test drives the random distribution object classes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

myrvtype = @DistNormalMu        % define a variable whose value is a pointer to the random object creator's class
tst = myrvtype()                % generate a default instance: this will also allow one to see what parameters are defined for the class 
myrv = myrvtype(0, 10, 20)      % create an instance of the random variable. with different values for its hyperparameters
TestCase(myrv)              % run some diagnostic tests for this class of random variables

myrvtype = @DistNormalMuSig        % define a variable whose value is a pointer to the random object creator's class
tst = myrvtype()                % generate a default instance: this will also allow one to see what parameters are defined for the class 
myrv = myrvtype(0, 5, 10, 100)      % create an instance of the random variable. with different values for its hyperparameters
SampleThetaBayes(myrv,10000);
mean(myrv.Thetavec')
var(myrv.Thetavec')

TestCase(myrv)              % run some diagnostic tests for this class of random variables

dirData = dir();      %# Get the data for the current directory
dirIndex = [dirData.isdir];  %# Find the index for directories
fileList = {dirData(dirIndex).name};  % Get list of directories
definedRVs={fileList{strncmp(fileList, '@Dist', 4)}}  % set of RVs available should include those in directories starting with '@Dist'

for i=1:length(definedRVs)
    myrvtype = str2func(definedRVs{i})        % loop through each type of RV which is available
    myrv = myrvtype()               % create an instance of that RV
    TestCase(myrv)                  % run some diagnostics for that RV
    if i<length(definedRVs)
        pause
    end
end


