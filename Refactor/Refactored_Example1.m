%This script will walk you though how to use the IRbcd function to set up
%and run a test problem. We begin by intializing the problem size.
n=64;
%PRset_var is a variation on PRset. It works the same way. Here we are
% setting a new variable Rvar, which contains the Rvalues. See the file
% Documentation.m for the other new parameters that can be passed in.
ProbOptions = PRset_var('Rvar',2 * ones(1,4));

%The new PRtomo_var only returns the matrix b and probInfo. See its
%documentation for everything contained in probInfo
[b,probInfo] = PRtomo_var(n,ProbOptions);

% Here we initialize the IR options. See the file Documentation.m or IRbcd
% for a list of the options avaliable. 
iterOptions = IRset('nonlinSolver','imfil','BCDlsSolver','cgls','BCDmaxIter',10);

%This function runs the BCD loop. See the IRbcd documentation for what
%iterInfo contains.
[x,iterInfo] = IRbcd(b,iterOptions,probInfo);

%Lets see how it did!
PRshowx(x,probInfo.ProbParams)