%This script will walk you though how to use the IRbcd function to set up
%and run a test problem. We begin by intializing the problem size.
n=64;
%PRset_var is a variation on PRset. It works the same way. Here we are
% setting a new variable Rvar, which contains the Rvalues. See the file
% Documentation.m for the other new parameters that can be passed in.

rng(5) %Set the random number generator for reproducibility. 

% Here in declare an image string avaliable in PRset, or one of your own.
% Must be square.
image = imresize(double(imread('cameraman.tiff')),[64,64]);

%PRset sets the options for the tomography problem. 
ProbOptions = PRset('Rvar',2 * ones(1,4),'Rpert',0.25,'anglespert',0.25,'phantomImage',image);

%The new PRtomo_var only returns the matrix b and probInfo. See its
%documentation for everything contained in probInfo
[b,probInfo] = PRtomo_var(n,ProbOptions);

% Here we initialize the IR options. See the file Documentation.m or IRbcd
% for a list of the options avaliable. 
iterOptions = IRset('nonlinSolver','imfil','BCDmaxIter',20,...
    'Rbounds',0.1250,'angleBounds',0.1250,'BCDlsSolver','fista');

%This function runs the BCD loop. See the IRbcd documentation for what
%iterInfo contains.
[x,iterInfo] = IRbcd(b,iterOptions,probInfo);

%Lets see how it did!
PRshowbcd(iterInfo,probInfo)