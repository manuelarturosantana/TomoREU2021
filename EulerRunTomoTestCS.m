tic
%
% This script sets up test problem with changes in R and theta. 
% To run different versions of the problem you just need to change the values
% in the first few lines. Here are some comments about the variables:
%  m               = number of times R changes 
%                   (in PRtomo_var this is also called m, and is the length 
%                    of vector R, or number of columns in the array angles)
%  Rguess          = guess that R is always equal to this value for all angles
%  Rtrue           = true values of R, which may be slight perturbations of Rguess
%  Rnoise          = scalar constant on the amount of perturbation added to Rguess
%  Rnoise_guess    = The intial guess for Rnoise.
%  Rpert           = A vector containing the actual perterbations on R.
%  ang_noise_guess = The initial guess for the angles parameter used in the
%                    BCD.
%  ang_noise       = scalar constant on the amount of perturbation added to
%                    ang_true
%  ang_guess       = The guess of what the projection angles are.
%  angles_true     = true projection angles, a slight perturpation of
%                    ang_guess.
%  p               = [total number of projection angles]/m, should be an 
%                    integer (see also PRtomo_var, where p is the number of
%                    rows in the array angles)
%  span            = scalar that determines the angular span of the rays, in 
%                    degrees. I think this should be constant, even if the
%                    the source moves (and hence the geometry parameters
%                    change). We will make this constant based on the 
%                    maximum value of Rtrue.
% budget           = Number of functions calls that the optimization
%                    function is allowed to call before stopping. Matlab
%                    uses 100 * 2m by default.
% func_delt        = A stopping tolerance for the optimization function.
%                    When the change in function value reaches this
%                    tolerance the optimization stops.
% isImfil          = If true uses Imfil, if false uses lsqnonlin.
% optIter          = The numer of times the BCD loop will run.
%
% The next two are Anderson acceloration parameters.
%
% max_stored_residuals = Determines the maximum amount of previous
%                        iterations that can be used to calculate a new 
%                        solution in the anderson acceleration.
%                         
% dropTol               = If the matrix of residuals in the anderson
%                       acceleration has a condition number beyond dropTol, 
%                       a column is droped until the condition number is 
%                       below drop tol.


%These are the basic variables to change for most runs. The advanced
%parameters are below.
Rnoise          = 1;
ang_noise       = 1;
optIter         = 20;
isImfil         = true;
% Use the syntax below if you want to use your own image. Change
% 'spine.tif' to the image you want to use.
 image           = imresize(double(imread('spine.tif')),[256,256]);
%image           = 'sheppLogan';
runDescription   = "This is the standard test using Irons-Tuck Acceloration."

% Don't Forget to change the filename at the bottom.



%Advanced Parameters to change.
randomSeed      = 4; rng(randomSeed);
n               = 256;
m               = 180;
Rnoise_guess    = 0;
Rguess          = 2;
RPert           = Rnoise*(rand(1,m) - 0.5);
Rtrue           = Rguess*ones(1,m) + RPert;
angles_guess    = (0:2:358);
ang_noise_guess = 0;
p               = length(angles_guess)/m; 
span            = 2*atand(1/(2*max(Rtrue)-1));
ProbOptions     = PRset('CTtype', 'fancurved', 'span', span,'phantomImage',image);
budget          = 100 * 2 * m;
func_delt       = 1e-15;
R_lower = -0.5 * Rnoise;
R_upper = 0.5 * Rnoise;
angle_lower = -0.5 * ang_noise;
angle_upper = 0.5 * ang_noise;


% Structure containing all parameters for information and reproducibility of
% the run.
runInputs = struct('randomSeed',randomSeed,'n',n,'m',m,'Rnoise',Rnoise,...
    'Rnoise_guess',Rnoise_guess,'Rguess',Rguess,'Rtrue',Rtrue,...
    'angles_guess',angles_guess,'ang_noise_guess',ang_noise_guess,'ang_noise',...
    ang_noise,'p',p,'span',span,'ProbOptions',ProbOptions,'budget',budget,'func_delt',...
    func_delt,'optIter',optIter,'R_LOWER',R_lower,...
    'R_upper',R_upper,'angle_lower',angle_lower,'angle_upper',angle_upper,...
    'runDescription',runDescription);

%
% From this point on you probably don't need to understand exactly what the
% code is doing unless you are tying to debug it, just know it runs the bcd
% loop for the specified number of iterations.
%
% There are a couple of variables that may be useful in you investigation.
%      xs       : Contains the solution vector x at each iteration as 
%                 a column vector.
%      x_k      : The final solution vector.
%      paramTrue: A vector with the true R parameter followed by the true
%                 theta parameters.
%      p_0      : The final parameter estimation vector
%      xErrors  : The relative error of the x vector at each iteration
%      pErrors  : The relative error of the parameter vector at each
%                 iteration.

% Check to make sure p = is an integer  
if p ~= fix(p)
    error('p = Nangles/m needs to be an integer')
end
%
% if p is an integer, reshape true and guess vector angles into an array 
% with p rows and m columns, each column corresponds to an entry in vector 
% R. Additionally add the noise to each column of the true angles.
%
angles_true = angles_guess;
angle_pert = ang_noise * (rand(1,m) - 0.5);
angles_true  = reshape(angles_true, p, m) + angle_pert;
angles_guess = reshape(angles_guess,p,m);

[Atrue, btrue, xtrue, ProbInfo] = PRtomo_var(n, Rtrue, angles_true, ProbOptions);
b = PRnoise(btrue);
paramTrue = [RPert angle_pert];

disp("The size of A is");
disp(size(Atrue));
%
% Now setup a similar problem, using the guess for R and guess angles
% (note that using PRtomo_var with a scalar input for R and corresponding
% vector for angles, should produce the same as calling PRtomo)
%
[A, ~, ~, ~] = PRtomo_var(n, Rguess, angles_guess(:), ProbOptions);

%
% Now solve using IRhybrid_lsqr to see the difference in solutions
% using the true A and the guess for A
%

IRoptions   = IRset('Iterbar','off');
[x1, info1] = IRhybrid_lsqr(Atrue, b,IRoptions);
[x0, info2] = IRhybrid_lsqr(A, b,IRoptions);


%
%Now we enter into the BCD loop. First we initialize our guess for the
%parameters on R and theta.
%
RParams = ones(1,m) * Rnoise_guess;
angleParams = ones(1,m) * ang_noise_guess;
% Here we collect the error norm for plotting later.
xErrors = [norm(x0 - xtrue)/norm(xtrue)];
pErrors = [norm(paramTrue - [RParams angleParams])/norm(paramTrue)];
RErrors = [norm(RPert - RParams) / norm(RPert)];
angErrors = [norm(angleParams - angle_pert)/norm(angle_pert)];

% This saves the x_k solutions incase the BCD only shows semi-convergence.
xs = x0;


%Here we set the function options depending on which optimization function
%was decided on being used.
if isImfil
    imOptions = imfil_optset('least_squares',1,'simple_function',1, ...
    'function_delta', func_delt);
    bounds = [ones(1,m) * R_lower ones(1,m)* angle_lower; ...
    ones(1,m) * R_upper ones(1,m) * angle_upper]';
    BCDinfo.bounds = bounds; BCDinfo.imOptions = imOptions;
else %Set up the parameters for the lsqnonlin
    optOptions = optimoptions('lsqnonlin','MaxFunctionEvaluations',budget,...
        'FunctionTolerance',func_delt, 'UseParallel',false);
    lb = [ones(1,m) * R_lower ones(1,m) * angle_lower];
    ub = [ones(1,m) * R_upper ones(1,m) * angle_upper];
    BCDinfo.lb = lb; BCDinfo.ub = ub; BCDinfo.optOptions = optOptions
end


%Intialize the info for the Block Coordinate Descent. 
BCDinfo.n = n; BCDinfo.RParams = RParams; BCDinfo.angleParams = angleParams;
BCDinfo.angles_guess = angles_guess; BCDinfo.ProbOptions = ProbOptions; 
BCDinfo.b = b; BCDinfo.x = x0; BCDinfo.Rguess = Rguess;
BCDinfo.isImfil = isImfil; BCDinfo.budget = budget; BCDinfo.m = m;
BCDinfo.IRoptions = IRoptions;

%Intialization for Crossed Secant
x_curr = x0;
BCDinfo.x = x_curr;
%Initialize old parameters to begin Crossed Secant loop
[G_old, iterInfo] = fpBCD(BCDinfo);
delta_x_old = G_old - x_curr;

%Update the current x guess in the BCD info.
BCDinfo.RParams = iterInfo.RParams;
BCDinfo.angleParams = iterInfo.angleParams;
%Store the errors for plotting at each iteration.
xs = [xs, x_curr];
xErrors = [xErrors,norm(x_curr - xtrue) / norm(xtrue)];
pErrors = [pErrors,norm(paramTrue - iterInfo.p_0)/norm(paramTrue)];
RErrors = [RErrors,norm(RPert - iterInfo.RParams) / norm(RPert)];
angErrors = [angErrors,norm(iterInfo.angleParams - angle_pert)/norm(angle_pert)];

%This enters the AABCD optimization loop.
for k = 1:optIter - 1
    [G_curr, iterInfo] = fpBCD(BCDinfo);
    delta_x_curr = G_curr - x_curr;
    %Change the BCD info for calculating composed G
    %The second difference quoient of x
    dG = G_curr - G_old;
    d_delta_x = delta_x_curr - delta_x_old;
    %The better version
    x_curr = G_curr - dot(dG,d_delta_x)/dot(d_delta_x,d_delta_x) * delta_x_curr;
    %Save the current G and delta_x values for later.
    G_old = G_curr;
    delta_x_old = delta_x_curr;
    %Update the current x guess in the BCD info.
    BCDinfo.x = x_curr;
    BCDinfo.RParams = iterInfo.RParams;
    BCDinfo.angleParams = iterInfo.angleParams;
    %Store the errors for plotting at each iteration.
    xs = [xs, x_curr];
    xErrors = [xErrors,norm(x_curr - xtrue) / norm(xtrue)];
    pErrors = [pErrors,norm(paramTrue - iterInfo.p_0)/norm(paramTrue)];
    RErrors = [RErrors,norm(RPert - iterInfo.RParams) / norm(RPert)];
    angErrors = [angErrors,norm(iterInfo.angleParams - angle_pert)/norm(angle_pert)];
end    

runTime = toc
p_0 = iterInfo.p_0;
x2 = x0;
x_k = x_curr;
%The first aurgument is the name of the file.
save CScontrola runInputs RParams angleParams runTime p_0 x_k ...
    xErrors pErrors RErrors angErrors xs angle_pert RPert paramTrue x1 x2 ...
    xtrue isImfil ProbInfo