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
                rng(1);
n               = 64;
m               = 4;   
Rnoise          = 1;
Rguess          = 2;   
Rtrue           = Rguess*ones(1,m) + Rnoise*(rand(1,m) - 0.5);
angles_guess    = (0:2:358);
ang_noise_guess = 0;
ang_noise       = 1;
p               = length(angles_guess)/m; 
span            = 2*atand(1/(2*max(Rtrue)-1));
ProbOptions     = PRset('CTtype', 'fancurved', 'span', span,'phantomImage','sheppLogan');
budget          = 100 * 2 * m;
func_delt       = 1e-6;
optIter         = 3;
isImfil         = false;

%
% Here we set the bounds for the optimization function as a column vector. R_LOWER
% is set as a constant, as 0.5 * sqrt(2) is the lowest value that PR tomo
% will allow for R.
%

R_LOWER = 0.5 * sqrt(2) + .001;
R_upper = 2.5;
angle_lower = -.5;
angle_upper = .5;


%
% From this point on you probably don't need to understand exactly what the
% code is doing unless you are tying to debug it, just know it runs the bcd
% loop for the specified number of iterations.
%

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
paramTrue = [Rtrue angle_pert];

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
[x1, info1] = IRhybrid_lsqr(Atrue, b);
[x2, info2] = IRhybrid_lsqr(A, b);

figure(1), clf
PRshowx(xtrue, ProbInfo)
title('True Solution','fontsize', 20)

figure(2), clf
PRshowx(x1, ProbInfo)
title('Solution with True A','fontsize', 20)

figure(3), clf
PRshowx(x2, ProbInfo)
title('Solution with Noisey A','fontsize', 20)

%
%Now we enter into the BCD loop. First we initialize our guess for the
%parameters on R and theta.
%
RParams = ones(1,4) * Rguess;
angleParams = ones(1,4) * ang_noise_guess;
%We use x2 as our first x calculation.
x_k = x2;
% Here we collect the error norm for plotting later.
xErrors = [norm(x2 - xtrue)/norm(xtrue)];
pErrors = [norm(paramTrue - [RParams angleParams])/norm(paramTrue)];

% This saves the x_k solutions incase the BCD only shows semi-convergence.
xs = [x_k];


%Here we set the function options depending on which optimization function
%was decided on being used.
if isImfil
    imOptions = imfil_optset('least_squares',1,'simple_function',1, ...
    'function_delta', func_delt);
    bounds = [ones(1,m) * R_LOWER ones(1,m)* angle_lower; ...
    ones(1,m) * R_upper ones(1,m) * angle_upper]';
else %Set up the parameters for the lsqnonlin
    optOptions = optimoptions('lsqnonlin','MaxFunctionEvaluations',budget,...
        'FunctionTolerance',func_delt, 'UseParallel',true);
    lb = [ones(1,m) * R_LOWER ones(1,m) * angle_lower];
    ub = [ones(1,m) * R_upper ones(1,m) * angle_upper];
    
end

%This enters the BCD optimization loop.
for i = 2:optIter
    %Here we perform the non-linear least squares solution
    if isImfil
        p_0 = lsqAp_var(n,RParams,angleParams,angles_guess,bounds, budget,...
        ProbOptions,imOptions,b,x_k);
    else
        p_0 = lsqAp(n,RParams,angleParams,angles_guess,lb,ub,ProbOptions,...
            optOptions,b,x_k);
    end
    %Here we split in to the optimized solution to then build a better A
    %matrix
    RParams = p_0(1:length(p_0) / 2);
    angleParams = p_0((length(p_0) / 2) + 1:end);
    Theta_k = angles_true + angleParams;
    [A3,~,~,~] = PRtomo_var(n,RParams,Theta_k,ProbOptions);
    %After building A we minimize in the x block coordinate.
    [x_k, info_k] = IRhybrid_lsqr(A3,b);
    xs = [xs,x_k];
    xErrors = [xErrors,norm(x_k - xtrue) / norm(xtrue)];
    pErrors = [pErrors,norm(paramTrue - p_0)/norm(paramTrue)];
end

figure(4), clf
PRshowx(x_k,ProbInfo)
if isImfil
    title('Solution After BCD (Imfil)','fontsize', 20)
else 
    title('Solution After BCD (lsqnonlin)','fontsize', 20)
end

figure(5), clf
plot(xErrors)
hold on 
plot(pErrors)
hold off
legend('xError Norms','p Error Norms');