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
                rng(5);
n                = 64;
m                = 4;
noise_scale      = 0.5;
noise_guess      = 0;
perturbations    = noise_scale * (rand(1,m) - 0.5);
Rguess           = 2;
angles_guess     = (0:2:358);
p                = length(angles_guess)/m; 
ProbOptions      = PRset('CTtype', 'fancurved','phantomImage','sheppLogan');
optIter          = 10;
lb              = -0.5;
ub              = 0.5;
isR             = true;


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

angles_guess = reshape(angles_guess,p,m);
angles_true = angles_guess;
Rtrue = Rguess * ones(1,m);
if isR
    Rtrue = Rtrue + perturbations;
else
    angles_true  = angles_true + perturbations;
end
span             = 2*atand(1/(2*max(Rtrue)-1));
ProbOptions.span = span;
%
% if p is an integer, reshape true and guess vector angles into an array 
% with p rows and m columns, each column corresponds to an entry in vector 
% R. Additionally add the noise to each column of the true angles.
%


[Atrue, btrue, xtrue, ProbInfo] = PRtomo_var(n, Rtrue, angles_true, ProbOptions);
b = PRnoise(btrue);

paramTrue = perturbations;

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


% Now we enter the BCD loop using x2 as our first x calculation.
x_k = x2;
% Here we collect the error norm for plotting later.
xErrors = [norm(x2 - xtrue)/norm(xtrue)];
pErrors = [abs(paramTrue - noise_guess)/abs(paramTrue)];

% This saves the x_k solutions incase the BCD only shows semi-convergence.
xs = [x_k];


%This enters the BCD optimization loop.
for i = 2:optIter
    %Here we perform the non-linear least squares solution
    p_0 = singOptParallel(n,m,angles_guess,lb,ub,ProbOptions,b,x_k,Rguess,isR);
    if isR
        Rvals = ones(1,m) * Rguess + p_0;
        Theta_k = angles_true;
    else 
        Rvals = Rtrue;
        Theta_k = angles_guess + p_0;
    end
    [A3,~,~,~] = PRtomo_var(n,Rvals,Theta_k,ProbOptions);
    %After building A we minimize in the x block coordinate.
    [x_k, info_k] = IRhybrid_lsqr(A3,b);
    xs = [xs,x_k];
    xErrors = [xErrors,norm(x_k - xtrue) / norm(xtrue)];
    pErrors = [pErrors,abs(paramTrue - p_0)/abs(paramTrue)];
    disp(i)
    disp(p_0)
end

figure(4), clf
PRshowx(x_k,ProbInfo)
if isR
    title('Solution After BCD on R','fontsize', 20)
else 
    title('Solution After BCD on angles','fontsize', 20)
end

figure(5), clf
plot(xErrors,'-o');
hold on 
plot(pErrors,'-*');
hold off
legend('xError Norms','p Error Norms');
xlabel('Number of Iterations','fontsize',15);
ylabel('Relative error','fontsize',15);