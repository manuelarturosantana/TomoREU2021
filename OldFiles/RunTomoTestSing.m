%
% This script sets up test problem with changes in R or Theta.
% To run different versions of the problem you just need to change the values
% in the first few lines. Here are some comments about the variables:
%  m               = number of times R changes 
% noise_scale      = the value to scale the amount of uniform random noise
%                    from the interval [-0.5,0.5];
% noise_guess      = guess on the intial amount of noise. Can be a scalar
%                    or 1 X m vector. It is only used for the inital
%                    parameter error.
% perturbations    = The the actual perturbations.
% Rguess           = The inital guess for R. Used as true R if not testing
%                    R.
% angles_guess     = The initial guess angles. Used as the true angles if
%                    not testing the angles.
%  p               = [total number of projection angles]/m, should be an 
%                    integer (see also PRtomo_var, where p is the number of
%                    rows in the array angles)
% optIter          = The number of times the BCD loop will run.
% lb,ub            = Lower and upper bounds on the amount size of
%                    perturbations during the optimization.
% isR              = Set to True if testing R values, false if testing theta
%                    values
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
optIter          = 40;
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
%      paramTrue: A vector with the true perturbations
%      p_0      : The final parameter estimation vector
%      xErrors  : The relative error of the x vector at each iteration
%      pErrors  : The relative error of the parameter vector at each
%                 iteration.

% Check to make sure p = is an integer  
if p ~= fix(p)
    error('p = Nangles/m needs to be an integer')
end


%
% If p is an integer, reshape the angles into the right form for
% PRtomo_var. Then add perturbations to the appropriate parameter to create
% a true value different than the guess.
%
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


% Set up the true values problem.
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


%
% Plot the best x approximation, as well as the error norms.
%

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