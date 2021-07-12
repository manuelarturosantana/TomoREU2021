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

                rng(5);
n               = 64;
m               = 4;   
Rnoise          = 0.5;
Rnoise_guess    = 0;
Rguess          = 2;
RPert           = Rnoise*(rand(1,m) - 0.5);
Rtrue           = Rguess*ones(1,m) + RPert;
angles_guess    = (0:2:358);
ang_noise_guess = 0;
ang_noise       = 0.5;
p               = length(angles_guess)/m; 
span            = 2*atand(1/(2*max(Rtrue)-1));
image           = imresize(double(imread('spine.tif')),[256 256]);
ProbOptions     = PRset('CTtype', 'fancurved', 'span', span,'phantomImage',image);
budget          = 100 * 2 * m;
func_delt       = 1e-6;
optIter         = 10;
isImfil         = true;
max_stored_residuals = 3;
dropTol = 1e6;

%
% Here set the bounds for the optimization function to constrain R and the
% angle perterbation to.
%

R_lower = -0.5;
R_upper = 0.5;
angle_lower = -0.5;
angle_upper = 0.5;


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
[x0, info2] = IRhybrid_lsqr(A, b);

figure(1), clf
PRshowx(xtrue, ProbInfo)
title('True Solution','fontsize', 20)

figure(2), clf
PRshowx(x1, ProbInfo)
title('Solution with True A','fontsize', 20)

figure(3), clf
PRshowx(x0, ProbInfo)
title('Solution with Noisey A','fontsize', 20)

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


%Intialization for Anderson Acceloration.
x_curr = x0;
BCDinfo.x = x_curr;
G = [ ];
num_stored_residuals = 0;
%Number of previous solutions to use with Anderson Acceloration


%This enters the AABCD optimization loop.
for k = 1:optIter
    [g_curr, iterInfo] = fpBCD(BCDinfo);
    f_curr = g_curr - x_curr;
    if k > 1
        %Form G on the second iteration.
        delta_f = f_curr - f_old;
        delta_g = g_curr - g_old;
   
        %Update the columns on G
        if num_stored_residuals < max_stored_residuals
            G = [G delta_g];
        else
            G = [G(:,2:end) delta_g];
        end
        %We update this even though the second step doesn't actually
        %increase the number of columns to G. This is a work saving step
        %because we have to redo the QR factorization of F, but will need G
        %in this size later.
        num_stored_residuals = num_stored_residuals + 1;
    end
    %Set f and g old for calculating delta f the next iteration.
    f_old = f_curr; g_old = g_curr;
    %Sets the x value on the first iteration
    if num_stored_residuals == 0
        x_curr = g_curr;
    else
        %one vector QR factorization
        if num_stored_residuals == 1 
            Q(:,1) = delta_f / norm(delta_f);
            R(1,1) = norm(delta_f);
        else
            %delete a QR column if too many.
            if num_stored_residuals > max_stored_residuals 
                [Q,R] = qrdelete(Q,R,1);
                num_stored_residuals = num_stored_residuals - 1;
                %This deals with a qrdelete usage explained below.
                if size(R,1) ~= size(R,2)
                    Q = Q(:,1:num_stored_residuals -1);
                    R = R(1:num_stored_residuals - 1,:);
                end
                %Explination: If Q is not square then the matlab function
                %QR delete removes one column of Q, and and column and row
                %of R. If Q is square then the column demenion of Q is not
                %reduced, and R only has a column removed. This behavior is
                %to account for thick QR decompoitions, but since we are 
                %using a thin one we must correct if this happens.
            end
            %One pass of modified Gram-Schmidt to update the QR
            %factorization. Recall this uses the fact that since F = QR
            %then the column of F is a linear combination of columns of Q
            %with coefficents coming from the last column of R.
            for i = 1:num_stored_residuals -1 
                R(i,num_stored_residuals) = Q(:,i)' * delta_f;
                delta_f = delta_f - R(i,num_stored_residuals) * Q(:,i); 
            end
            %Completing the Graham-Scmidt Iteration
            R(num_stored_residuals,num_stored_residuals) = norm(delta_f);
            Q = [Q,delta_f / norm(delta_f)]; 
        end
        %Here we delete more columns in the QR factorization to deal with
        %poor conditioning of the R matrix that may occur.
        while (cond(R)) > dropTol && num_stored_residuals > 1
           [Q,R] = qrdelete(Q,R,1);
           num_stored_residuals = num_stored_residuals - 1;
           %If we change the size of the QR, we must also change the size
           % of G.
           G = G(2:end);
           %See above explination for this step
           if size(R,1) ~= size(R,2)
                    Q = Q(:,1:num_stored_residuals -1);
                    R = R(1:num_stored_residuals - 1,:);
           end
        end
       gamma = R \ (Q' * f_curr);
       %updating the x info based on the anderson acceleration.
       x_curr = g_curr - G * gamma;
    end
    
%Update the current x guess in the BCD info.
BCDinfo.x = x_curr;
%Store the errors for plotting at each iteration.
xs = [xs, x_curr];
xErrors = [xErrors,norm(x_curr - xtrue) / norm(xtrue)];
pErrors = [pErrors,norm(paramTrue - iterInfo.p_0)/norm(paramTrue)];
RErrors = [RErrors,norm(RPert - iterInfo.RParams) / norm(RPert)];
angErrors = [angErrors,norm(iterInfo.angleParams - angle_pert)/norm(angle_pert)];
end    

    
    
figure(4), clf
PRshowx(x_curr,ProbInfo)
if isImfil
    title('Solution After BCD (Imfil)','fontsize', 20)
else 
    title('Solution After BCD (lsqnonlin)','fontsize', 20)
end

figure(5), clf
plot(xErrors,'-o','LineWidth',2);
hold on 
plot(pErrors,'-*','LineWidth',2);
plot(RErrors,'-d','LineWidth',2);
plot(angErrors,'-^','LineWidth',2);
hold off
legend('xError Norms','p Error Norms', 'R error Norms','Angle Error Norm');
xlabel('Number of Iterations','fontsize',15);
ylabel('Relative error','fontsize',15);