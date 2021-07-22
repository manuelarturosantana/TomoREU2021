function [x,iterInfo] = AAbcd(b,options,probInfo)
% Performs the block coordinate descent using Anderson Acceleration
%  [x,iterInfo] = AAbcd(b,options,probInfo)
%  Inputs:
%       b: Sinogram vector from PRtomo_var
%  options: structure from IRset. Contains the usual fields that can be
%          passed into the IRset options, which will then be passed into
%          the linear least squares solver. Additionally it contains these
%          fields
%          IRset new options
%
%     BCDstopTol  - Stopping tolerance for the relative norm change between
%                   iterations in the BCD loop. Different than the one provided
%                   in IRset to allow the user to still set a stopping tol for
%                   the least squares solver function.
%                   [positive scalar | {1e-3}]
% 
%     Rbounds    - Percentage of R to use as a bound above and below.
%                  Lowerbound is max(Rvar) - max(Rvar(i)) * Rbounds, and
%                  upperbound is max(Rvar(i)) + max(Rvar(i)) * Rbounds.
%                  [positive scalar | {0.5}]
% 
%     angleBounds  - Percentage of difference between angles to be used for the
%                  bounds. Assumes angles are equally spaced.
%                  [positve scalar | {0.5}]
% 
%     nonlinSolver - Which solver to use for the non-linear least squares
%                   problem. imfil requires the imfil package and lsqnonlin
%                   requires the matlab optmization toolbox.
%                   ['imfil' | {'lsqnonlin'}]
% 
%     accel        - Determines which fixed point acceleration technique is
%                    used.
%                    ['anderson','ironstuck','secant',{'none'}]
% 
%     BCDmaxIter   - Maximum number of iterations of the BCDloop. Seperate from
%                    the standard max iterations in IRset so the user can use
%                    IRset for the least squares solver too.
%                    [positive integer | {10}]              
% 
%     BCDlsSolver  - Which solver from IR tools to use to solve the
%                    regularization problem. Seperate from standard version in
%                    IRset to allow user to use IR set for least squares solver
%                    too.
%                    [{'lsqr'},'cgls','irn','fista']
%     dispIter     - Display the iteration during each loop of the BCD
% 
%     Anderson Acceleration Only
%          maxRes - The maximum number of residuals allowed in the anderson
%                   acceleration.
%                   [positive integer | {3}]
% 
%          dropTol - Tolerance to drop a column in the QR factorization is the
%                    condition number of the residual matrix is above dropTol.
%                    [positive scalar | {1e6}]
% 
%     Optimization Parameters
%         budget - Maximum number of function evaluations that can be used in
%                  evaluating finding the minimum. 
%                  [positive integer | {2 *100 * number of perturbations}]
% 
%        funcDelt - Value for stopping criterion of size of change in function 
%                   while doing the optimization
%                  [positive scalar | {1e-6}]
%
% probInfo: structure from PRtomo_var
%    
% Outputs
%       x: The image vector after the BCD 
%       iterInfo: Structure containg the following fields
%               xErrors: Relative image error at each iteration.
%               pErrors: Relative parameter errors at each iteration.
%               RErrors: Relative R Errors at each iteration.
%               angErrors: Relative angle Errors at each iteration.
%               numIter : Number of Iterations completed before stopping
%               runTime : Total run time in seconds

%Intialize true parameters for keeping track
numPert = length(probInfo.true.Rpert);
paramTrue = [probInfo.true.Rpert probInfo.true.anglePert];
xtrue     = probInfo.true.x;

% Create the A matrix from the initial parameter guess
A = createA(probInfo.n, probInfo.TomoInfo.Rvar, probInfo.anglesvar, probInfo.TomoInfo);

%Select the Tikhonov Regularization solver 
if strcmp(options.BCDlsSolver,'lsqr')
    lsSolver = @IRhybrid_lsqr;
elseif strcmp(options.BCDlsSolver,'cgls')
    lsSolver = @IRcgls;
elseif strcmp(options.BCDlsSolver,'irn')
    lsSolver = @IRirn;
elseif strcmp(options.BCDlsSolver,'fista')
    lsSolver = @IRfista;
else
    error('BCD least squares solve not recognized')
end

%time it
tic

%Solve for the first iteration of x, passing in options which contains
%IRoptions.
[x_curr, ~] = lsSolver(A, b,options);
if strcmp(options.dispIter,'on')
        fprintf('BCD iteration %d completed \n',1)
end

%Initialize the Rparameters and angle parameters, assuming intial guess of
%no perturbations.
RParams = zeros(1,length(probInfo.TomoInfo.Rvar));
angleParams = zeros(1,length(probInfo.TomoInfo.Rvar));

% Here we collect the error norm for plotting later.
xErrors = norm(x_curr - xtrue)/norm(xtrue);
pErrors = norm(paramTrue - [RParams angleParams])/norm(paramTrue);
RErrors = norm(probInfo.true.Rpert - RParams) / norm(probInfo.true.Rpert);
angErrors = norm(angleParams - probInfo.true.anglePert)/norm(probInfo.true.anglePert);

% This saves the x_k solutions for later plotting
xs = x_curr;


% Calculating the bounds.
Rmax = max(probInfo.TomoInfo.Rvar);
R_lower = - Rmax * options.RBounds; %Optimization is over the perturbations.
R_upper =  Rmax * options.RBounds;

%The next block is to account for the case of perterbations at every angle.
[n,~] = size(probInfo.anglesvar); 
if n == 1
   angleDiff = probInfo.anglesvar(2) - probInfo.anglesvar(1);
else
    angleDiff = probInfo.anglesvar(2,1) - probInfo.anglesvar(1,1);
end

angle_lower = -angleDiff * options.angleBounds;
angle_upper = angleDiff * options.angleBounds;

%Here we set the function options depending on which optimization function
%was decided on being used.
if strcmp(options.nonlinSolver,'imfil')
    optOptions = imfil_optset('least_squares',1,'simple_function',1, ...
    'function_delta', options.funcDelt);
    bounds = [ones(1,numPert) * R_lower ones(1,numPert)* angle_lower; ...
    ones(1,numPert) * R_upper ones(1,numPert) * angle_upper]';
elseif strcmp(options.nonlinSolver,'lsqnonlin') 
    optOptions = optimoptions('lsqnonlin','MaxFunctionEvaluations',options.budget,...
        'FunctionTolerance',options.funcDelt);
    bounds = [ones(1,numPert) * R_lower ones(1,numPert)* angle_lower; ...
    ones(1,numPert) * R_upper ones(1,numPert) * angle_upper]';
else
    error('BCD nonlinear optimization function not recognized')
end

%Initialization for Anderson Acceleration
G = [ ];
num_stored_residuals = 0;

%This enters the AABCD optimization loop.
for k = 2:options.BCDmaxIter 
    x_old = x_curr;

    [g_curr, p] = fpBCD(options,probInfo,RParams,angleParams,optOptions,b,x_curr,bounds,lsSolver);
    f_curr = g_curr - x_curr;
    if k > 2
        %Form G on the second iteration.
        delta_f = f_curr - f_old;
        delta_g = g_curr - g_old;

        %Update the columns on G
        if num_stored_residuals < options.maxRes
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
            if num_stored_residuals > options.maxRes
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
        while (cond(R)) > options.dropTol && num_stored_residuals > 1
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

    %Separate the parameter vector.
    RParams = p(1:length(p) / 2);
    angleParams = p((length(p) / 2) + 1:end);
    
    %Store the errors for plotting at each iteration.
    xs = [xs,x_curr];
    xErrors = [xErrors norm(x_curr - xtrue)/norm(xtrue)];
    pErrors = [pErrors norm(paramTrue - [RParams angleParams])/norm(paramTrue)];
    RErrors = [RErrors norm(probInfo.true.Rpert - RParams) / norm(probInfo.true.Rpert)];
    angErrors = [angErrors norm(angleParams - probInfo.true.anglePert)/norm(probInfo.true.anglePert)];   
    
    if strcmp(options.dispIter,'on')
        fprintf('BCD iteration %d completed \n',k) 
    end
    %Stopping condition
    if norm(x_curr - x_old) / norm(x_old) < options.BCDStopTol
        break
    end
end
    runTime = toc;
    x = x_curr;
    %Store things to return
    iterInfo.xErrors = xErrors;
    iterInfo.pErrors = pErrors;
    iterInfo.RErrors = RErrors;
    iterInfo.angleErrors = angErrors;
    iterInfo.numIter   = k;
    iterInfo.runTime = runTime;
    iterInfo.xsols = xs;
    
    %Create solution with true parameters and return it.
    Rvals =   probInfo.TomoInfo.Rvar + probInfo.true.Rpert;
    angleVals = probInfo.anglesvar + probInfo.true.anglePert; 
    Atrue = createA(probInfo.n,Rvals,angleVals,probInfo.TomoInfo);
   [xtrueparam, ~] = lsSolver(Atrue, b,iterOptions);
   iterInfo.xtrueparam = xtrueparam;
end