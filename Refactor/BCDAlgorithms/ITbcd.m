function [x,iterInfo] = ITbcd(b,options,probInfo)
% Performs the block coordinate descent without any acceleration technique
%  [x,iterInfo] = NAbcd(b,options,probInfo)
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
angleDiff = probInfo.anglesvar(2,1) - probInfo.anglesvar(1,1);
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

%This enters the BCD optimization loop.
for i = 2:options.BCDmaxIter
    x_old = x_curr;
    [G, p] = fpBCD(options,probInfo,RParams,angleParams,optOptions,b,x_curr,bounds,lsSolver);
    delta_x = G - x_curr;
    %This is G(G(x_n)) now x_curr passed in is G
    [GoG,~] = fpBCD(options,probInfo,RParams,angleParams,optOptions,b,G,bounds,lsSolver);
    delta_G = GoG - G;
    %The second difference quoient of x
    delta_x2 = delta_G - delta_x;
    %The accelerated x
    x_curr = GoG - (dot(delta_G, delta_x2)/norm(delta_x2)^2) * delta_G;
    
    RParams = p(1:length(p) / 2);
    angleParams = p((length(p) / 2) + 1:end);
    xs = [xs,x_curr];
    
    xErrors = [xErrors norm(x_curr - xtrue)/norm(xtrue)];
    pErrors = [pErrors norm(paramTrue - [RParams angleParams])/norm(paramTrue)];
    RErrors = [RErrors norm(probInfo.true.Rpert - RParams) / norm(probInfo.true.Rpert)];
    angErrors = [angErrors norm(angleParams - probInfo.true.anglePert)/norm(probInfo.true.anglePert)];
    
    if strcmp(options.dispIter,'on')
        fprintf('BCD iteration %d completed \n',i)
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
    iterInfo.angErrors = angErrors;
    iterInfo.numIter   = i;
    iterInfo.runTime = runTime;
    iterInfo.xsols = xs;
end