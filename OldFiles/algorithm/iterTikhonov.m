function [x,iterCount] = iterTikhonov(A,b,q,tol,maxIter)
    % Function to compute an Iterated Tikhonov regularization
    % Input:
    %      A  : Matrix in the damped least squares equation
    %      b  : b vector in the dampled least squares equatoin
    %      q  : Set to 1 if you want a stationary iterated Tikhonov
    %          regularization, or set to a number between 0 and 1 for a
    %          non-stationary version with a geomtric sequence.
    %      tol: Tolerance to stop the iteration when the change in x is
    %          below this tolerance.
    %      maxIter: The maximum number of iterations to run the
    %      regularizatoin.
    %
    % Output:  
    %      x         : The vector after the regularization.
    %      iterCount : The number of iterations actually used.
    
    %Generate an initial x vector and best regularization parameter.
    [x_0, info] = IRhybrid_lsqr(A,b);
    alpha = info.StopReg.RegP;
    
    %Intialize error and iterations
    iterCount = 1;
    error = 10 * tol;
    x_old = x_0;
    
    %Set up the normal equations.
    A_trans_b = A' * b;
    ATA = A' * A;
    [n,~] = size(ATA);
    while error > tol && iterCount < maxIter
        alpha = alpha * q;
        %Finishing normal equations set up for Tikhonov Iteration
        b_n = A_trans_b + alpha * x_old;
        B = ATA + alpha * speye(n);
        %Matlab built in Conjugate Gradient method with Preconditioners.
        x_new = pcg(B,b_n);
        error = norm(x_new - x_old);
        x_old = x_new;
        iterCount = iterCount + 1;
    end
    x = x_old;
end

