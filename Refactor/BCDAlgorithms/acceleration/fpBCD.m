function [x_curr,p] = fpBCD(options,probInfo,RParams,angleParams,optOptions,b,x_curr,bounds,lsSolver)
%{ 
Performs one iteration of the block coordinate descent as if it was a fixed
point iteration.

Inputs:
    options    : Structure created by IRset.
    probInfo   : Structure created by PRset_var.
    RParams    : Current guess for the perturbations on R.
    angleParams: Current guess for the perturbations on the angles.
    optOptions : Structure for options on imfil or lsqnonlin.
    b          : RHS vector also called the sinogram.
    x_curr     : Best approximation for the image vector currently.
    bounds     : numPerts x 2 array containing the lower bounds in the first
                 column and the upper bounds in the second column.
    lsSolver   : Function handle to solve the least square problem.

Outputs:
    x_curr     : The new image vector approximation.\
    p          : The new parameter vector approximation of the form
                [Rparams angleParams]
%} 

%Here we perform the non-linear least squares solution
    if strcmp(options.nonlinSolver,'imfil')
        p = optParamParallel_var(RParams,angleParams,probInfo,options.budget,...
            optOptions,b,x_curr,bounds);
    else
        %We default to imfil style bounds passed in, so here we split into
        %the correct style for lsqnonlin
        lb = bounds(:,1);
        ub = bounds(:,1);
        p = optParamParallel(RParams,angleParams,probInfo,optOptions,b,x_curr,lb,ub);
    end
    %Here we split in to the optimized solution to then build a better A
    %matrix
    RParams = p(1:length(p) / 2);
    Rvals = probInfo.TomoInfo.Rvar + RParams;
    angleParams = p((length(p) / 2) + 1:end);
    Theta_k = probInfo.anglesvar + angleParams;
    A = createA(probInfo.n,Rvals,Theta_k,probInfo.TomoInfo);
    %After building A we minimize in the x block coordinate.
    [x_curr, ~] = lsSolver(A,b,options);
end