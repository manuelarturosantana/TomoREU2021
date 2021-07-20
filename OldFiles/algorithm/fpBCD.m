function [x_k, iterInfo] = fpBCD(inputs)
    %This function computes one iteration of the BCD loop, as if it was a
    %fixed point iteration. 
    % inputs: A structure containing the following information.
    %      n          : The size of one side of the image.
    %      RParams    : The parameters representing the perturbations on R.
    %      angleParams: The Parameters representing the current
    %      perturbations on the angles.
    %      angles_guess: The starting angles which the parameters are added
    %      to.
    %      bounds      : array for imfil bounds.
    %      Proboptions : Info for PRtomo_var
    %      imOptions   : Input options for imfil
    %      b           : right hand side vector.
    %      x           : current guess for x iteration.
    %      Rguess      : starting value for R to which perturbations are added.
    %      lb,ub       : Lower and upper bounds for lsqnonlin
    %      m           : Number of Rparameters.
    % Outputs:
    %      x       : Next x iteration vector.
    %      IterInfo: Structure containing the following fields.
    %            RParams: New R Parameters.
    %            R_vals: New R Values
    %            p_0   : New Parameter vectors.
    %            angleParams: New angle Parameters.

    %This part evaluates g_curr = g(x)
    %Here we perform the non-linear least squares solution
    if inputs.isImfil
        p_0 = optParamParallel_var(inputs.n,inputs.RParams,inputs.angleParams,inputs.angles_guess,...
            inputs.bounds, inputs.budget,inputs.ProbOptions,inputs.imOptions,...
            inputs.b,inputs.x,inputs.Rguess);
    else
        p_0 = optParamParallel(inputs.n,inputs.RParams,inputs.angleParams,inputs.angles_guess,...
            inputs.lb, inputs.ub,inputs.ProbOptions,inputs.optOptions,inputs.b,...
            inputs.x,inputs.Rguess);
    end
    %Here we split in to the optimized solution to then build a better A
    %matrix
    RParams = p_0(1:length(p_0) / 2);
    Rvals = ones(1,inputs.m) * inputs.Rguess + RParams;
    angleParams = p_0((length(p_0) / 2) + 1:end);
    Theta_k = inputs.angles_guess + angleParams;
    [A3,~,~,~] = PRtomo_var(inputs.n,Rvals,Theta_k,inputs.ProbOptions);
    %After building A we minimize in the x block coordinate.
    [x_k, ~] = IRhybrid_lsqr(A3,inputs.b);
    
    %Out put the info for iteration.
    iterInfo.RParams = RParams;
    iterInfo.Rvals = Rvals;
    iterInfo.p_0 = p_0;
    iterInfo.angleParams = angleParams;
   
end