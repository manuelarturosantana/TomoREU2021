function p = lsqAp_var(n,Rvals,thetaVals,angles,bounds,budget,PRoptions,imOptions,b,xk)
% Nesting function for performing the minimization. of A(p)x_k - b. Performs
% the optimization using the imfil function and package. 
%
% Input:
%      n: problem size. The size of the image should be n x n.
%      Rvals: Array of parameters to representing R parameters
%      thetaVals: Array of parameters representing noise to add to each
%      angle theta
%      bounds: a N x 2 array where the first column is the lower bounds on
%      the constraints and the second is the upper. Rlower shouldn't be
%      lower than half square root of 2.
%      budget: Same as budget for imfill, assuming simple_function=1
%      PRoptions: Options for PRtomo
% Output: 
%      p: a column vector with the first half containing the minimization
%      parameters for R, and the second half containing those for theta.
    
    param = [Rvals,thetaVals]';
    
    p = imfil(param,@multAndSub,budget,bounds,imOptions);
    p = p';
    %This is a nested function so the pass in parameters are avaliable to 
    % it, making it easier to use imfil.
    function vec = multAndSub(x)
        %Since imfil only takes a single vector we stack the R's and 
        %theta noise on top of each other on the x. 
        Rparams = x(1:length(x) / 2); 
        thetaParams = x((length(x) / 2) + 1:end);
        %To each column of angles this add the corresponding theta
        %parameter.
        angles2 = angles + thetaParams';
        Ap = PRtomo_var(n,Rparams,angles2,PRoptions);
        vec = Ap * xk - b;
    end
end

