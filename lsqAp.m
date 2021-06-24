function p = lsqAp(n,Rparam,thetaParam,angles,lb,ub,PRoptions,optOptions,b,xk,Rstart)
%Nesting function for performing the minimization. of A(p)x_k - b
% using the matlab optimization toolbox function lsqnonlin.
%
% Input:
%      n         : Problem size. The size of the image should be n x n.
%      RParam    : Array of parameters to representing R parameters
%      thetaParam: Array of parameters representing noise to add to each
%                  angle theta
%      angles    : Assumed angles for tomography problem. Should be in the
%                  same format as PRtomo_var.
%      lb        : Same lowerbound as lsqnonlin.
%      ub        : Same upperbound as lsqnonlin.
%      PRoptions : Structure used in PRtomo
%      optOptions: Structure used for options in lsqnonlin  
%      b         : The RHS vector in A(p)x = b
%      xk        : The current approximation of the x vector.
%      Rstart    : The value which to add the Rparamters to.
%
%      Output:
%      p: A row vector such that the first half of the elements and the
%      optimized R parameters, and the second half are the optimized theta
%      parameters.
    param = [Rparam,thetaParam];
    
    p = lsqnonlin(@multAndSub,param,lb,ub,optOptions);
    %This is a nested function so the pass in parameters are avaliable to 
    % it, but it only takes on parameter for the lsnonlin.
    function vec = multAndSub(x)
        %Since lsqnonlin only takes a single vector we stack the R's and 
        %theta noise on top of each other on the x. 
        Rparams = x(1:length(x) / 2);
	    Rvals = Rstart * ones(1,length(Rparams)) + Rparams;
        thetaParams = x((length(x) / 2) + 1:end);
        %To each column of angles this add the corresponding theta
        %parameter. This assuming that thetaParams is a row vector.
        angles2 = angles + thetaParams;
        Ap = PRtomo_var(n,Rvals,angles2,PRoptions);
        vec = Ap * xk - b;  
    end
end
