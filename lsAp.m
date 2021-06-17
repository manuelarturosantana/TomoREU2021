function p = lsAp(n,Rvals,thetaVals,inOptions,b,xk)
%Nesting function for performing the minimization. of A(p)x_k - b.
%
% Input:
%      n: problem size. The size of the image should be n x n.
%      Rvals: Array of parameters to representing R parameters
%      thetaVals: Array of parameters representing noise to add to each
%      angle theta
%      inOptions: struct containg the following fields angles, image
%      %%UPDATE THIS AS TIME GOES ON.
%      b: The RHS vector in A(p)x = b
%      xk: The current approximation of the x vector.
    
    param = [Rvals,thetaVals];
    
    p = lsqnonlin(@multAndSub,param);
    %This is a nested function so the pass in parameters are avaliable to 
    % it, but it only takes on parameter for the lsnonlin.
    function vec = multAndSub(x)
        %Since lsqnonlin only takes a single vector we stack the R's and 
        %theta noise on top of each other on the x. 
        Rs = x(1:length(x) / 2); 
        thetas = x((length(x) / 2) + 1:end);
        Ap = makeAp(n,Rs,thetas,inOptions);
        vec = Ap * xk - b;  
    end

end

