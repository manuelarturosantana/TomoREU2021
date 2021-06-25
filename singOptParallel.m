function paramVec = singOptParallel(n,params,angles,lb,ub,PRoptions,b,xk,Rstart,isR)
% singOptParallel Wrapper function to perform the nonlinear optimization
% problem in parallel for either R or Theta
%   Inputs: 
%         n           : The size of one part of the image being tested.
%         params      : initial guess on the noise.
%         angles      : The initial guess angles. Should be in a form such that
%                       there are length(RParams) columns, and each column 
%                       corresponds to the the angles for that parameter.
%         lb          : Lower bound for parameter.
%         ub          : Upper bound for parameter.
%         PRoptions   : Options for PRtomo
%         b           : Right hand side vector
%         xk          : Current guess for the x vector.
%         Rstart      : R value to add R params to.
%         isR         : True if R is the parameter we are minimizing, false
%                       if it is not.
%
%    Output:
%         paramVec    : A vector wit the optimized R parameters in the first
%                       half, and angle parameters in the second.
    
    m = length(params);
    %
    %Here we reshape b so that each column contains the correct values of b
    %for each sub problem. The syntax says we want length(b) / m  equally
    % sized rows and matlab will then fill in the number of columns
    %
    b = reshape(b,length(b) / m,[]);
    paramVec = zeros(1,m);
    
    %This reassignment is to prevent broadcasting in the parfor loop.
    parfor i = 1:m
        paramVec(i) = lsaSing(n,param,angles,lb,ub,PRoptions,b(:,i,...
            xk,Rstart,isR);
    end
end
