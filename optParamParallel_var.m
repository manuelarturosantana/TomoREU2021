function paramVec = optParamParallel_var(n,RParams,angleParams,angles,bounds,budget,PRoptions,optOptions,b,xk,Rstart)
% optParamParallel Wrapper function to perform the nonlinear optimization
% problem in parallel using imfil.
%   Inputs: 
%         n           : The size of one part of the image being tested.
%         RParams     : Current guess for the R parameters.
%         angleParams : Current guess for the angle parameters.
%         angles      : The initial guess angles. Should be in a form such that
%                     there are length(RParams) columns, and each column 
%                     corresponds to the the angles for that parameter.
%         bounds      : Same as bounds for imfil.
%         budget      : Same as budget for imfil.
%         PRoptions   : Options for PRtomo
%         optOptions  : Options for lsqnonlin
%         b           : Right hand side vector
%         xk          : Current guess for the x vector.
%         Rstart      : R value to add R params to.
%
%    Output:
%         paramVec    : A vector wit the optimized R parameters in the first
%                       half, and angle parameters in the second.
    
    m = length(RParams);
    %
    %Here we reshape b so that each column contains the correct values of b
    %for each sub problem. The syntax says we want length(b) / m  equally
    % sized rows and matlab will then fill in the number of columns
    %
    b = reshape(b,length(b) / m,[]);
    Rvals = zeros(1,m);
    angleVals = zeros(1,m);
    
    %This reassignment is to prevent broadcasting in the parfor loop.
    Rbounds = bounds(1:m,:);
    angleBounds = bounds(1 + m:end,:);
    parfor i = 1:m
        probBounds = [Rbounds(i,:); angleBounds(i,:)];
        x = lsqAp_var(n,RParams(i),angleParams(i),angles(:,i),probBounds,...
            budget,PRoptions,optOptions,b(:,i),xk,Rstart);
        Rvals(i) = x(1);
        angleVals(i) = x(2);
    end
    paramVec = [Rvals,angleVals];
end


