function paramVec = optParamParallel(n,RParams,angleParams,angles,lb,ub,PRoptions,optOptions,b,xk)
% optParamParallel
%   Inputs: 
%         n           : The size of one part of the image being tested.
%         RParams     : Current guess for the R parameters.
%         angleParams : Current guess for the angle parameters.
%         angles      : The initial guess angles. Should be in a form such that
%                     there are length(RParams) columns, and each column 
%                     corresponds to the the angles for that parameter.
%         lb          : Lower bound vector for the optimization function.
%         ub          : Upper bound vector for the optimization function.
%         PRoptions   : Options for PRtomo
%         optOptions  : Options for lsqnonlin
%         b           : Right hand side vector
%         xk          : Current guess for the x vector.
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
    
    % This reshape makes two columns each with the first m items.
    lb = reshape(lb,m,2);
    ub = reshape(ub,m,2);
    parfor i = 1:m
        x = lsqAp(n,RParams(i),angleParams(i),angles(:,i),lb(i,:),ub(i,:),...
            PRoptions,optOptions,b,xk);
        Rvals(i) = x(1);
        angleVals(i) = x(2);
    end
    paramVec = [Rvals,angleVals];
end

