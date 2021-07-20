function paramVec = optParamParallel(RParams,angleParams,probInfo,optOptions,b,x_curr,lb,ub)
% optParamParallel: Wrapper function to perform nonlinear least squares
%                   function using lsqnonlin.
%   Inputs: 
%         RParams     : Current guess for the R parameters.
%         angleParams : Current guess for the angle parameters.
%         probInfo    : Structure from PRtomo_var.
%         optOptions  : Options for lsqnonlin
%         b           : Right hand side vector
%         xk          : Current guess for the x vector.
%         lb          : Same as lb in lsqnonlin.
%         ub          : Same as ub in lsqnonlin.
%
%    Output:
%         paramVec    : A vector wit the optimized R parameters in the first
%                       half, and angle parameters in the second.
    
    numPert = length(probInfo.TomoInfo.Rvar);
    %
    %Here we reshape b so that each column contains the correct values of b
    %for each sub problem. The syntax says we want length(b) / m  equally
    % sized rows and matlab will then fill in the number of columns
    %
    b = reshape(b,length(b) / numPert,[]);
    
    %Initializing the parameter vectors
    Rvals = zeros(1,numPert);
    angleVals = zeros(1,numPert);
    
    % This reshape makes two columns each with the first m items.
    lb = reshape(lb,numPert,2);
    ub = reshape(ub,numPert,2);
    
    %These initializations are to make broadcast variables smaller
    n = probInfo.n;
    Rstart = probInfo.TomoInfo.Rvar;
    angles = probInfo.anglesvar;
    PRoptions = probInfo.TomoInfo;
    
    parfor i = 1:numPert
        x = lsqAp(n,RParams(i),angleParams(i),angles(:,i),lb(i,:),ub(i,:),...
            PRoptions,optOptions,b(:,i),x_curr,Rstart(i));
        Rvals(i) = x(1);
        angleVals(i) = x(2);
    end
    paramVec = [Rvals,angleVals];
end

