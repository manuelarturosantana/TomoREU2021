function paramVec = optParamParallel_var(RParams,angleParams,probInfo,budget,optOptions,b,x_curr,bounds)
% optParamParallel Wrapper function to perform the nonlinear optimization
% problem in parallel using imfil.
%   Inputs: 
%         RParams     : Current guess for the R parameters.
%         angleParams : Current guess for the angle parameters.
%         budget      : Same as budget for imfil.
%         probInfo    : Structure created by PRtomo_var
%         optOptions  : Options for lsqnonlin
%         b           : Right hand side vector
%         xk          : Current guess for the x vector.
%         bounds      : Same as bounds for imfil.
%
%    Output:
%         paramVec    : A vector with the optimized R parameters in the first
%                       half, and angle parameters in the second.
    
    numPerts = length(probInfo.TomoInfo.Rvar);
    %
    %Here we reshape b so that each column contains the correct values of b
    %for each sub problem. The syntax says we want length(b) / m  equally
    % sized rows and matlab will then fill in the number of columns
    %
    b = reshape(b,length(b) / numPerts,[]);
    
    %Initialize the parameter vector.
    Rvals = zeros(1,numPerts);
    angleVals = zeros(1,numPerts);
    
    %This reassignment is to prevent broadcasting in the parfor loop.
    Rbounds = bounds(1:numPerts,:);
    angleBounds = bounds(1 + numPerts:end,:);
    
    %These intializations are to make broadcast variables smaller
    n = probInfo.n;
    angles = probInfo.anglesvar;
    PRoptions = probInfo.TomoInfo;
    Rstart = probInfo.TomoInfo.Rvar;
    
    for i = 1:numPerts
        probBounds = [Rbounds(i,:); angleBounds(i,:)];
        x = lsqAp_var(n,RParams(i),angleParams(i),angles(:,i),probBounds,...
            budget,PRoptions,optOptions,b(:,i),x_curr,Rstart(i));
        %Each sub problem only depends on one R and angle val.
        Rvals(i) = x(1);
        angleVals(i) = x(2);
    end
    paramVec = [Rvals,angleVals];
end


