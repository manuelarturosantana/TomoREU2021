function A = makeAp(n,Rvals,thetaVals,inOptions)
%   This function takes in parameters for the R values and theta Values and
%   returns the matrix A via the PRtomo function.
%
%   Inputs:
%   n: problem size. The desired image should be n x n
%   Rvals: The value of R for interval. 
%   thetaVals: The amount of noise to add to each angle theta, not the
%   angles themselves.
%   inOptions: Structure containg the angles, image
%
%   Output
%   A: A matrix built by the values given from Rvals and thetaVals

if isempty(inOptions.angles)
    angles = 0:2:358;
else
    angles = inOptions.angles;
end

if isempty(inOptions.phantomImage)
    image = 'sheppLogan';
else
    image = inOptions.phatomImage;
end

numAngles = length(angles);
A = [ ];

numPerIter = numAngles / length(Rvals);

for i = 0:(length(Rvals) - 1) %changing to zero index
    iterAngles = angles((i * numPerIter) + 1: (i + 1) * numPerIter) + thetaVals(i + 1); %add noise to angles
    PRoptions = PRset('CTtype','fancurved','phantomImage',image,...
        'angles',iterAngles,'R',Rvals(i + 1));
    [An,~,~,~] = PRtomo(n,PRoptions);
    A = [A;An];
end


