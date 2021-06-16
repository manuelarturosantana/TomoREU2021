function A = makeAp(n,Rvals,thetaVals,inOptions)
%UNTITLED Summary of this function goes here
%   thetaVals is currently the amount of noise to add, not the angles
%   themselves


angles = inOptions.angles;
numAngles = length(angles);
A = [ ];

numSwitches = numAngles / length(Rvals);

for i = 0:numSwitches - 1 %changing to zero index
    iterAngles = angles((i * 90) + 1: (i + 1) * 90) + thetaVals(i); %add noise to angles
    PRoptions = PRset('CTtype','fancurved','phantonImage',inOptions.phantomImage,...
        'angles',iterAngles,'R',Rvals(i));
    [An,~,~,~] = PRtomo(n,PRoptions);
    A = [A;An];
end

