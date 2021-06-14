function [A,b,x,ProbInfo] = RAnoise(varargin)
% RAnoise Adds noise to the R parameter and the angles of a PR tomo problem
%
% [A,b,x,ProbInfo] = RAnoise(n,Rstart,angles)
% [A,b,x,ProbInfo] = RAnoise(n,Rstart,angles,Ntheta)
% [A,b,x,ProbInfo] = RAnoise(n,Rstart, angles, Ntheta, Rnoise, AngNoise)
% 
% Input:
%   n: The size of the image such that the size is n x n
%   Rstart: The value of R to add noise to.
%   angles: The angels to add noise to in the problem. 
%   Ntheta: The number of times to change the noise on R at evenly spaced 
%       increments of angles. Default length(angles).
%   Rnoise: The scalaring constant on the amount of noise added to R. 
%       The amount of noise added is Rnoise * (rand() - 0.5). Default 0.5.
%   AngNoise: The scalaring constant on the amount of noise added to the angles. 
%       The amount of noise added is Rnoise * (rand() - 0.5). Default 0.5.

if length(varargin) < 3
    error("ERROR: RAnoise not enough aurguments")
end

n = varargin{1};
Rstart = varargin{2};
angles = varargin{3};

switch length(varargin)
    case 3
        Ntheta = length(angles);
        Rnoise = 0.5;
        AngNoise = 0.5;
    case 4
        Ntheta = varargin{4};
        Rnoise = 0.5;
        AngNoise = 0.5;
    case 6
        Ntheta = varargin{4};
        Rnoise = varargin{5};
        AngNoise = varargin{6};
end

%See the documentation on PRtomo to see why this is
numRows = round(sqrt(2) * n) * length(angles);

A = zeros(numRows, n^2);
b = zeros(numRows,1);
x = zeros(n^2,1);

rowsPerIter = round(sqrt(2) * n);
numIter = 0;

changeRN = round(length(angles) / Ntheta );

RN = Rstart + Rnoise * (rand() - 0.5);
for i = angles
    anglesN = i + AngNoise * (rand() - 0.5);
    options = PRset('CTtype','fancurved','angles',anglesN, 'R', RN);
    [An,bn,xn,ProbInfo] = PRtomo(n,options);
    A(numIter * rowsPerIter + 1: (numIter + 1) * rowsPerIter,:) = An;
    b(numIter * rowsPerIter + 1: (numIter + 1) * rowsPerIter) = bn;
    x = x + xn;
    numIter = numIter + 1;
    if mod(numIter, changeRN) == 0
        RN = Rstart + Rnoise * (rand() - 0.5);
    end
end
ProbInfo.bSize = size(b);

end

