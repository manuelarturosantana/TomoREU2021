function [A,b,x,ProbInfo] = RAnoise(varargin)
% RAnoise Adds noise to the R parameter and the angles of a PR tomo problem
%
% [A,b,x,ProbInfo] = RAnoise(n,inOptions)
% [A,b,x,ProbInfo] = RAnoise(n,inOptions,Ntheta)
% [A,b,x,ProbInfo] = RAnoise(n,inOptions,Ntheta, Rnoise)
% [A,b,x,ProbInfo] = RAnoise(n,inOptions,Ntheta, Rnoise, AngNoise)
%
% Input:
%   n: The size of the image such that the size is n x n
%   inOptions: Structure contains the following fields. A sufficient
%   structure can be generated with PRset.
%       R - Assumed Distance from the source to center. Defaults to 2.     
%       angles - Vector of projection angles. Defaults to 0:2:358
%       phantomImage - 2D image or string compatable with PRset. Defaults
%       to 'sheppLogan'.
%
%   Optional Aurguments
%   Ntheta: The number of times to change the noise on R at evenly spaced 
%       increments of angles. Default length(inOptions.angles).
%   Rnoise: The scalaring constant on the amount of noise added to R. 
%       The amount of noise added is Rnoise * (rand() - 0.5). Default 0.5.
%   AngNoise: The scalaring constant on the amount of noise added to the angles. 
%       The amount of noise added is Rnoise * (rand() - 0.5). Default 0.5.
%
% Output
%   A: Matrix of true data with perturbation
%   b: Right hand side vector of true data
%   x: True x vector
%   ProbInfo: struct for use in displaying images in PRshow


% varargin allows for a different amount of aurguments to be passed into
% the function.
if length(varargin) < 2
    error("ERROR: RAnoise not enough aurguments")
end

n = varargin{1};
inOptions = varargin{2};

if isempty(inOptions.R)
    Rstart = 2;
else
    Rstart = inOptions.R;
end

if isempty(inOptions.angles)
    angles = 0:2:358;
else
    angles = inOptions.angles;
end

if isempty(inOptions.phantomImage)
    inOptions.phantomImage = 'sheppLogan';
end

%switch statement sets default values depending on the number of aurguments
%passed in.
switch length(varargin)
    case 2
        Ntheta = 1;
        Rnoise = 0.5;
        AngNoise = 0.5;
    case 3
        Ntheta = varargin{3};
        Rnoise = 0.5;
        AngNoise = 0.5;
    case 4
        Ntheta = varargin{3};
        Rnoise = varargin{4};
        AngNoise = 0.5;
    case 5
        Ntheta = varargin{3};
        Rnoise = varargin{4};
        AngNoise = varargin{5};
end

%See the documentation on PRtomo to see why this is
numRows = round(sqrt(2) * n) * length(angles);

A = [ ];

b = zeros(numRows,1);

%See PRtomo documentation and code for why this number
rowsPerIter = round(sqrt(2) * n);


changeRN = round(length(angles) / Ntheta );

RN = Rstart + Rnoise * (rand() - 0.5);

numIter = 0;

for i = angles
    anglesN = i + AngNoise * (rand() - 0.5);
    options = PRset('CTtype','fancurved','phantomImage', inOptions.phantomImage,...
        'angles',anglesN, 'R', RN);
    [An,bn,xn,ProbInfo] = PRtomo(n,options);
    A = [A; An]; %The fastest way despite reallocation in each loop.
    b(numIter * rowsPerIter + 1: (numIter + 1) * rowsPerIter) = bn;
    numIter = numIter + 1;
    if mod(numIter, changeRN) == 0
        RN = Rstart + Rnoise * (rand() - 0.5);
    end
end
ProbInfo.bSize = size(b);
x = xn;

end

