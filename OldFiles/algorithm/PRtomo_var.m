function [A,b,x,ProbInfo] = PRtomo_var(n, R, angles, ProbOptions)
% Writen by Dr. James Nagy
% Modified by Manuel Sanata
%
% PRtomo_var Generates data for X-ray tomographic reconstruction problems
% with (possibly) variable geometry information
%
% [A,b,x,ProbInfo] = PRtomo_var(n,R,angles)
% 
% Input:
%   n:      The size of the image such that the size is n x n
%   R:      the distance from the source to the center (may not be constant)
%           of the phantom is R*N.
%           - can be a scalar, in which case R is the same for every
%             projection angle. In this case, angles will be a vector.
%           - can be a vector, in which case R will change according to
%             the structure of angles.
%   angles: projection angles, in degrees.
%           - if R is a scalar, then this is a vector with all of the
%             projection angles
%           - if R is a vector with m entries, then this is an array
%             with m columns, where it is assumed that R(j) remains
%             constant for all angles(:,j)
%  ProbOptions: structure used in PRtomo
%          
% Visualizing input R and angles: 
%           R = [   R(1)        R(2)       R(3)   ...    R(m)   ]
%      angles = [theta(1,1) theta(1,2) theta(1,3) ... theta(1,m)]
%               [theta(2,1) theta(2,2) theta(2,3) ... theta(2,m)]
%               [theta(3,1) theta(3,2) theta(3,3) ... theta(3,m)]
%                    .          .         .               .
%                    .          .         .               .
%                    .          .         .               .
%               [theta(p,1) theta(p,2) theta(p,3) ... theta(p,m)]
%
%   R: can be a scalar, so that R is always the same for every projection
%      angle
%
% Output:
%  A : Sparse matrix or function handle for forward/adjoint problem.
%  b : Vector with projection data (the sinogram).
%  x : Vector with image (i.e., exact image with stacked columns).
%  ProbInfo : structure containing some information about problem
%      problemType - kind of test problem (in this case: 'tomography')
%      xType       - solution type (in this case 'image2D')
%      bType       - data type (in this case 'image2D')
%      xSize       - size of image x
%      bSize       - size of sinogram b
%
%


if isempty(ProbOptions.phantomImage)
    image = 'sheppLogan';
else
    image = ProbOptions.phantomImage;
end

[~, m] = size(angles);
if m ~= length(R)
    error('Number columns in ''angles'' must be the same as length of ''R''')
end

ProbOptions = PRset(ProbOptions, 'R', R(1),'angles',angles(:,1),...
    'phantomImage', image);
[A, b, x, ProbInfo] = PRtomo(n, ProbOptions);
for j = 2:m
    ProbOptions = PRset(ProbOptions, 'R', R(j),'angles',angles(:,j),...
        'phantomImage', image);
    [Aj, bj, ~, ~] = PRtomo(n, ProbOptions);
    A = [A; Aj];
    b = [b; bj];
end
ProbInfo.bSize = size(b);
