function A = createA(n,Rvals,angVals,ProbOptions)
%
% createA creates the matrix A to be used in the optimization functions.
%  A = createA(Rvals,angVals,ProbOptions)
%  Note, one reason for this many pass ins is to avoid large broadcast 
%  variables in the parallelization.
%
% Input:
%   n:       The problem size such that the image is n x n;
%   Rvals:   the distance from the source to the center (may not be constant)
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
%

one = 1;

%Set up the info passed into each PRset
PRoptions = ProbOptions;

%The number of time perturbations are added
numPerts = length(Rvals);

PRoptions = PRset(PRoptions, 'R', Rvals(1),'angles',angVals(:,1));
[A,~,~,~] = PRtomo(n, PRoptions);
for j = 2:numPerts
    PRoptions = PRset(PRoptions, 'R', Rvals(j),'angles',angVals(:,j));
    [Aj, ~, ~, ~] = PRtomo(n, PRoptions);
    A = [A; Aj];
end

end
