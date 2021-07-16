% This is the beginning of the documentation for the new PRset and IRset
% function.
%
% PRset new options 
%
% Rvar     -  vector containing the r value for each time the angles are
% divided in the matrix. The length of Rvar must divide the the length of
% the angles.
%
% NOTE: anglesvar was not added since the length of Rvar determines how
% many columns the angles are divided into. Therefore the forming of
% anglesvar is done under the hood. If the user wishes to change the angles
% they may do so in the normal way for PRset.
%
%            [vector of positive scalars | {[2,2,2,2]}]
% Rpert     -  Scalar multiple on the perturbations on R added as Rpert * (rand() - 0.5)  
%               [scalar | {0.5}]
% anglespert     -  Scalar multiple on the perturbations on the angles added as anglespert * (rand() - 0.5)  
%               [scalar | {0.5}]
% bnoise     - Noise level to add to such that
%              bn = b + noise and ||noise(:)||_2 / ||b(:)||_2 = NoiseLevel
%              [scalar | {0.1}]
%
% Image may be the only other one that they would want to change for this.

% IRset new options
%
% BCDstopTol  - Stopping tolerance for the relative norm change between
%               iterations in the BCD loop. Different than the one provided
%               in IRset to allow the user to still set a stopping tol for
%               the least squares solver function.
%               [positive scalar | {1e-3}]
%
% Rbounds    - Percentage of R to use as a bound above and below.
%              Lowerbound is max(Rvar) - max(Rvar(i)) * Rbounds, and
%              upperbound is max(Rvar(i)) + max(Rvar(i)) * Rbounds.
%              [positive scalar | {0.5}]
%
% angleBounds  - Percentage of difference between angles to be used for the
%              bounds. Assumes angles are equally spaced.
%              [positve scalar | {0.5}]
%
% nonlinSolver - Which solver to use for the non-linear least squares
%               problem. imfil requires the imfil package and lsqnonlin
%               requires the matlab optmization toolbox.
%               ['imfil' | {'lsqnonlin'}]
%
% accel        - Determines which fixed point acceleration technique is
%                used.
%                ['anderson','ironstuck','secant',{'none'}]
%
% BCDmaxIter   - Maximum number of iterations of the BCDloop. Seperate from
%                the standard max iterations in IRset so the user can use
%                IRset for the least squares solver too.
%                [positive integer | {10}]              
%
% BCDlsSolver  - Which solver from IR tools to use to solve the
%                regularization problem. Seperate from standard version in
%                IRset to allow user to use IR set for least squares solver
%                too.
%                [{'lsqr'},'cgls','irn','fista']
% dispIter     - Display the iteration during each loop of the BCD
%
% Anderson Acceleration Only
%      maxRes - The maximum number of residuals allowed in the anderson
%               acceleration.
%               [positive integer | {3}]
%       
%      dropTol - Tolerance to drop a column in the QR factorization is the
%                condition number of the residual matrix is above dropTol.
%                [positive scalar | {1e6}]
%
% Optimization Parameters
%     budget - Maximum number of function evaluations that can be used in
%              evaluating finding the minimum. 
%              [positive integer | {2 *100 * number of perturbations}]
%
%    funcDelt - Value for stopping criterion of size of change in function 
%               while doing the optimization
%              [positive scalar | {1e-6}]