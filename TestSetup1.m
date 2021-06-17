%
% Let's set up a test problem that has "true" R values slightly 
% perturbed from our constant guess for R. Some comments about the 
% first few lines:
%  m           = number of times R changes 
%                (in PRtomo_var this is also called m, and is the length 
%                 of vector R, or number of columns in the array angles)
%  Rguess      = guess that R is always equal to this value for all angles
%  Rtrue       = true values of R, which may be slight perturbations of Rguess
%  Rnoise      = scalar constant on the amount of perturbation added to Rguess
%  angles_true = true projection angles
%  p           = [total number of projection angles]/m, should be an 
%                integer (see also PRtomo_var, where p is the number of
%                rows in the array angles)
%  span        = scalar that determines the angular span of the rays, in 
%                degrees. I think this should be constant, even if the
%                the source moves (and hence the geometry parameters
%                change). We will make this constant based on the 
%                maximum value of Rtrue.
%
n            = 128;
m            = 4;   
Rnoise       = 0.5;
Rguess       = 2;   
Rtrue        = Rguess*ones(1,m) + Rnoise*(rand(1,m) - 0.5);
angles_true  = (0:2:358); 
p            = length(angles_true)/m; 
span         = 2*atand(1/(2*max(Rtrue)-1));
ProbOptions  = PRset('CTtype', 'fancurved', 'span', span);
                          
% Check to make sure p = is an integer  
if p ~= fix(p)
    error('p = Nangles/m needs to be an integer')
end
%
% if p is an integer, reshape vector angles into an array with p rows
% and m columns, each column corresponds to an entry in vector R
%
angles_true  = reshape(angles_true, p, m);

[Atrue, btrue, xtrue, ProbInfo] = PRtomo_var(n, Rtrue, angles_true, ProbOptions);
b = PRnoise(btrue);

%
% Now setup a similar problem, using the guess for R and true angles
% (note that using PRtomo_var with a scalar input for R and corresponding
% vector for angles, should produce the same as calling PRtomo)
%
[A, ~, ~, ~] = PRtomo_var(n, Rguess, angles_true(:), ProbOptions);

%
% Now solve using IRhybrid_lsqr to see the difference in solutions
% using the true A and the guess for A
%
[x1, info1] = IRhybrid_lsqr(Atrue, b);
[x2, info2] = IRhybrid_lsqr(A, b);

figure(1), clf
PRshowx(xtrue, ProbInfo)

figure(2), clf
PRshowx(x1, ProbInfo)

figure(3), clf
PRshowx(x2, ProbInfo)
