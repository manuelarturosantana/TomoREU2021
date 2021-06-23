% Written by: Dr. James Nagy
% Modified by: Mai Phuong Pham Huynh

% This following code focuses on analyzing the difference between adding
% constant and random perturbations into the angles values.

% Let's set up a test problem that has "true" R values slightly 
% perturbed from our constant guess for R. Some comments about the 
% first few lines:
%  m               = number of times R changes 
%                   (in PRtomo_var this is also called m, and is the length 
%                    of vector R, or number of columns in the array angles)
%  Rguess          = guess that R is always equal to this value for all angles
%                    (since we only want to see the difference when we
%                    change the angles only, Rguess = Rtrue, i.e R remains
%                    constant for every angles)
%  Rnoise          = scalar constant on the amount of perturbation added to Rguess
%  ang_noise       = scalar constant on the amount of perturbation added to
%                    ang_true
%  ang_guess       = The guess of what the projection angles are.
%  angles_true     = true projection angles, a slight perturpation of
%                    ang_guess.
%  p               = [total number of projection angles]/m, should be an 
%                    integer (see also PRtomo_var, where p is the number of
%                    rows in the array angles)
%  span            = scalar that determines the angular span of the rays, in 
%                    degrees. I think this should be constant, even if the
%                    the source moves (and hence the geometry parameters
%                    change). We will make this constant based on the 
%                    maximum value of Rtrue.

n               = 64;
m               = 4;   
Rnoise          = 0.25;
Rguess          = 2;   
Rtrue           = Rguess*ones(1,m); 
angles_guess    = (0:2:358);
ang_noise       = 0.5;
p               = length(angles_guess)/m; 
span            = 2*atand(1/(2*max(Rtrue)-1));
ProbOptions     = PRset('CTtype', 'fancurved', 'span', span);
                          
% Check to make sure p = is an integer  
if p ~= fix(p)
    error('p = Nangles/m needs to be an integer')
end
%
% if p is an integer, reshape true and guess vector angles into an array 
% with p rows and m columns, each column corresponds to an entry in vector 
% R. Additionally add the noise to each column of the true angles.
% Note that here, instead of adding random perturbations into the angles, 
% we will add constant noise of 0.5
angles_true = angles_guess;
% Use the 2 lines below if adding random perturbations
%angle_pert = ang_noise * (rand(1,m) - 0.5);
%angles_true  = reshape(angles_true, p, m) + angle_pert;
angles_true  = reshape(angles_true, p, m) + ang_noise;
angles_guess = reshape(angles_guess,p,m);
[Atrue, btrue, xtrue, ProbInfo] = PRtomo_var(n, Rtrue, angles_true, ProbOptions);
b = PRnoise(btrue);

%
% Now setup a similar problem, using Rtrue and guess angles
% (note that using PRtomo_var with a scalar input for R and corresponding
% vector for angles, should produce the same as calling PRtomo)
%
[A, ~, ~, ~] = PRtomo_var(n, Rtrue, angles_guess, ProbOptions);

%
% Now solve using IRhybrid_lsqr to see the difference in solutions
% using the true A and the guess for A
%
[x1, info1] = IRhybrid_lsqr(Atrue, b);
[x2, info2] = IRhybrid_lsqr(A, b);

figure(1), clf
PRshowx(xtrue, ProbInfo)
title('True Solution','fontsize', 18)

figure(2), clf
PRshowx(x1, ProbInfo)
title('Solution with true angles','fontsize', 18)

figure(3), clf
PRshowx(x2, ProbInfo)
title('Solution with guessed angles','fontsize', 18)

figure(4), clf
PRshowx(abs(x2-x1), ProbInfo)
title('Difference between solution with true and guessed angles', ...
    'fontsize', 18)
