% Written by Mai Phuong Pham Huynh
% Test set-up to answer questions from Section 6.1 of the Overleaf
% guideline from Dr. Jim Nagy

cameraman = double(imresize(imread('cameraman.tiff'), [256, 256]));

% Set seed - keep the result consistent
rng(7719)

% Set up R
Rguess = 2;
Rtrue = Rguess + 0.5*(rand(1,1)-0.5);

% Using Rtrue
options1 = PRset('phantomImage', cameraman, 'CTtype', 'fancurved', 'R', Rtrue);
[Atrue, btrue, xtrue, ProbInfo] = PRtomo(options1);
b = PRnoise(btrue);

% Using Rguess
options2 = PRset('phantomImage', cameraman, 'CTtype', 'fancurved', 'R', Rguess);
[A, ~, ~, ~] = PRtomo(options2);

% Solve for x
[x1, info1] = IRhybrid_lsqr(Atrue, b);
[x2, info2] = IRhybrid_lsqr(A, b);

% Show result
figure(1), clf
PRshowx(x1,ProbInfo)
title('Solution using Rtrue', 'fontsize', 18)

figure(2), clf
PRshowx(x2,ProbInfo)
title('Solution using Rguess', 'fontsize', 18)

