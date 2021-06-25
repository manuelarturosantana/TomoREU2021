% Written by: Mai Phuong Pham Huynh
% Set-up test for Section 6.1.2

% Set seed
rng(7719)

% Parameters
angles = 0:2:358;
Rguess = 2;
Rtrue = Rguess*ones(180,1) + 0.5*(rand(180,1) - 0.5);

% Using Rtrue
options1 = PRset('phantomImage', cameraman, 'CTtype', 'fancurved', 'R', ...
    Rtrue(1), 'angles', angles);
[Atrue, btrue, xtrue, ProbInfo] = PRtomo_var(options1);
for i=2:length(Rtrue)
    options1 = PRset('phantomImage', cameraman, 'CTtype', 'fancurved', ...
        'R', Rtrue(i), 'angles', angles);
    [Atrue_iter, btrue_iter, xtrue_iter, ProbInfo] = PRtomo(options1);
    Atrue = [Atrue; Atrue_iter];
    btrue = [btrue; btrue_iter];
end
b = PRnoise(btrue);

% Using Rguess
options2 = PRset('phantomImage', cameraman, 'CTtype', 'fancurved', 'R', ...
    Rguess, 'angles', angles);
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
