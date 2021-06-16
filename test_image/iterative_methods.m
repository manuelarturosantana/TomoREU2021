

% set up
spine = double(imresize(imread('spine.tif'), [256, 256]));
theta = 0:1:179;
ProbOptions = PRset('phantomImage', spine, ...
    'CTtype', 'fancurved', ...
    'angles', theta);
[A, b_true, x_true, ProbInfo] = PRtomo(ProbOptions);
[b, NoiseInfo] = PRnoise(b_true);
figure(1), clf
PRshowb(b, ProbInfo)
figure(2), clf
PRshowx(x_true,ProbInfo)
x = IRcgls(A,b);
figure(3), clf
PRshowx(x, ProbInfo)

CGLSotpions = IRcgls('defaults');
CGLSoptions = IRset(CGLSoptions, 'x_true', x_true);
[x, IterInfo] = IRcgls(A,b,CGLSoptions);

IterInfo
IterInfo.BestReg

FS = 18;
MS = 10;
LW = 2;
% plot relative error in each iteration (minimum around 11th iteration)
figure(4), clf
axes('FontSize', FS), hold on
plot(IterInfo.Enrm, 'b-', 'LineWidth', LW)
xlabel('Iteration')
ylabel('Relative error')
plot(IterInfo.BestReg.It, IterInfo.BestReg.Enrm, 'ro',...
    'MarkerSize', MS, 'LineWidth', LW), hold off
% this plot is unsatisfied because we don't know the true solution
% plot the best solution
figure(5), clf
PRshowx(IterInfo.BestReg.X, ProbInfo)

% if we know the norm of the noise level
CGLSoptions2 = IRset(CGLSoptions, 'NoiseLevel', 0.01);
[x2, IterInfo2] = IRcgls(A,b,CGLSoptions2);
IterInfo2.StopReg
figure(6), clf
axes('FontSize', FS), hold on
plot(IterInfo2.Enrm, 'b-', 'LineWidth', LW)
xlabel('Iteration')
ylabel('Relative error')
plot(IterInfo2.BestReg.It, IterInfo2.BestReg.Enrm, 'ro',...
    'MarkerSize', MS, 'LineWidth', LW), hold off
% the more noise, the less iterations because The more noise, the less 
% “true” information you can squeeze out of the data. -> less iterations
% plot the best solution
figure(7),clf
PRshowx(IterInfo2.BestReg.X, ProbInfo)

CGLSoptions3 = IRset(CGLSoptions, 'RegParam', 5.7863);
[x3, IterInfo3] = IRcgls(A,b, CGLSoptions3);



