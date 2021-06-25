% To use this script enter a .mat file created by EulerRunTomo, and it will
% plot the same information as RunTomoTestParallel would. The mat file
% should contain at least the following variables.
%
% xtrue    : The true solution of the tomography problem
% x1       : The solution with the true parameters on A.
% x2       : The solution with the initially guessed parameters.
% x_k      : The solution after BCD has been performed.
% isImfil  : Boolean indicating if imfil or nonlinlsqr was used.
% xErrors  : Relative errors of x at each iteration.
% pErrors  : Relative errors of p at each iteration.
% RErrors  : Relative errors of R parameters at each iteration.
% angErrors: Relative errors of the angles at each iteration.

filename = input("Please input filename as a string or variable: ");

results = load(filename);

figure(1), clf
PRshowx(results.xtrue, ProbInfo)
title('True Solution','fontsize', 20)

figure(2), clf
PRshowx(results.x1, ProbInfo)
title('Solution with True A','fontsize', 20)

figure(3), clf
PRshowx(results.x2, ProbInfo)
title('Solution with Noisey A','fontsize', 20)

figure(4), clf
PRshowx(results.x_k,ProbInfo)
if results.isImfil
    title('Solution After BCD (Imfil)','fontsize', 20)
else 
    title('Solution After BCD (lsqnonlin)','fontsize', 20)
end

figure(5), clf
plot(results.xErrors,'-o');
hold on 
plot(results.pErrors,'-*');
plot(results.RErrors,'-d');
plot(results.angErrors,'-^');
hold off
legend('xError Norms','p Error Norms', 'R error Norms','Angle Error Norm');
xlabel('Number of Iterations','fontsize',15);
ylabel('Relative error','fontsize',15);

