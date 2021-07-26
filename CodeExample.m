%This script is to be sneaky and show a code example for the poster
%presentation
pause

disp(">> n = 256")
pause

disp(">> ProbOptions = PRset('Rvar',2 * ones(1,180),'Rpert',0.5,'anglespert',0.5);")
pause

disp(">> [b,probInfo] = PRtomo_var(n,ProbOptions);")
pause

fprintf(">> iterOptions = IRset('nonlinSolver','imfil','BCDmaxIter',20,...\n'Rbounds',0.1250,'angleBounds',0.1250,'BCDlsSolver','irn');\n")
pause

load IRNcontrol.mat

disp(">> [x,iterInfo] = IRbcd(b,iterOptions,probInfo);")
pause

disp(">> PRshowbcd(iterInfo,probInfo)")

PRshowbcd(iterInfo,probInfo)