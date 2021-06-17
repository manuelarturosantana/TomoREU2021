big = 110; small = 1.0256; %calculated singular values
cgOptions = IRset('RegParam',15);
n = 64;
AOptions = PRset('R', 2, 'CTtype', 'fancurved','angles',0:2:358, 'phantomImage','sheppLogan');
[A,b,x,ProbInfo] = RAnoise(n,AOptions,4,1,1);
% [A2,~,~,~] = PRtomo(n,AOptions);
Rvals = [2,2,2,2];
ThetaVals = [.1,.1,.1,.1];
A2 = makeAp(n,Rvals,ThetaVals,AOptions);
x_0 = IRcgls(A2,b,cgOptions);
disp(norm(x_0 - x)) %switch to relative norm
PRshowx(x_0,ProbInfo)
for i =1:10
    p_0 = lsAp(n,Rvals,ThetaVals,AOptions,b,x_0);
    Rvals = p_0(1:length(p_0) / 2);
    ThetaVals = p_0((length(p_0) / 2) + 1:end);
    A3 = makeAp(n,Rvals,ThetaVals,AOptions); %Is makeAP correct, and is it the correct thing to do, that is the question.
    x_0 = IRcgls(A3,b,cgOptions);
    disp(norm(x_0 - x));
end
figure(2)
PRshowx(x_0,ProbInfo);

%question 2, what about b? What about when b has noise in it?

% p_0 = lsAp(n,Rvals,ThetaVals,AOptions,b,x_0);
% Rvals = p_0(1:length(p_0) / 2);
% ThetaVals = p_0((length(p_0) / 2) + 1:end);
% A3 = makeAp(n,Rvals,ThetaVals,AOptions);
% x1 = IRcgls(A3,b,cgOptions);
% figure(4)
% PRshowx(x1,ProbInfo);
