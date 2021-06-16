n = 50; options = PRset('CTtype', 'fancurved');
[A,b,x,ProbInfo] = PRtomo(n,options);
NoiseLevel = 0.1; bn = PRnoise(b, NoiseLevel);
xinv = A\bn; PRshowx(xinv, ProbInfo)
Af = full(A);
[U,S,V] = svd(Af);

SingVal = diag(S);
error = zeros(length(SingVal),1);

for i = 1:length(SingVal)
    lambda = SingVal(i);
    xtik = [A; lambda*speye(size(A,2))]\[bn;zeros(size(A,2),1)];
    error(i) = norm(xtik-x)/norm(x);
end

optimized_error = min(error)
optimized_lambda = SingVal(error == optimized_error)
optimized_xtik = [A; optimized_lambda*speye(size(A,2))]\[bn;zeros(size(A,2),1)];

PRshowx(optimized_xtik, ProbInfo)

