function [x,iterCount] = iterTikhonov(A,b,q,tol,maxIter)
    [x_0, info] = IRhybrid_lsqr(A,b);
    alpha_0 = info.StopReg.RegP;
    iterCount = 1;
    error = 10 * tol;
    x_old = x_0;
    A_trans_b = A' * b;
    ATA = A' * A;
    [n,~] = size(ATA);
    while error > tol && iterCount < maxIter
        b_n = A_trans_b + alpha_0 * q ^ iterCount * x_old;
        B = ATA + speye(n);
        x_new = B \ b_n;
        error = norm(x_new - x_old);
        x_old = x_new;
        iterCount = iterCount + 1;
    end
    x = x_old;
end

