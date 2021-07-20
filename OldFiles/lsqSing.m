function p = lsqSing(n,Rstart,angles,lb,ub,PRoptions,b,xk,isR)
%Nesting function for performing the minimization. of A(p)x_k - b
% using the fminbd with only one parameter.
%
% Input:
%      n         : Problem size. The size of the image should be n x n.
%      param     : The parameter guess for either R or theta
%      angles    : Assumed angles for tomography problem. Should be in the
%                  same format as PRtomo_var.
%      lb        : Lowerbound on parameter
%      ub        : Upperbound on parameter.
%      PRoptions : Structure used in PRtomo
%      b         : The RHS vector in A(p)x = b
%      xk        : The current approximation of the x vector.
%      Rstart    : The value which to add the R parameter to.
%      isR       : If true minimizes the R noise. If false minimizes the
%                  theta noise.
%
%      Output:
%      p: The minimized parameter.
    p = fminbnd(@multAndSub,lb,ub);
    %This is a nested function so the pass in parameters are avaliable to 
    % it, but it only takes on parameter for the fminbnd
    function paramVal = multAndSub(x)
        %This if/else block allows the function to choose if it is
        % optimizing the R noise or the angle noise.
        if isR
            Rval = Rstart  + x;
            angles2 = angles;
        else
            Rval = Rstart;
            angles2 = angles + x;
        end
        Ap = PRtomo_var(n,Rval,angles2,PRoptions);
        %To minimize this as a function of one variable we take the norm.
        paramVal = norm(Ap * xk - b);  
    end
end