function [x,IterInfo] = IRbcd(b,options,probInfo)
    
    
    % The first part of the BCD this function should do is come up with th
    % initial guess, first by generating A
    %[A, ~, ~, ~] = PRtomo_var(n, Rguess, angles_guess(:), ProbOptions);
    
    % Then The function needs to used the passed in b to make the intial
    % guess for x
    %[x0, info2] = options.Functiontouse(A, b);
    
    % Since this is the first round of BCD in generating the guess
    % options.maxIter = maxIter - 1
    
    %Next the initialization of R Params and Theta Params as zero guess
    %RParams = ones(1,m) * Rnoise_guess;
    %angleParams = ones(1,m) * ang_noise_guess; 
    %Probably just initialize to 0 for the rest of of the problem
    
    % Then create the options and bounds from the passed in information
    %
    %if isImfil
    %     imOptions = imfil_optset('least_squares',1,'simple_function',1, ...
    %     'function_delta', func_delt);
    %     bounds = [ones(1,m) * R_lower ones(1,m)* angle_lower; ...
    %     ones(1,m) * R_upper ones(1,m) * angle_upper]';
    % else %Set up the parameters for the lsqnonlin
    %     optOptions = optimoptions('lsqnonlin','MaxFunctionEvaluations',budget,...
    %         'FunctionTolerance',func_delt, 'UseParallel',false);
    %     lb = [ones(1,m) * R_lower ones(1,m) * angle_lower];
    %     ub = [ones(1,m) * R_upper ones(1,m) * angle_upper];
    %     
    % end
    
    % Here we collect the error norms for the info section.
%     xErrors = [norm(x0 - xtrue)/norm(xtrue)];
%     pErrors = [norm(paramTrue - [RParams angleParams])/norm(paramTrue)];
%     RErrors = [norm(RPert - RParams) / norm(RPert)];
%     angErrors = [norm(angleParams - angle_pert)/norm(angle_pert)];

    % Then a big if else block for the bcd
%     if strcmp(options.accel,'anderson')
%     elseif strcmp(options.accel,'ironstuck')
%     elseif strcmp(options.accel,'secant')
%     elseif strcmp(options.accel,'none')
%     else
%         error
%     end

%Finally, update the iter info and return it.
end