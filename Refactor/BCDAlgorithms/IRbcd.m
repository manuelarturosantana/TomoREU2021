function [x,iterInfo] = IRbcd(b,iterOptions,probInfo)
%   This function perform the block coordinate descent for tomography problems.
%         [x,iterInfo] = IRbcd(b,iterOptions,probInfo)
%
%   Inputs: 
%         b      : RHS vector or sinogram.
%    interOptions: Structure created by IRset. Used for both the linear least
%                  squares solver and the BCD loop. 
%                  The following parameters are relevant to the BCD. 

%
%   BCDstopTol  - Stopping tolerance for the relative norm change between
%                 iterations in the BCD loop. Different than the one provided
%                 in IRset to allow the user to still set a stopping tol for
%                 the least squares solver function.
%                 [positive scalar | {1e-3}]
%
%   Rbounds    - Percentage of R to use as a bound above and below.
%                Lowerbound is max(Rvar) - max(Rvar(i)) * Rbounds, and
%                upperbound is max(Rvar(i)) + max(Rvar(i)) * Rbounds.
%                 [positive scalar | {0.5}]
%
%   angleBounds  - Percentage of difference between angles to be used for the
%                  bounds. Assumes angles are equally spaced.
%                 [positve scalar | {0.5}]
%
%   nonlinSolver - Which solver to use for the non-linear least squares
%                   problem. imfil requires the imfil package and lsqnonlin
%                   requires the matlab optmization toolbox.
%                   ['imfil' | {'lsqnonlin'}]
%
%   accel        - Determines which fixed point acceleration technique is
%                   used.
%                   ['anderson','ironstuck','secant',{'none'}]
%
%   BCDmaxIter   - Maximum number of iterations of the BCDloop. Seperate from
%                   the standard max iterations in IRset so the user can use
%                   IRset for the least squares solver too.
%                   [positive integer | {10}]              
%
%   BCDlsSolver  - Which solver from IR tools to use to solve the
%                   regularization problem. Seperate from standard version in
%                   IRset to allow user to use IR set for least squares solver
%                   too.
%                   [{'lsqr'},'cgls','irn','fista']
%
%   dispIter     - Display the iteration during each loop of the BCD
%                    [{'on'},{'off}]
%
%
%    Optimization Parameters
%        budget - Maximum number of function evaluations that can be used in
%              evaluating finding the minimum. 
%              [positive integer | {2 *100 * number of perturbations}]
%
%       funcDelt - Value for stopping criterion of size of change in function 
%                 while doing the optimization
%                [positive scalar | {1e-6}]
%
%    Anderson Acceleration Only
%         maxRes - The maximum number of residuals allowed in the anderson
%                  acceleration.
%                  [positive integer | {3}]
%          
%         dropTol - Tolerance to drop a column in the QR factorization is the
%                   condition number of the residual matrix is above dropTol.
%                   [positive scalar | {1e6}]
%
%    ProbInfo   - structure output by PRtomo
%
%    Output:
%          x - The x vector solution after BCD.
%    iterInfo - structure containg the following fields.
%
%             xErrors - Array of relative image error at each iteration.
%             paramErrors - Array of relative error of the perturbations at
%                           each iteration.
%             angErrors   - Array of relative error of the angles
%                           perturbations at each iteration.
%             RErrors     - Array of the relative error the the R
%                           perterbations at each iteration.
%             numIter     - The number of iterations performed before stopping.
%             runTime     - Total run time of the BCD as measured by
%                           tic/toc
%             xsols       - array with the x vector after each iteration.
%             xtrueparam  - vector containing the x solution with the true
%                           parameters, but noise added to b.
    
    %Set up default parameters
    iterOptions = setIterOptions(iterOptions,probInfo);
    
    %Runbcd loop based on acceleration technique
    if strcmpi(iterOptions.accel,'none')
       [x,iterInfo] = NAbcd(b,iterOptions,probInfo);
    elseif strcmpi(iterOptions.accel,'anderson')
        [x,iterInfo] = AAbcd(b,iterOptions,probInfo);
    elseif strcmpi(iterOptions.accel,'ironstuck')
        [x,iterInfo] = ITbcd(b,iterOptions,probInfo);
    elseif strcmpi(iterOptions.accel,'secant')
        [x,iterInfo] = CSbcd(b,iterOptions,probInfo);
    else
        error('Acceleration technique not recognized.')
    end

end


%SubFunction--------------------------------------------------------------
function options = setIterOptions(iterOptions,probInfo)
%This function sets the defaults for the iterOptions
    if isempty(iterOptions.BCDStopTol)
        iterOptions.BCDStopTol = 1e-3;
    end
    if isempty(iterOptions.RBounds)
        iterOptions.RBounds = 0.5;
    end
    if isempty(iterOptions.angleBounds)
        iterOptions.angleBounds = 0.5;
    end
    if isempty(iterOptions.nonlinSolver)
        iterOptions.nonlinSolver = 'lsqnonlin';
    end
    if isempty(iterOptions.accel)
        iterOptions.accel = 'none';
    end
    if isempty(iterOptions.BCDmaxIter)
        iterOptions.BCDmaxIter = 10;
    end
    if isempty(iterOptions.BCDlsSolver)
        iterOptions.BCDlsSolver = 'lsqr';
    end
    if isempty(iterOptions.dispIter)
        iterOptions.dispIter = 'on';
    end
    if isempty(iterOptions.maxRes)
        iterOptions.maxRes = 3;
    end
    if isempty(iterOptions.dropTol)
        iterOptions.dropTol = 10e6;
    end
    if isempty(iterOptions.budget)
        %Imfil recommended number function evaluations. 
        iterOptions.budget = 2 * 100 *length(probInfo.TomoInfo.Rvar);
    end
    if isempty(iterOptions.funcDelt)
        iterOptions.funcDelt = 1e-6;
    end
    options = iterOptions;
end
