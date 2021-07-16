function [x,iterInfo] = IRbcd(b,iterOptions,probInfo)
    
    %TODO functionality to return x from true parameters
    %[x1, info1] = IRhybrid_lsqr(Atrue, b);
    
    %Set up default parameters
    iterOptions = setIterOptions(iterOptions,probInfo);
    
    %Runbcd loop based on acceleration technique
    if strcmpi(iterOptions.accel,'none')
       [x,iterInfo] = NAbcd(b,iterOptions,probInfo);
    elseif strcmpi(iterOptions.accel,'anderson')
        [x,iterInfo] = AAbcd(b,iterOptions,probInfo);
%     elseif strcmpi(iterOptions.accel,'ironstuck')
%         [x,iterInfo] = ITbcd(b,iterOptions,probInfo);
%     elseif strcmpi(iterOptions.accel,'secant')
%         [x,iterInfo] = CSbcd(b,iterOptions,probInfo);
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
