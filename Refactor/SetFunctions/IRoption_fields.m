function allfields = IRoption_fields
%IRoption_fields Define option field names for IR methods 
%
% This function is used by IRset and IRget to create a cell array of
% all possible fields for the input options structure to IR methods.

% Silvia Gazzola, University of Bath
% Per Christian Hansen, Technical University of Denmark
% James G. Nagy, Emory University
% April, 2018.

% This file is part of the IR Tools package and is distributed under the 
% 3-Clause BSD License. A separate license file should be provided as part 
% of the package.

allfields = {'x0'; 
    'MaxIter'; 
    'MaxIterIn'; 
    'MaxIterOut';
    'InProjIt';
    'x_true';
    'AllX';
    'NoStop'; 
    'NoStopIn'; 
    'NoStopOut'; 
    'NE_Rtol'; 
    'IterBar'; 
    'BestEnrm'; 
    'RegMatrix'; 
    'RegParam'; 
    'Reorth'; 
    'restart'; 
    'TotIterMax'; 
    'ktotcount'; 
    'Ktot'; 
    'stopCrit'; 
    'stopGCV'; 
    'resflatTol';
    'jbdTol';
    'jbdTolx';
    'sirt_method'; 
    'GCVweight'; 
    'GCVflatTol';
    'regPflatTol';
    'GCVminTol'; 
    'tolX'; 
    'stopOut'; 
    'stabOut';
    'thr0'; 
    'LSQRtols'; 
    'omega'; 
    'stoprule'; 
    'taudelta'; 
    'transf'; 
    'none';
    'enrichment'; 
    'stdCGLS_out'; 
    'noiseType'; 
    'noiseParam'; 
    'MatUpdate';
    'tollalpha'; 
    'verbosity'; 
    'nonneg'; 
    'Ubound'; 
    'relaxParam'; 
    'RR_Rtol';
    'NoiseLevel';
    'eta'; 
    'RegParam0';
    'RegParamVect0';
    'xMin';
    'xMax';
    'xEnergy';
    'trunc';
    'inSolver'; 
    'adaptConstr';
    'nonnegativity';
    'DecompOut';
    'shrink';
    'stepsize';
    'backtracking';
    'backit';
    'backscalar';
    'SparsityTrans';
    'hybridvariant';
    'wname';
    'wlevels';
    'warmrestart';
    'weight0';
    'weight1'; %%% to be possibly removed
    'svdbasis';
    'qnorm';
    'dimension';
    'reginskaExp';
    'plotty';
    'RegParamRange';
    'RegParRegRange';
    'discrflatTol';
    'discrbilStopTol';
    'regbilStopTol';
    %Newly Added
    'BCDStopTol'; %1e-3
    'RBounds'; %0.5
    'angleBounds'; %0.5
    'nonlinSolver'; %'imfil' or 'nonlinsolver'
    'accel'; %'anderson''ironstuck',secant,'none'
    'BCDmaxIter'; %This is here so we can pass this function into the solver as well
    'BCDlsSolver'; %"gcls","lsqr","fista","irn"
    'dispIter';
    %Anderson Acceleration Only
    'maxRes'; %'3'
    'dropTol'; %'1e6'
    %BCD solver parameters
    'budget'; %200 * numParameters
    'funcDelt'; %1e-6
    };
