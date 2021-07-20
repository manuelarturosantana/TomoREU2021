function BCDsetup
%This function adds all the codes for BCD to the file path

addpath(genpath(fileparts(mfilename('fullpath'))));
addpath('./SetFunctions','./BCDAlgorithms','./BCDAlgorithms/algorithms','./BCDAlgorithms/acceleration')
status = savepath;
if status == 1
    warning('IRbcd tools was added to the MATLAB search path for the current session only. Adding it permanently failed, probably due to a write permission issue. It is possible to manually add IR Tools permanently to your search path, but may require consulting your system administrator. Alternatively, you can re-run this installation function in each new MATLAB session where you want to use IR Tools.')
end
