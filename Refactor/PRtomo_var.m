function [b,ProbInfo] = PRtomo_var(n, ProbOptions)
% PRtomo_var Generates data for X-ray tomographic reconstruction problems
% with  variable geometry information
%
% [b,ProbInfo] = PRtomo_var(n,ProbOptions)
% 
% Input:
%   n          :The size of the image such that the size is n x n
%   ProbOptions: Structure containing the following optional fields. Note
%               some fields avaliable in PRtomo (e.g span, CTtype) are 
%               automatically set by this problem.
%          phantomImage - user supplied test image of size n x n, of type numeric,
%               2-D only, or character string indicating
%               'shepplogan' : Shepp-Logan phantom
%               'smooth' : a smooth image
%               'binary' : a binary image
%               'threephases' : random image with pixel values 0, 0.5, 1
%                  arranged in domains
%               'threephasessmooth' : similar to threephases, but the
%                  domains have smoothly varying pixel values and there is
%                  a smooth background
%               'fourphases' : similar to 'binary' but with three phases
%                  separated by (thin) structures that form the fourth phase
%               'grains' : a random image with Voronoi cells
%               'ppower' : a random image with patterns of nonzero pixels
%               Default: 'shepplogan'.
%               This image is then stored in the output vector ProbInfo.true.xtrue.
%
%   angles:    vector of projection angles, in degrees 
%               [vector of positive scalars |{0:2:358}]
%
%   Rvar:       vector containing the r value for each time the angles are
%                divided in the matrix. The length of Rvar must divide the 
%                 length of the angles.
%                [vector of positive scalars | {[2,2,2,2]}]
%
%   Rpert     -  Scalar multiple on the perturbations on R added as Rpert * (rand() - 0.5)  
%               [scalar | {0.5}]
%
%   anglespert     -  Scalar multiple on the perturbations on the angles 
%                added as anglespert * (rand() - 0.5)  
%               [scalar | {0.5}]
%
%    bnoise     - Noise level to add to such that bn = b + noise and 
%                 ||noise(:)||_2 / ||b(:)||_2 = NoiseLevel  
%                 [scalar | {0.1}]
%          
% Output:
%  b : Vector with projection data (the sinogram), with noise added.
%  ProbInfo: Structure with the following fields
%          n: Problem size such that the problem is n x n.
%  anglesvar: Passed in angles resized into length(Rvar) columns
%       true: Structure containing the true problem data in the followin
%             fields
%            anglePert: True perturbations on the angles.
%            angles   : True angles
%            Rpert    : True perturbations on the R vals.
%            Rvar     : True R values.
%            x        : True images x
%       ProbParams: Structure for PRshowx containing the following fields.
%                     problemType - kind of test problem (in this case: 'tomography')
%                     xType       - solution type (in this case 'image2D')
%                     bType       - data type (in this case 'image2D')
%                     xSize       - size of image x
%                     bSize       - size of sinogram b
%       TomoInfo  : Structure passed in as Proboptions, to be used by the
%                   IRbcd function.

%Set the default options
ProbOptions = PRtomo_varDefaults(ProbOptions);

ProbInfo.n = n;

%Make sure that the number of R values divides number of  
 numAngles = length(ProbOptions.angles);
 numPerts  = length(ProbOptions.Rvar);
 if numAngles / numPerts ~= fix(numAngles / numPerts)
     error('Number of Rvals must divide the number of angles.')
 end

% Reshape the angles to have one for each R column
anglesvar = reshape(ProbOptions.angles,[],numPerts);
ProbInfo.anglesvar = anglesvar;

% Add the perturbations and save the true information
anglePertTrue = ProbOptions.anglespert * (rand(1,numPerts) - 0.5);
anglesTrue    = anglesvar + anglePertTrue;
Rperttrue     = ProbOptions.Rpert * (rand(1,numPerts) - 0.5);
Rtrue         = ProbOptions.Rvar + Rperttrue;

%Save the true information
ProbInfo.true.anglePert = anglePertTrue; 
ProbInfo.true.angles    = anglesTrue;
ProbInfo.true.Rpert     = Rperttrue;
ProbInfo.true.Rvar      = Rtrue;


createOpts = PRset_var(ProbOptions, 'R', Rtrue(1),'angles',anglesTrue(:,1));
[~,b, xTrue, ProbParams] = PRtomo(n, createOpts);
%save the true image.
ProbInfo.true.x = xTrue;
 for j = 2:numPerts
     createOpts = PRset_var(ProbOptions, 'R', Rtrue(j),'angles',anglesTrue(:,j));
     [~, bj, ~, ~] = PRtomo(n, createOpts);
     b = [b; bj];
 end
ProbParams.bSize = size(b);
b = PRnoise(b,ProbOptions.bnoise);

%Save the information needed for PRshowx later
ProbInfo.ProbParams = ProbParams;
ProbInfo.TomoInfo = ProbOptions;

end


%Subfunction  ----------------------------------
function options = PRtomo_varDefaults(ProbOptions)
%This function takes in the ProbOptions and sets the defaults and several
%other parameters needed.
    if isempty(ProbOptions.Rvar)
        ProbOptions.Rvar = [2,2,2,2];
    end
    if isempty(ProbOptions.Rpert)
       ProbOptions.Rpert = 0.5; 
    end
    if isempty(ProbOptions.anglespert)
        ProbOptions.anglespert = 0.5;
    end
    if isempty(ProbOptions.angles)
        ProbOptions.angles = 0:2:358;
    end
    if isempty(ProbOptions.bnoise)
        ProbOptions.bnoise = 0.1;
    end
    ProbOptions.CTtype = 'fancurved';
    %intialize the span
    span = 2*atand(1/(2*max(ProbOptions.Rvar)-1));
    ProbOptions.span = span;
    options = ProbOptions;
end