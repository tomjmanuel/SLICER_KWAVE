% Call nii2MRARFI
% Tom
% 5/13/20

% This script calls the function nii2MRARFI
% using the uigetfile, select the arfi magnitude .nii file

% prior to using this script, either download the files Gstudy and use the 
%   "on the fly" converter, or install and use dcm2niix converter on your
%   own machine to conver the parrec file. The output of either of these 
%   options will be 4 .nii images; mag,real,phase,imaginary
%   All four need to be in same directory

% select the arfi magnitude file
[file, path] = uigetfile('*.nii','Select Arfi magnitude .nii file');

% create full path file name to magnitude file
fni = strcat(path,file);

% create full path output file name (in this case, append _arfi)
fno = strcat(path,file(1:end-4),'_arfi.nii');

%relative value (0 to 1) to mask displacement image with using magnitude
maskVal = 0.3; 

% call function (it writes out a nifti, nothing is returned)
nii2MRARFI(fni,fno,maskVal);