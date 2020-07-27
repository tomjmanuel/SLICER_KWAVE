% Call nii2MRARFI
% Tom
% 5/13/20

% This script calls the function nii2thermometry
% using the uigetfile, select the thermometry magnitude .nii file

% prior to using this script, either download the files Gstudy and use the 
%   "on the fly" converter, or install and use dcm2niix converter on your
%   own machine to conver the parrec file. The output of either of these 
%   options will be 4 .nii images; mag,real,phase,imaginary
%   The mag and phase image need to be in the same directory when you call
%   this
% select the thermometry magnitude file
[file, path] = uigetfile('*.nii','Select thermometry magnitude .nii file');

% create full path file name to magnitude file
fni = strcat(path,file);

% create full path output file name (in this case, append _arfi)
fno = strcat(path,file(1:end-4),'_Thermom.nii');

%relative value (0 to 1) to mask displacement image with using magnitude
maskVal = 0.05; 

%% call function (it writes out a nifti, nothing is returned)
dynamicSlice = 12; %The dynamic to use for image formation, to give temp contrast

% a figure will pop up visualizing the slice selected, rerun the function
% with a new dynamic selected if you don't like it.
nii2thermometry(fni,fno,maskVal,dynamicSlice);