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


%% generate arfi image and info to write nifti

% call function
[arfIm, info, magIm] = nii2MRARFI(fni);
close all

%% interactively scale contrast
% click the black white icon to interactively scale contrast
% when done hit adjust data 
% then in figure file->export to workspace to save the image data as a 
% variable (I save as arfi2)
imtool(arfIm,'InitialMagnification','fit')

%% select mask value
% make sure 2 change line 37 to match whatever variable you exported
newarf = arfi2;
maskVal = .4;
newarf = newarf.*(magIm>maskVal);
imagesc(newarf)

%% write newarf out to nifti
ni = zeros(info.ImageSize);
ni(:,:,1) = newarf;
ni(:,:,2) = newarf;

niftiwrite(ni,fno,info);



