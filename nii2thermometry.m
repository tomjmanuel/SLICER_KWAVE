% nii2thermometry
% Tom
% 8/7/2020

% This script depends on PhaseUnwrapping2D
% this can be found on Tom's github within the slicerkwave repo
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


%% create rough versions of temperature images
dynamicSlice = 12; %The dynamic to use for image formation, to give temp contrast

% modify fni to phase image
dyn = dynamicSlice;
phase_fn = strcat(fni(1:end-4),'_ph.nii');

%load phase image
Pim = niftiread(phase_fn);
Pinfo = niftiinfo(phase_fn);

% get Phase back into its original values (-pi to pi)
Pim = double(Pim);
Pim = Pim.*Pinfo.MultiplicativeScaling;
Pim = Pim+Pinfo.AdditiveOffset;
Pim = Pim.*.001; % they scaled up by 1000 to keep precision with int16


% load magnitude image
Im = niftiread(fni);
ImInfo = niftiinfo(fni);
Im = double(Im);

% phase image dimensions are:
% [X Y slice time]
% I'm developing using images with 5 slices
nSlices = size(Pim,3);
nDyn    = size(Pim,4);

% % We want to save a volume which represents maximum dynamic contrast
% % the input value selects a slice to use
% % display this slice for visualization
% imagesc(squeeze(Pim(:,:,1,dyn)))
% title('Slice being compared to baseline for temperature')


% grab baseline from 3rd to 5th images
bL = squeeze(mean(Pim(:,:,:,3:5),4));

% get echo time in milliseconds (only works for 2digit echotimes)
TEstring = Pinfo.raw.descrip(4:5);
TEms = str2double(TEstring);

if isnan( TEms)
    'warning: TE unknown, not temperature'
    TEms = 10;
end

% factor to convert phase to temp 
alpha = 0.01;
B0 = 7.0;
factor = 1.0 / (42.576*alpha*B0*(TEms*1e-3)*2*pi);

% get dynamic frame
dynFrame = squeeze(Pim(:,:,:,dyn));
tempIm = dynFrame-bL; % temperature image (still in phase units)

% click through these to see if any slices need unwrapping
for i=1:nSlices
    figure
    imagesc(tempIm(:,:,i),[-1.5 1.5])
    title(i)
    axis image
end

%% tune mask using mask Val
maskVal = 0.1; % set between 0 and 1 (0.1 is working in phantoms)
mask = squeeze(mean(Im(:,:,:,1:3),4));
mask = mask./max(mask(:));
mask= mask>maskVal;

close all
imagesc(mask(:,:,1))

%% unwrap phase?
% skip this section if no slices have wrapping at focus
% if a slice needs unwrapping
% if a multiple slices need it, just run this section multiple times
frame2unwrap = 3;

IM_phase=dynFrame(:,:,frame2unwrap);                         %Phase image
IM_mask = mask(:,:,frame2unwrap);
%  Set parameters
max_box_radius=30;                           %Maximum search box radius (pixels)
threshold_std=5;                            %Number of noise standard deviations used for thresholding the magnitude image

% Unwrap
residue_charge=PhaseResidues(IM_phase, IM_mask);                            %Calculate phase residues
branch_cuts=BranchCuts(residue_charge, max_box_radius, IM_mask);            %Place branch cuts
[IM_unwrapped_dyn, rowref, colref]=FloodFill(IM_phase, branch_cuts, IM_mask);   %Flood fill phase unwrapping

% Unwrap
IM_phase = bL(:,:,frame2unwrap);
residue_charge=PhaseResidues(IM_phase, IM_mask);                            %Calculate phase residues
branch_cuts=BranchCuts(residue_charge, max_box_radius, IM_mask);            %Place branch cuts
[IM_unwrapped_bL, rowref, colref]=FloodFill(IM_phase, branch_cuts, IM_mask);   %Flood fill phase unwrapping

figure
subplot(121)
imagesc(tempIm(:,:,frame2unwrap),[-1.5 1.5])
title('not unwrapped')

tempIm(:,:,frame2unwrap)=IM_unwrapped_dyn-IM_unwrapped_bL;
subplot(122)
imagesc(tempIm(:,:,frame2unwrap),[-1.5 1.5])
title('unwrapped')


%% convert from phase to temp
tempIm = factor.*tempIm.*mask;

%%now put temperature values into an image with appropriate orientation
info = Pinfo;
info.ImageSize = [info.ImageSize(1) info.ImageSize(2) nSlices];
info.PixelDimensions = [info.PixelDimensions(1) info.PixelDimensions(2) info.PixelDimensions(3)];
info.Datatype = 'double';
info.AdditiveOffset = 0;
info.MultiplicativeScaling = 1;

niftiwrite(tempIm,fno,info);
close all