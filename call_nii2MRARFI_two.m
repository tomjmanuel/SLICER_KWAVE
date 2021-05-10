% Call nii2MRARFI
% Tom + Huiwen
% 10/13/20
% we had to update this script to work with 3T and to use bipolar 
% gradients

% This script calls the function nii2MRARFI_two
% using the uigetfile, select the arfi magnitude .nii file

% prior to using this script, either download the files Gstudy and use the 
%   "on the fly" converter, or install and use dcm2niix converter on your
%   own machine to conver the parrec file. The output of either of these 
%   options will be 4 .nii images; mag,real,phase,imaginary
%   All four need to be in same directory

%For Huiwen's ARFI sequence, interleaved fashion is on.
%There are two mixes for each scan: 1) FUS on, 2) FUS off.
% naming convention_p1 is +grad, p2 is -grad

% sometimes we don't collect both directions if image is good enough
% if only collecting positive use this
pos_only = 0;

% load the scan with positive motion-encoding gradients (MEG)
[file, path] = uigetfile('*.nii','Select Arfi magnitude +grad .nii file');
% create full path file name to magnitude file
fni_p = strcat(path,file);

if (pos_only)
    fni_m = fni_p;
else
    % load the scan with negative motion-encoding gradients (MEG)
    [file, path] = uigetfile('*.nii','Select Arfi magnitude -grad .nii file');
    % create full path file name to magnitude file
    fni_m  = strcat(path,file);
end

fni = {fni_p, fni_m};
% create full path output file name (in this case, append _arfi)
fno = strcat(path,file(1:end-4),'_arfi.nii');

%%
% call function (it writes out a nifti, nothing is returned)
% 3rd slice
gradStrength = 40; %mT/m

MEGdur = 8; %ms
[arfi, info, magIm] = nii2MRARFI_two(fni,MEGdur,gradStrength,pos_only);
magIm = magIm(:,:,:,1);
%% permute arfi if smallest dimension isn't last
sliceDim = 3;   % set this variable to hold slice dimension
                % for instance, if dim(arfi) = 160x5x160, sliceDim=2
                % if dim(arfi) = 160x160x5, sliceDim = 3

perm_flag = 0; % flag to denote if image was permuted
if sliceDim ==2
    arfi = permute(arfi,[1,3,2]);
    magIm = permute(magIm,[1,3,2]);
    perm_flag=1;
elseif sliceDim == 1
    arfi = permute(arfi,[2,3,1]);
    magIm = permute(magIm,[2,3,1]);
    perm_flag=1;
end


%% subtract out background phase
% This portion of the code was written by Sumeeth
% it subtracts background phase...
% make sure to click a region in the brain or phantom but not at the focus
sliceOfInterest = 2; % set this to be the slice you want to choose region in
figure
imagesc(-arfi(:,:,sliceOfInterest),[-1 1])
title('click on region to use for phase background')
pos = ginput(1);
Nx = size(arfi,1); Ny = size(arfi,2);
[xx,yy] = ndgrid((1:Nx)-pos(2),(1:Ny)-pos(1)); % mask location can have a big effect on results
mask = (xx.^2 + yy.^2)<6^2;
nSlice = size(arfi,3);
for i = 1 : nSlice % subtract out background phase from a circular ROI ("mask") near the focus; repeat for each dynamic
    tmp1 = arfi(:,:,i);
    arfi(:,:,i) = tmp1 - mean(tmp1(mask));
end
close all
%% run this to pull up imtool and scale contrast
imtool(arfi(:,:,sliceOfInterest),'InitialMagnification','fit');
%% next rescale intensity values using imtool for better contrast
% click the black and white circle in the imtool window
% use the vertical bars in the histogram to adjust contrast
% using max and min in the imtool window section, write those down here
% at the end of this section, foo should have values from 0 to 1
close all
mm = 0.;
MM = 1.5;
% scale the image such that mm goes to zero and MM goes to 1
foo = arfi;
foo = foo-mm; %mm values become zero with this shift, MM become MM-mm
foo(foo<0)=0; %all things below zero can be zero now
foo(foo>(-mm+MM))=(-mm+MM); 
foo = foo./(-mm+MM);
imagesc(foo(:,:,sliceOfInterest))

%% finally crop image with a magnitude mask
magIm = magIm- min(magIm(:)); % scale 0 to 1
magIm = magIm./max(magIm(:));
%%
maskVal = .05; % 0 to 1
mask = magIm>maskVal;
%mask = repmat(mask,[1,1,nSlice]); % replicate mask through all slices
foo2 = foo.*mask; % apply mask
% visualize result (try new maskVal if not satisfied)
close all
imagesc(foo2(:,:,sliceOfInterest));

%% now repermute if that was applied
arfi = foo2;
if perm_flag
    if sliceDim ==2
        arfi = permute(arfi,[1,3,2]);
    end
    if sliceDim ==1
        arfi = permute(arfi,[2,3,1]);
    end
end

%% sometimes info.Imagesize has a fourth dimension
dim = size(arfi);
vox = info.PixelDimensions(1:3);
info.ImageSize = dim;
info.PixelDimensions = vox;

%% write arfi out to nifti using header from nifti's read into this code

info.Datatype = 'double';
ni = arfi;
niftiwrite(ni,fno,info);
