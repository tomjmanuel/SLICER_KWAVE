% nii2thermometry
% Tom
% 8/7/2020

% This script depends on PhaseUnwrapping2D if you need phase unwrapping
% (often you do not)
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
dynamicSlice = 13; %The dynamic to use for image formation, to give temp contrast
B0 = 7.0; %Magnet field, 3 or 7 T?
dyn = dynamicSlice;

% use real and imaginary images to compute phase
% for whatever reason, just using the phase is often noisier
real_fn = strcat(fni(1:end-4),'_real.nii');
im_fn = strcat(fni(1:end-4),'_imaginary.nii');
Rim = double(niftiread(real_fn));
Iim = double(niftiread(im_fn));
ImInfo = niftiinfo(im_fn);

Rim = Rim.*ImInfo.MultiplicativeScaling+ImInfo.AdditiveOffset;
Iim = Iim.*ImInfo.MultiplicativeScaling+ImInfo.AdditiveOffset;

%%
% construct imaginary image
ImComp = complex(Rim,Iim);

% image dimensions are:
% [X Y slice time]
% I'm developing using images with 5 slices
nSlices = size(ImComp,3);
nDyn    = size(ImComp,4);

% grab baseline from 3rd to 5th images
bL = squeeze(mean(ImComp(:,:,:,3:6),4));

% get echo time in milliseconds (only works for 2digit echotimes)
TEstring = ImInfo.raw.descrip(4:5);
TEms = str2double(TEstring);

if isnan( TEms)
    'warning: TE unknown, not temperature'
    TEms = 10;
end

% factor to convert phase to temp 
alpha = 0.01;
factor = 1.0 / (42.576*alpha*B0*(TEms*1e-3)*2*pi);
tempIm = angle(bL .* conj(squeeze(ImComp(:,:,:,dynamicSlice))));

%% tune mask using mask Val
maskVal = .08; % set between 0 and 1 (0.1 is working in phantoms)
mask = squeeze(mean(abs(ImComp(:,:,:,3:6)),4));
mask = mask./max(mask(:));
mask= mask>maskVal;

close all
imagesc(mask(:,:,3))

%% convert from phase to temp and find good contrast
tempIm = factor.*tempIm.*mask;
%%
% use the imtool to pick out the best data range
% click the black and white icon in the window that pops up
% this lets you scale the histogram
% tune the histogram and then copy the values from the "window"
% section into the mm and MM variables below
frame = 3; %frame with focus in it to scale contrast
imtool(tempIm(:,:,frame),'InitialMagnification','fit');

%% grab and visualize new contrast scaled temperature
% using max and min in the imtool window section, write those down here
% at the end of this section, foo should have values from 0 to 1
close all
mm = 0; % minimum deg C
MM = 2.2; % maximum deg C
% scale the image such that mm goes to zero and MM goes to 1
foo = tempIm;
foo = foo-mm; %mm values become zero with this shift, MM become MM-mm
foo(foo<0)=0; %all things below zero can be zero now
foo(foo>(-mm+MM))=(-mm+MM); 
foo = foo./(-mm+MM);
imagesc(foo(:,:,frame))

%% now put temperature values into an image with appropriate orientation
tempIm = foo;
info = ImInfo;
info.ImageSize = [info.ImageSize(1) info.ImageSize(2) nSlices];
info.PixelDimensions = [info.PixelDimensions(1) info.PixelDimensions(2) info.PixelDimensions(3)];
info.Datatype = 'double';
info.AdditiveOffset = 0;
info.MultiplicativeScaling = 1;

niftiwrite(tempIm,fno,info);
close all