% Given a nifti file, compute ARFI displacement
% output this to a nifti file with the appropriate imagespace information
% so that it can be loaded into Slicer and overlaid on anatomical image

% Input:
%   fni: input filename, path + filename
%       points at a .nii file which is one of 4 .nii output by dcm2niix
%       - dcm2niix can be used directly from gstudy or
%       it can be installed and run from the command line on your pc to
%       convert par rec files into niftis
%       - Given 1 par file, dcm2niix returns 4 seperate nii
%       - These all need to be in same path 
%   fno: output filename
%   maskVal:percentage (0 to 1)
%           mask displacement image with magnitude image>maskVal
%           just set to 0 if you don't want to mask
%   dyn:    dynamic number to use for temperature contrast
%           i.e. what slice to compare to baseline


function nii2thermometry(fni,fno,maskVal,dyn)

    % modify fni to phase image
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
    
    % We want to save a volume which represents maximum dynamic contrast
    % the input value selects a slice to use
    % display this slice for visualization
    imagesc(squeeze(Pim(:,:,1,dyn)))
    title('Slice being compared to baseline for temperature')
    
    
    % grab baseline from 3rd to 5th images
    bL = squeeze(mean(Pim(:,:,:,3:5),4));
    
    % get echo time in milliseconds (only works for 2digit echotimes)
    TEstring = Pinfo.raw.descrip(4:5);
    TEms = str2double(TEstring);

    % factor to convert phase to temp 
    alpha = 0.01;
    B0 = 7.0;
    factor = 1.0 / (42.576*alpha*B0*(TEms*1e-3)*2*pi);
    
    % get dynamic frame
    dynFrame = squeeze(Pim(:,:,:,dyn));
    
    % convert from phase to temp
    tempIm = factor.* (dynFrame-bL);
    
    % mask using mask Val
    mask = squeeze(mean(Im(:,:,:,1:3),4));
    mask = mask./max(mask(:));
    mask(mask<maskVal)=0;
    tempIm = tempIm.*mask;
    
    % blur to remove noise
    %tempIm = imgaussfilt(tempIm,.5);

    % now put temperature values into an image with appropriate orientation
    info = Pinfo;
    info.ImageSize = [info.ImageSize(1) info.ImageSize(2) nSlices];
    info.PixelDimensions = [info.PixelDimensions(1) info.PixelDimensions(2) info.PixelDimensions(3)];
    info.Datatype = 'double';
 
    niftiwrite(tempIm,fno,info);

end