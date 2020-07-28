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


function [arfi, info, mag] = nii2MRARFI(fni)

    % modify fni to point at real and imaginary images
    realfn = strcat(fni(1:end-4),'_real.nii');
    imagfn = strcat(fni(1:end-4),'_imaginary.nii');
    
    %load real image
    R = niftiread(realfn);
    Rinfo = niftiinfo(realfn);
    
    % get R back into its original values
    R = squeeze(double(R));
    R = R.*Rinfo.MultiplicativeScaling+Rinfo.AdditiveOffset;
    
    % load imaginary image
    Im = niftiread(imagfn);
    ImInfo = niftiinfo(imagfn);
    
    % get R back into its original values
    Im = squeeze(double(Im));
    Im = Im.*ImInfo.MultiplicativeScaling+ImInfo.AdditiveOffset;
    
    data = complex(R,Im);

    % compute ARFI using a crytpic equation that "subtracts dynamics"
    arfi = angle( data(:,:,2:4:end) .* conj(data(:,:,1:4:end)) .* conj( ...
    data(:,:,4:4:end) .* conj(data(: ,:,3:4:end)) )); % subtraction of dynamics (used timing and averaging)
   
    % compute ARFI using a crytpic equation that "subtracts dynamics"
%    arfi = angle( data(:,:,4:4:end) .* conj(data(:,:,3:4:end)) .* conj( ...
%    data(:,:,6:4:end) .* conj(data(: ,:,5:4:end)) )); % subtraction of dynamics (used timing and averaging)
   

    Nx = size(arfi,1);
    Ny = size(arfi,2);
    Navg = size(arfi,3);

    % convert from radians to pi
    arfi = arfi ./ (-2*2*pi*42.58*1e6*40*1e-3*3*1e-3) * 1e6; % radians -> microns

    % This portion of the code was written by Sumeeth
    % it subtracts background phase...
    figure
    imagesc(-arfi(:,:,1),[-1 1])
    title('click on region to use for phase background')
    pos = ginput(1);
    [xx,yy] = ndgrid((1:Nx)-pos(2),(1:Ny)-pos(1)); % mask location can have a big effect on results
    mask = (xx.^2 + yy.^2)<6^2;
    %imagesc(mask)
   % mask2 =  (xx.^2 + yy.^2)<4^2;
    %maskConc = mask-mask2;
    %maskConc = maskConc>0;
    for i = 1 : Navg % subtract out background phase from a circular ROI ("mask") near the focus; repeat for each dynamic
        tmp1 = arfi(:,:,i);
        arfi(:,:,i) = tmp1 - mean(tmp1(mask));
    end
    
    % combine multiple arfi acquisitions if they exist
    arfi = mean(arfi,3);
    
    % mask the displacement image with the magnitude image
    % actually, masking with zeros isn't the best since
    % set masked values = to -10 um;
    mag = abs(mean(data,3)); mag = mag ./ max(mag(:));
%     mask = mag > maskVal;
%     arfi(~mask)=-10;
    
    % now put arfi values into an image with appropriate orientation
    % I'm only able to make this work as a multislice
    % so i make it 2 identical slices of half their actual thickness
    info = Rinfo;
    info.ImageSize = [info.ImageSize(1) info.ImageSize(2) 2];
    info.PixelDimensions = [info.PixelDimensions(1) info.PixelDimensions(2) info.PixelDimensions(3)/2];
    info.Datatype = 'double';
    info.AdditiveOffset = 0;
    info.MultiplicativeScaling = 1000;
    
%     ni = zeros(info.ImageSize);
%     ni(:,:,1) = arfi;
%     ni(:,:,2) = arfi;
    
    % rescale values to be easier to visualize in slicer (not quant)
 
%     niftiwrite(ni,fno,info);

end