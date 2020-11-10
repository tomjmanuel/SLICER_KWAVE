% Given a nifti file, compute ARFI displacement
% return arfi images along with header info to write nifti in the call
% function

% Input:
%   fni: input filename, path + filename
%       points at a .nii file which is one of 4 .nii output by dcm2niix
%       - dcm2niix can be used directly from gstudy or
%       it can be installed and run from the command line on your pc to
%       convert par rec files into niftis
%       - Given 1 par file, dcm2niix returns 4 seperate nii
%       - These all need to be in same path 
%   fno: output filename

% output:
%   arfi: arfi image
%   info: header info that can be used to write nifti
%   magIm: magnitude image to use for masking in parent function


function [arfi, info, magIm] = nii2MRARFI_two(fni, gradStrength, MEGdur, pos_only)

    % modify fni to point at real and imaginary images
    
    % here, depends on what type of nii files you get from g-study.
    % function "load_data_mix": the nii files downloaded from g-study
    % - FUS-ON:(mix 1: "*real.nii","imaginary.nii"), FUS-OFF(mix 2: "*reala.nii", ,
    % "imaginarya.nii").
    %
    % function "load_data": the nii files downloaded from g-study
    % - ("*real.nii","imaginary.nii"). the data dimension for each nii file
    % should be [nx,ny,nz,2]. 
    % 
    data = [];
    for ii = 1:2:3
        [data(:,:,:,ii:ii+1),info] = load_data_mix(fni{1,(ii+1)/2});
    end
    
    
    % compute ARFI using a crytpic equation that "subtracts dynamics"
    if ~pos_only
        arfi = angle( data(:,:,:,1:4:end) .* conj(data(:,:,:,2:4:end)) .* conj( ...
        data(:,:,:,3:4:end) .* conj(data(:,:,:,3:4:end)) )); % subtraction of dynamics (used timing and averaging)
        % convert from radians to pi
        arfi = arfi ./ (-2*2*pi*42.58*1e6*gradStrength*1e-3*MEGdur*1e-3) * 1e6; % radians -> microns
    end
    if pos_only
        arfi = angle( data(:,:,:,1:4:end) .* conj(data(:,:,:,2:4:end)));
        % convert from radians to pi
        arfi = arfi ./ (-1*2*pi*42.58*1e6*gradStrength*1e-3*MEGdur*1e-3) * 1e6; % radians -> microns
    end
    % get magnitude image for masking later
    magIm = abs(data);
    
end

function [data,Rinfo] = load_data(fni)
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

end

function [data,Rinfo2] = load_data_mix(fni)
    data = [];
    realfn = strcat(fni(1:end-4),'_reala.nii'); %off 
    imagfn = strcat(fni(1:end-4),'_imaginarya.nii'); %off

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

    data(:,:,:,1) = complex(R,Im);

    realfn = strcat(fni(1:end-4),'_real.nii'); %on
    imagfn = strcat(fni(1:end-4),'_imaginary.nii'); %on

    %load real image
    R2 = niftiread(realfn);
    Rinfo2 = niftiinfo(realfn);

    % get R back into its original values
    R2 = squeeze(double(R2));
    R2 = R2.*Rinfo2.MultiplicativeScaling+Rinfo2.AdditiveOffset;

    % load imaginary image
    Im = niftiread(imagfn);
    ImInfo = niftiinfo(imagfn);

    % get R back into its original values
    Im = squeeze(double(Im));
    Im = Im.*ImInfo.MultiplicativeScaling+ImInfo.AdditiveOffset;

    data(:,:,:,2) = complex(R,Im);
end