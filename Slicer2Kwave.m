% Tom 4/22/20
% script for running kwave using simspace setup by slicer
% this script works with H101_makeNii but not for multielement arrays
clear all
%% define inputs
fnout = 'h101_prms';

% ct filename 
% CT that has been cropped and co-registered to xdcr mask in slicer
ctfn = 'CT_h101.nii.gz';

% xdcr filename
% point to a nifti generated by makeNiftiTransducer
% this nifti should have been used to resample CT image in slicer
xdfn = 'H101mask_hard.nii.gz';

f0 = 1E6;                 % frequency [MHz]
Amp = 1000;                    % amplitude at source.p_mask
phase = 0;               % phase of source wave
sensor.record = {'p_rms'};    %set to p_max, or p_rms depending on what you want

%% load data and set final sim size (to include padding for pml & fft efficiency)
ct.data = niftiread(ctfn);
ct.info = niftiinfo(ctfn);
xdcr = niftiread(xdfn);                          
if size(ct.data)~=size(xdcr)
    'WARNING: CT simspace and transducer mask are different sizes'
end

                           
% kwave needs padding around the grid for the "perfectly matched layer"
% also, kwave is much faster if you use numbers with small prime factors
dim = ct.info.ImageSize;
'Your simpace size'
dim
checkFactors(min(dim),max(dim)+100)

%% look at the dimensions of your image and select a final grid size that fits
% your simspace and is in one of the first 3 rows output by checkFactors
padsize = [162 128 128];    

%% setup and run simulation
% pad arrays
padpre = floor((padsize-dim)/2);
padpost = ceil((padsize-dim)/2);
ct.data = padarray(ct.data,padpre,0,'pre');
ct.data = padarray(ct.data,padpost,0,'post');
xdcr = padarray(xdcr,padpre,0,'pre');
xdcr = padarray(xdcr,padpost,0,'post');
dim = padsize;
% get medium properties (density and speed of sound for now)
medium = getAcousticProperties(ct.data);

% create kgrid
vox = 1E-3.*ct.info.PixelDimensions;
kgrid = kWaveGrid(dim(1), vox(1), dim(2), vox(2), dim(3), vox(3));
[kgrid.t_array, dt_SIM] = makeTime(kgrid, medium.sound_speed);

% create pressure vector
source.p = createCWSignals(kgrid.t_array, f0, Amp, phase);
source.p_mask = xdcr==255; %only get spherical cap part of xdcr file
sensor.mask = ones(dim);
%%
sensor_data = kspaceFirstOrder3DG(kgrid,medium,source,sensor);


%% put output into nifti
% grab data
if sensor.record{1} == 'p_rms'
    data = sensor_data.p_rms;
elseif sensor.record{1} == 'p_max'
    data = sensor_data.p_max;
else
    'this script only supports recording p and p_rms'
end
%%
data = reshape(data,dim);

% unpad
data = data(padpre(1)+1:end-padpost(1),padpre(2)+1:end-padpost(2),...
    padpre(3)+1:end-padpost(3));

%%
info = niftiinfo(xdfn);
% if file won't write, you may need to change some properties in info
info.Datatype = 'single';
% write nifti using metadata from ctfile
niftiwrite(data,fnout,info,'compressed',true);



