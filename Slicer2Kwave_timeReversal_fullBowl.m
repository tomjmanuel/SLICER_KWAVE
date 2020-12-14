% Tom 5/15/20
% script for running time reversal in kwave using simspace
% setup with slicer
%
% This script lets you pick a point in the grid to place source at
% and then returns the phases to drive the multi-element array at
% through skull (if present) to focus at that point

% all xdcr files are generated by makeNiftiTransducer
clear all
%% define inputs
% .mat output holding pressure at sensor (vectorized)
fnout = 'patsesnsor_Emel_fullbowl.mat';

% ct filename 
% CT that has been cropped and co-registered to xdcr mask in slicer
ctfn = 'CT_simspace_20201209.nii.gz';

% xdcr filenames
% point to a nifti generated by makeNiftiTransducer
% this nifti should have been used to resample CT image in slicer
xdfn = 'IGT.nii.gz';

% point to a .mat file with pixel coordinates
elem_locs = load('IGT_pixcoords.mat');

f0 = 0.65E6;               % frequency [MHz]
Amp = 1000;                % amplitude at source.p_mask

% amount to steer (axial dim is third dim)
axialsteer =  0;    %[m] negative steers towards xdcr, positive away
latsteer1   = 0;        %[m] lat steering in first dimension
latsteer2   = 0; %5E-3;        %[m] lat steering in second dimension

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
padsize = [432 432 384];    

%% setup simulation
% put element locations into a grid
% note their vectorized order with ordervar and orderV
Egrid = zeros(dim);
ordervar = zeros(dim);
A = elem_locs.A;
for i=1:length(A)
    Egrid(A(i,1),A(i,2),A(i,3))=1;
    ordervar(A(i,1),A(i,2),A(i,3))=i;
end
orderV = ordervar(:);
orderV(orderV == 0)=[];

% pad arrays
padpre = floor((padsize-dim)/2);
padpost = ceil((padsize-dim)/2);
ct.data = padarray(ct.data,padpre,0,'pre');
ct.data = padarray(ct.data,padpost,0,'post');
xdcr = padarray(xdcr,padpre,0,'pre');
xdcr = padarray(xdcr,padpost,0,'post');
Egrid = padarray(Egrid,padpre,0,'pre');
Egrid = padarray(Egrid,padpost,0,'post');
dim = padsize;
% get medium properties (density and speed of sound for now)
if min(ct.data(:))>=0
    ct.data = double(ct.data)-1024;
    'shifting ct data'
end

medium = getAcousticProperties(ct.data);


%% test in water
% foo = [];
% foo.sound_speed = 1480.*ones(dim);
% foo.density = 997.*ones(dim);
% medium = foo;

%% create kgrid
vox = 1E-3.*ct.info.PixelDimensions;
kgrid = kWaveGrid(dim(1), vox(1), dim(2), vox(2), dim(3), vox(3));
[kgrid.t_array, ~] = makeTime(kgrid, medium.sound_speed);

%% create pressure vector
source.p = createCWSignals(kgrid.t_array, f0, Amp, 0);
foo = round(1/(f0*kgrid.dt)); %samp/cyc
% zero out source after 5 cyc;
source.p(foo*5:end)=0;

%% choose point to focus to
[~,I] = max(xdcr(:));
[l, m, n]= ind2sub(dim,I);
focus_pos = [l, m, n];

% move focus pos according to steering parameters set at top
axsteerpx = round(axialsteer/vox(3));
latsteerpx1 = round(latsteer1/vox(1));
latsteerpx2 = round(latsteer2/vox(2));

focus_pos(1) = focus_pos(1)+latsteerpx1;
focus_pos(2) = focus_pos(2)+latsteerpx2;
focus_pos(3) = focus_pos(3)+axsteerpx;

% check to make sure focus_pos is inside dim
if min(dim-focus_pos)<1||min(focus_pos)<1
    error('focus position is not in the grid')
end


%% setup source
source.p_mask = zeros(dim);
source.p_mask(focus_pos(1),focus_pos(2),focus_pos(3))=1;


%% setup sensor
% old code just used Egrid
% sensor.mask = Egrid;
% here use bowls from xdcr 
sensor.mask = xdcr==1;
sensor.record = {'p'};    

sensor_data = kspaceFirstOrder3DG(kgrid,medium,source,sensor);
%  save output data
save(fnout,'sensor_data','medium','kgrid','focus_pos','orderV','padpre','padpost','ctfn','f0');



